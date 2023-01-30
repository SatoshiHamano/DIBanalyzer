# -*- coding:utf-8 -*-

__version__ = "2.2"

import sys, os, datetime, time
import glob
import shutil
import numpy, math
from astropy.io import fits
# from pyraf import iraf
import scipy.optimize
import scipy.ndimage
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from regression import confband

from Spec1Dtools import openspecfits
import mysql.connector
from urllib.parse import urlparse
from open_mysql_project import openproject


# Description:
#
#   The shift amount between the first input file and the other input files are searched and shift the other spectrum files to align with the first spectrum file.
#
#
# Usage:
#
#   $ python ccwaveshift_v1_3.py inputlist
#
#   inputlist is the list of the lists, in which the spectrum files with the same echelle orders are listed.
#
#
# Output:
#
#   The shifted spectrum files. The name of the output file is input file name + "s".
#
# Updates:
#
#   ver1.0-1.3 made and updated by Hamano.
#   ver1.0-1.2 are for debug. They cannot be used. ver1.3 is the first version.
#
#   ver2.0 updated by Hamano in Nov. 5th 2017
#       - The resampling method is changed. IRAF scombine is used for the resampling.
#       - PySpecshift is updated. Resampling is done after shifting.
#       - The functions which are not used in this and other scripts are left for the memorial of numerous try-and-error.
#
#   ver2.1 updated by Hamano in Apr. 19th 2018
#       - Python 3
#       - Import updated ver. of spectools (ver1.1).
#       - Not used functions are removed. (They can be seen in old version.)
#       - Interpolation method is replaced from IRAF/scombine with scipy.ndimage.interpolation.shift.
#           With this change, the computational time is shortened.
#
#   ver2.2 updated by Hamano in Sep. 22nd 2019
#       - the version number is removed from the name of the script.
#       - Evaluation function for estimating the shift is changed from absolute to squared.
#
##


def cc_spec_subp_sampling(sp1, sp2, cshift, width, step):
    # This function returns the shift with which the sp2 matches best with sp1.
    # The searching region of shift is from cshift-width to cshift+width.
    # The searching step is set by "step"
    # cshift, width, step are in pix.

    # open the spectra files

    [sp1x, sp1y] = sp1
    [sp2x, sp2y] = sp2

    # normalize them

    normalization_sp2 = numpy.median(sp2y)
    sp1y_norm = sp1y / numpy.median(sp1y)

    # set the parameters

    ec_short = 20  # the region size (in pix) to be removed in the calculation of difference of the spectrum.
    ec_long = 20  # the region size (in pix) to be removed in the calculation of difference of the spectrum.

    # calculate the region of sp1 to be used

    naxis1 = len(sp1x)
    start_sp1 = 0
    end_sp1 = naxis1

    # calculate the shift vector for sp2

    shiftnum = int(width / step * 2 + 1)
    shift = numpy.array(
        [cshift - width + float(i) * step for i in range(shiftnum)])  # calculate the shifts from the input param.

    # shift sp2, resample it and calculate their difference from sp1.

    dify = []
    for i in range(shiftnum):
        sp2y_shifted = scipy.ndimage.interpolation.shift(sp2y, shift[i], order=3)
        sp2y_norm_shifted = sp2y_shifted / normalization_sp2

        tmpydif = (sp1y_norm[start_sp1 + ec_short:end_sp1 - ec_long] - sp2y_norm_shifted[
                                                                       start_sp1 + ec_short:end_sp1 - ec_long]) ** 2.

        dify.append(numpy.average(tmpydif))
    # the shift with which the difference becomes lowest will be returned.

    return numpy.array(shift), numpy.array(dify)


def IRAF_PySpecshift(input, output, shift):
    imfits = "tmpfits_%s-%d%d%d%d.fits" % (
        datetime.date.today(), datetime.datetime.now().hour, datetime.datetime.now().minute,
        datetime.datetime.now().second,
        datetime.datetime.now().microsecond)
    if input.find(".fits") == -1:
        hdulist = fits.open(input + ".fits")
    else:
        hdulist = fits.open(input)
    prihdr = hdulist[0].header
    naxis1 = prihdr["NAXIS1"]
    cdelt1 = prihdr["CDELT1"]
    iraf.scopy(input, imfits)
    iraf.specshift(imfits, shift)
    iraf.scombine(imfits, output, w1=1., dw=cdelt1, nw=naxis1, logfile="null")
    os.remove(imfits)
    hdulist.close()


def atran_resampling(atranmodel, inputspectrum):
    atran_npz = numpy.load(atranmodel)
    wav = atran_npz["wav"]
    flux = atran_npz["flux"]
    id = atran_npz["id"]
    dlam_m = 1.e-2

    spx, spy, crval1, cdelt1, crpix1 = openspecfits(inputspectrum)
    dlam_o = cdelt1

    wav_res = spx
    flux_res = numpy.zeros(spy.shape)

    dlam_ave = (dlam_o + dlam_m) / 2.
    dlam_dif_om = (dlam_o - dlam_m) / 2.
    wav_m_dlamave = wav - dlam_ave
    wav_p_dlamave = wav + dlam_ave

    fulpixbool = [(i - dlam_dif_om <= wav) & (wav < i + dlam_dif_om) for i in spx]

    id_minmax = [[numpy.amin(id[i]) - 1, numpy.amax(id[i]) + 1] for i in fulpixbool]

    for i in range(len(spx)):
        [id_min, id_max] = id_minmax[i]
        dlam_low = (wav_p_dlamave[id_min] - spx[i]) / dlam_m
        dlam_high = (spx[i] - wav_m_dlamave[id_max]) / dlam_m

        flux_res[i] += (flux[id_min] * dlam_low + numpy.sum(flux[id_min + 1:id_max]) + flux[
            id_max] * dlam_high) / (dlam_low + dlam_high + (id_max - id_min - 1))

    return wav_res, flux_res


def waveshift(inputspec, output, m, atranmodel="atran.smo.11513_R28000_0ft.npz"):
    splitn = 10

    wf = open(output, "w")
    wf.write("order\tsplit\twavelength\tshift\tr_edge\tr_min\tr_dif\tr_ratio\n\n")

    # plt.figure()
    # pp2 = PdfPages("waveshift_measure_sample.pdf")

    for j in range(len(m)):
        spx, spy, _, cdelt1, _ = openspecfits(inputspec[j])
        wav_res, flux_res = atran_resampling(atranmodel, inputspec[j])

        naxis1 = len(spx)

        wav_c = []
        shift_c = []
        for i in range(splitn):
            stid = int(naxis1 / splitn * i)
            enid = int(naxis1 / splitn * (i + 1))

            shifts_minus, residue_minus = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                                [spx[stid:enid], spy[stid:enid]], -10.001, 5, 2.5)
            r_minus = numpy.average(residue_minus)
            # plt.scatter(shifts_minus, residue_minus, color="b")

            shifts_plus, residue_plus = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                              [spx[stid:enid], spy[stid:enid]], 10.001, 5, 2.5)
            r_plus = numpy.average(residue_plus)
            # plt.scatter(shifts_plus, residue_plus, color="b")

            r_edge = (r_minus + r_plus) / 2.

            shifts, residue = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], 0.001, 3.0, 0.5)
            # plt.scatter(shifts, residue, color="b")
            cshift = shifts[numpy.argmin(residue)]
            shifts, residue = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], cshift, 0.5, 0.1)
            # plt.scatter(shifts, residue, color="b")

            cshift = shifts[numpy.argmin(residue)]
            shifts, residue = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], cshift, 0.1, 0.02)
            # plt.scatter(shifts, residue, color="b")

            wav_c.append(numpy.average(wav_res[stid:enid]))
            shift_c.append(shifts[numpy.argmin(residue)])

            r_min = numpy.amin(residue)

            # i=0
            wf.write("%d\t%d\t%.3f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n" % (
                m[j], i, wav_c[i], shift_c[i] * cdelt1, r_edge, r_min, (r_edge - r_min), (r_edge - r_min) / r_edge))

            # plt.title("m=%d, %d/10, center=%.3f, delta lambda=%.3e" % (m[j], 5, wav_c[i], shift_c[i] * cdelt1))
            # plt.savefig(pp2, format="pdf")
            # plt.clf()

    wf.close()
    # pp2.close()

def read_waveshift(txtfile):
    rf = open(txtfile, "r")
    rl = rf.readlines()
    rf.close()
    param = numpy.array([[float(rl[i].split()[j]) for i in range(2, len(rl))] for j in range(8)])

    order = param[0]
    split = param[1]
    wavc = param[2]
    shift = param[3]
    r_edge = param[4]
    r_min = param[5]
    r_dif = param[6]
    r_ratio = param[7]

    return order, split, wavc, shift, r_edge, r_min, r_dif, r_ratio


def waveshift_plotter(txtfile, pdffile, fig, colors=["b", "orange"]):
    pp = PdfPages(pdffile)

    order, split, wavc, shift, r_edge, r_min, r_dif, r_ratio = read_waveshift(txtfile)
    m = list(set(order))
    m.sort()

    req = numpy.logical_or(r_edge > 0.5, r_edge - r_min < 1.e-3)
    # req = numpy.array([r_edge[i] > 0.5 for i in range(len(order))])

    clip = numpy.array([True for i in range(len(order))])
    for i in range(5):
        r_ratio_av = numpy.average(r_ratio[numpy.logical_not(req) & clip])
        r_ratio_std = numpy.std(r_ratio[numpy.logical_not(req) & clip])
        clip[(r_ratio - r_ratio_av) < - 2 * r_ratio_std] = False

    req2 = numpy.logical_not(numpy.logical_or(req, numpy.logical_not(clip)))

    fig.suptitle(txtfile.split("/")[-1])

    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    for i in range(len(order)):
        if req[i] or (not clip[i]):
            ax1.scatter(wavc[i], r_ratio[i], s=5, c="gray")
        else:
            ax1.scatter(wavc[i], r_ratio[i], c=colors[int(order[i]) % 2])
    ax1.plot([min(wavc), max(wavc)], [r_ratio_av, r_ratio_av], color="k", linestyle="solid")
    ax1.plot([min(wavc), max(wavc)], [r_ratio_av + 2 * r_ratio_std, r_ratio_av + 2 * r_ratio_std], color="k",
             linestyle="dashed")
    ax1.plot([min(wavc), max(wavc)], [r_ratio_av - 2 * r_ratio_std, r_ratio_av - 2 * r_ratio_std], color="k",
             linestyle="dashed")
    ax1.set_xlabel("Wavelength")
    ax1.set_ylabel("(r_edge - r_min)/r_edge")

    for i in range(len(order)):
        if req[i] or (not clip[i]):
            ax2.scatter(wavc[i], order[i], s=5, c="gray")
        else:
            ax2.scatter(wavc[i], order[i], c=colors[int(order[i]) % 2])
    ax2.set_xlabel("Wavelength")
    ax2.set_ylabel("Order")

    for i in range(len(order)):
        if req[i] or (not clip[i]):
            ax3.scatter(wavc[i], shift[i], s=5, c="gray")
        else:
            ax3.scatter(wavc[i], shift[i], c=colors[int(order[i]) % 2])

    for i in range(len(m)):
        order_req = order[req2]
        wavc_req = wavc[req2]
        shift_req = shift[req2]
        if len(order_req[order_req == m[i]]) > 2:
            a, b = numpy.polyfit(wavc_req[order_req == m[i]], shift_req[order_req == m[i]], 1)
            ax3.plot(wavc_req[order_req == m[i]], wavc_req[order_req == m[i]] * a + b, color="k")

    ax3.set_xlabel("Wavelength")
    ax3.set_ylabel("Shift")

    for i in range(len(order)):
        if req[i] or (not clip[i]):
            pass
        else:
            plt.scatter(wavc[i], r_dif[i], c=colors[int(order[i]) % 2])
    ax4.set_xlabel("Wavelength")
    ax4.set_ylabel("r_edge - r_min")
    plt.savefig(pp, format="pdf")
    plt.clf()
    pp.close()

def waveshift_dif(targetfile, reffile, fig, pp, clipsig=2.5, clipite=5):
    ax1 = fig.add_subplot(1,1,1)

    nite = 5
    lowsig = 2
    clipsig = 5

    order1, split1, wavc1, shift1, r_edge1, r_min1, r_dif1, r_ratio1 = read_waveshift(targetfile)
    order2, split2, wavc2, shift2, r_edge2, r_min2, r_dif2, r_ratio2 = read_waveshift(reffile)

    m = list(set(order1))
    m.sort()
    wavc_m = {}
    for i in m:
        wavc_m[i] = numpy.average(wavc1[order1==i])

    req1 = numpy.logical_or(r_edge1 > 0.5, r_edge1 - r_min1 < 1.e-3)
    # req1 = numpy.array([r_edge1[i] > 0.5 for i in range(len(order1))])

    clip1 = numpy.array([True for i in range(len(order1))])
    for i in range(nite):
        r_ratio_av1 = numpy.average(r_ratio1[numpy.logical_not(req1) & clip1])
        r_ratio_std1 = numpy.std(r_ratio1[numpy.logical_not(req1) & clip1])
        clip1[(r_ratio1 - r_ratio_av1) < - lowsig * r_ratio_std1] = False

    req_target = numpy.logical_not(numpy.logical_or(req1, numpy.logical_not(clip1)))

    req2 = numpy.array([r_edge2[i] > 0.5 for i in range(len(order2))])

    clip2 = numpy.array([True for i in range(len(order2))])
    for i in range(nite):
        r_ratio_av2 = numpy.average(r_ratio2[numpy.logical_not(req2) & clip2])
        r_ratio_std2 = numpy.std(r_ratio2[numpy.logical_not(req2) & clip2])
        clip2[(r_ratio2 - r_ratio_av2) < - lowsig * r_ratio_std2] = False

    req_ref = numpy.logical_not(numpy.logical_or(req2, numpy.logical_not(clip2)))

    req_com1 = numpy.logical_and(req_target, req_ref)
    stdshift = numpy.std(shift1[req_com1] - shift2[req_com1])
    aveshift = numpy.average(shift1[req_com1] - shift2[req_com1])
    clip3 = numpy.array([True for i in range(len(order2))])
    clip3[numpy.absolute(shift1 - shift2 - aveshift) > clipsig * stdshift] = False
    for i in range(nite):
        req_com = numpy.logical_and(req_com1, clip3)
        stdshift = numpy.std(shift1[req_com] - shift2[req_com])
        aveshift = numpy.average(shift1[req_com] - shift2[req_com])
        clip3[numpy.absolute(shift1 - shift2 - aveshift) > clipsig * stdshift] = False
    req_com = numpy.logical_and(req_com1, clip3)

    a, b = numpy.polyfit(wavc1[req_com], shift1[req_com] - shift2[req_com], 1)

    shift_m = {}
    wavc_arr = []
    shift_arr = []
    for i in m:
        shift_m[i] = a * wavc_m[i] + b
        wavc_arr.append(wavc_m[i])
        shift_arr.append(shift_m[i])
    wavc_arr = numpy.array(wavc_arr)
    shift_arr = numpy.array(shift_arr)

    lcb, ucb, wav = confband(wavc1[req_com], shift1[req_com] - shift2[req_com], a, b, wavc_arr)
    shift_err = (ucb - lcb) / 2.
    shift_err_m = {}
    for i in m:
        shift_err_m[i] = shift_err[m.index(i)]

    ax1.plot(wavc1[req_com], wavc1[req_com] * a + b, "k", lw=0.7, label=r"$\Delta s = (%.1e) * \lambda + (%.1e)$" % (a,b))
    ax1.plot(wav, lcb, "k--", lw=0.7, label=r"Conf. band ($\alpha = 0.95$)")
    ax1.plot(wav, ucb, "k--", lw=0.7)
    ax1.scatter(wavc1[req_com], shift1[req_com] - shift2[req_com], s=4)
    ax1.legend()
    ax1.set_xlabel(r"Wavelength ($\AA$)")
    ax1.set_ylabel(r"$\Delta s =$ s(target) - s(ref) ($\AA$)")
    ax1.set_ylim(min(-0.5, min(shift1[req_com] - shift2[req_com])), max(0.5, max(shift1[req_com] - shift2[req_com])))
    ax1.set_title("%s \n %s" % (targetfile.split("/")[-1], reffile.split("/")[-1]), fontsize=7)

    plt.savefig(pp, format="pdf")
    plt.clf()

    return m, shift_m, shift_err_m


if __name__ == "__main__":
    conn, cur = openproject()

    cur.execute(
        "select pipelineID, FrameNum, totalSNR, mode, obsdate, path, pipelinever from datareduction;")# where obsdate between '2016-09-01 00:00:00' and '2017-01-01 00:00:00';")
    rows = cur.fetchall()
    ppid = [i[0] for i in rows]
    fnum = [int(i[1]) for i in rows]
    snr = [float(i[2]) for i in rows]
    mode = [i[3] for i in rows]
    obsdate = [i[4] for i in rows]
    path = [i[5] for i in rows]
    pver = [i[6] for i in rows]

    fig = plt.figure(figsize=(16, 12))

    for i in range(len(ppid)):
        print(ppid[i])
        wsdir = path[i] + "waveshift_measure/"
        if not os.path.exists(wsdir):
            frames = ["NO%d" % (f + 1) for f in range(fnum[i])] + ["sum"]

            os.makedirs(wsdir)
            for fr in frames:
                os.makedirs(wsdir + fr)

            cur.execute("select wsave, wsstd, wsnum from reducedframe where pipelineID='%s';" % ppid[i])
            rows = cur.fetchall()
            wsave = [float(j[0]) for j in rows]
            wsstd = [float(j[1]) for j in rows]
            wsnum = [int(j[2]) for j in rows]

            for fr in frames:
                fsrdir = glob.glob(path[i] + "*_%s/VAC_norm/fsr*" % fr)
                fsrdir.sort()
                for fsr in fsrdir:
                    os.makedirs(wsdir + fr + "/" + fsr.split("/")[-1])
                    spfiles = glob.glob(fsr + "/*_norm.fits")
                    spfiles.sort()
                    m = [int(spf.rstrip("_VAC_norm.fits").rstrip(fsr.split("/")[-1]).rstrip("_").split("_m")[-1]) for spf in spfiles]
                    output = wsdir + fr + "/" + fsr.split("/")[-1] + "/" + ppid[i] + "_" + fsr.split("/")[
                        -1] + "_VAC_norm_%s.txt" % fr
                    outputpdf = output.replace(".txt", ".pdf")
                    waveshift(spfiles, output, m)
                    waveshift_plotter(output, outputpdf, fig)

    conn.close()
