# -*- coding:utf-8 -*-

__version__ = "2.2"

import sys, os, datetime, time
import numpy, math
from astropy.io import fits
from pyraf import iraf
import scipy.optimize
import scipy.ndimage
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time

from Spec1Dtools import openspecfits

iraf.noao()
iraf.onedspec()


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


def min_shift_search(shift, dify):
    dec_s = numpy.array([math.fabs(math.modf(shift[i])[0]) for i in range(len(shift))])
    dify_corr = dify - p_smooth[1] * (2. * (dec_s - 0.5) ** 2. + 0.5) ** 0.5

    return shift[numpy.argmin(dify_corr)]


def PySpecshift(input, output, shift):
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


def waveshift_oneorder(files):
    # the wavelength shift vector
    # the shift of 1st file from 1st file is set as 0.0.
    wshift_vec = [0.0]
    sp1x, sp1y, crval1, cdelt1, crpix1 = openspecfits(files[0])
    for i in range(1, len(files)):
        print("Calculating spectra shift between %s and %s" % (files[0], files[i]))

        # the first guess. center is 0. step is large (0.5pix). width is large (2pix)
        shift_1, dify_1 = cc_spec_subp_sampling(files[0], files[i], 0.001, 2.0, 0.5)
        cshift_p = shift_1[numpy.argmin(dify_1)]

        # the second guess. center is from first guess. the width is set as the step in first guess.
        shift_2, dify_2 = cc_spec_subp_sampling(files[0], files[i], cshift_p, 0.4, 0.1)
        cshift_sp = shift_2[numpy.argmin(dify_2)]

        # the last guess. center is from 2nd guess. the width is set as the step in second guess.
        shift_3, dify_3 = cc_spec_subp_sampling(files[0], files[i], cshift_sp, 0.07, 0.01)
        cshift_ssp = shift_3[numpy.argmin(dify_3)]

        wshift_vec.append(cshift_ssp * cdelt1)

    return wshift_vec


def waveshift_multiorder(files):
    # calculate the wavelength shift for each orders and make matrix of wavelength shift
    # ane axis is input files. the other axis is echelle orders.
    wshift_matrix = []
    for i in range(len(files)):
        wshift_matrix.append(waveshift_oneorder(files[i]))

    # calculate the median of wavelength shift.
    shift_median = []
    for j in range(len(files[0])):
        tmplist = []
        for i in range(len(files)):
            tmplist.append(wshift_matrix[i][j])
        shift_median.append(numpy.median(tmplist))

    return shift_median


def waveshiftClip(wshift_matrix, objnum, aplength, sigma_1st=1., sigma=2., iterate=5, stdthres=0.1):
    shift_average = []
    shift_calcnum = []
    shift_stddev = []
    for i in range(objnum):
        shiftvec = []
        shiftclip = numpy.array([0 for jj in range(aplength)])
        for j in range(aplength):
            shiftvec.append(wshift_matrix[j][i])
        shiftvecarray = numpy.array(shiftvec)

        for k in range(iterate):
            shift_average_ite = numpy.average(shiftvecarray[shiftclip == 0])
            shift_std_ite = numpy.std(shiftvecarray[shiftclip == 0])
            if k == 0:
                shiftclip[numpy.absolute(shiftvecarray - shift_average_ite) / shift_std_ite > sigma_1st] += 1
            elif k != iterate - 1:
                shiftclip[numpy.absolute(shiftvecarray - shift_average_ite) / shift_std_ite > sigma] += 1

        shift_average.append(shift_average_ite)
        shift_stddev.append(shift_std_ite)
        shift_calcnum.append(numpy.sum(shiftclip == 0))

    shift_stddev_arr = numpy.array(shift_stddev)
    shift_calcnum_arr = numpy.array(shift_calcnum)
    wshift_flag = numpy.zeros(objnum)
    wshift_flag[shift_stddev_arr > stdthres] += 1
    wshift_flag[shift_calcnum_arr < aplength / 2] += 1

    for i in range(objnum):
        if wshift_flag[i] != 0:
            shift_average[i] = 0.
            shift_calcnum[i] = 0

    return shift_average, shift_calcnum, shift_stddev


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


if __name__ == "__main__":
    plt.figure(figsize=(30,6))
    colors = ["k", "b"]
    splitn = 10

    pp = PdfPages("resample_shift_residue2.pdf")

    for m in range(42, 62, 1):
        inputspec = "HD183143_sum/VAC_norm/fsr1.30/HD183143_sum_m%d_fsr1.30_VAC_norm.fits" % m

        spx, spy, _, cdelt1, _ = openspecfits(inputspec)
        wav_res, flux_res = atran_resampling("atran.smo.11513_R28000_0ft.npz", inputspec)

        naxis1 = len(spx)

        wav_c = []
        shift_c = []
        for i in range(splitn):
            stid = int(naxis1 / splitn * i)
            enid = int(naxis1 / splitn * (i + 1))
            shifts_minus, residue_minus = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], -10.001, 5, 2.5)
            # plt.scatter(shifts_minus, residue_minus)
            r_minus = numpy.average(residue_minus)
            # plt.ylim(rmin - (rmax - rmin) * 0.1, rmax + (rmax - rmin) * 0.1)

            shifts_plus, residue_plus = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], 10.001, 5, 2.5)
            # plt.scatter(shifts_minus, residue_minus)
            r_plus = numpy.average(residue_plus)

            r_edge = (r_minus + r_plus)/2.

            shifts, residue = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], 0.001, 3.0, 0.5)
            # plt.scatter(shifts, residue)
            cshift = shifts[numpy.argmin(residue)]
            shifts, residue = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], cshift, 0.5, 0.1)
            # plt.scatter(shifts, residue)
            cshift = shifts[numpy.argmin(residue)]
            shifts, residue = cc_spec_subp_sampling([wav_res[stid:enid], flux_res[stid:enid]],
                                                    [spx[stid:enid], spy[stid:enid]], cshift, 0.1, 0.02)
            # plt.scatter(shifts, residue)
            # plt.savefig(pp, format="pdf")
            # plt.clf()

            # plt.scatter(shifts, residue, s=5)
            # plt.scatter(shifts[numpy.argmin(residue)], residue[numpy.argmin(residue)], s=5)

            wav_c.append(numpy.average(wav_res[stid:enid]))
            shift_c.append(shifts[numpy.argmin(residue)])

            r_min = numpy.amin(residue)

            print("%d\t%d\t%.3f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e" % (m, i, wav_c[i], shift_c[i] * cdelt1, r_edge, r_min, (r_edge-r_min), (r_edge - r_min) / r_edge))

            plt.step(wav_res[stid:enid], flux_res[stid:enid], where="mid", color=colors[0], label="model")
            plt.step(spx[stid:enid], spy[stid:enid] - 0.1, where="mid", color=colors[1], label="orig")
            plt.step(spx[stid:enid] + shift_c[i] * cdelt1, spy[stid:enid] + 0.1, where="mid", color="r", label="shifted")
            #
            # # plt.scatter(wav_c[i], shift_c[i], color=colors[i%2])
            #
            plt.title("m=%d (%d/%d): shift=%.3f" % (m, (i+1), splitn, shift_c[i] * cdelt1))
            plt.legend()
            plt.grid()
            plt.savefig(pp, format="pdf")
            plt.clf()


    pp.close()
        #
        # rf = open("atran.smo.1412.dat", "r")
        # rl = rf.readlines()
        # rf.close()
        #
        # wav = numpy.array([float(i.split()[1]) * 1.e+4 for i in rl])
        # flux = numpy.array([float(i.split()[2]) for i in rl])
        #
        # spx, spy, _, _, _ = openspecfits("HD183143_sum/VAC_norm/fsr1.30/HD183143_sum_m44_fsr1.30_VAC_norm.fits")
        #
        # plt.xlim(12600 - 5, 12600 + 15)
        #
        # plt.step(wav_res, flux_res, where="mid", label="resampled model")
        # plt.step(wav, flux - 0.1, where="mid", label="model")
        # plt.step(spx, spy + 0.1 , where="mid", label="input observed spectrum")
        #
        # plt.legend()
        #
        # plt.savefig("resample.png")

        # filename = sys.argv[1:]
        #
        # rfile = open(filename[0], "r")
        # rlines = rfile.readlines()
        # listfiles = [rlines[i].split()[0] for i in range(len(rlines))]
        # rfile.close()
        #
        # fileend = "s"
        #
        # files = []
        #
        # for i in range(len(listfiles)):
        #     rflist = open(listfiles[i], "r")
        #     rflines = rflist.readlines()
        #     rflist.close()
        #     files.append([rflines[j].split()[0] for j in range(len(rflines))])
        #
        # shift_median = waveshift_multiorder(files)
        #
        # for i in range(len(rlines)):
        #     for j in range(len(files[i])):
        #         PySpecshift(files[i][j], files[i][j].rstrip("fits").rstrip(".") + fileend, shift_median[j])
