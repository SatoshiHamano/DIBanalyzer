# -*- coding:utf-8 -*-

from waveshift_measure import waveshift_dif
import sys
import mysql.connector
from urllib.parse import urlparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.interpolate
import scipy.ndimage
import glob
from Spec1Dtools import openspecfits, savespecfits
from astropy.io import fits
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler
from specutils import Spectrum1D
from astropy import units as u
import numpy
import os
import shutil
# from pyraf import iraf
from spectra_plotter import spectrum_plot
import argparse
from open_mysql_project import openproject

# iraf.noao()
# iraf.onedspec()


def alternativequestion(question, anss, defans):
    flagans = 0
    while flagans == 0:
        flag = input(question)
        if flag in anss:
            flagans += 1
        else:
            print("Answers: ", anss)

    if flag != "":
        return flag
    else:
        return defans


def scale_measure(sp1, sp2, relshift, shifterr, w_low, w_high, fig, pp):
    # open the spectra files

    trans_thres = 0.8

    sp1x, sp1y, a1, b1, c1 = openspecfits(sp1)
    sp2x, sp2y, a2, b2, c2 = openspecfits(sp2)
    n1 = len(sp1x)
    n2 = len(sp2x)
    if n1 > n2:
        sp1x = sp1x[0:n2]
        sp1y = sp1y[0:n2]
    elif n2 > n1:
        sp2x = sp2x[0:n2]
        sp2y = sp2y[0:n2]

    if not numpy.all(numpy.absolute(sp1x - sp2x) < 1e-3):
        print("The spectrum wavelength parameters are not the same.")
        print(sp1, a1, b1, c1)
        print(sp2, a2, b2, c2)
        sys.exit()
    sp1_hdulist = fits.open(sp1)
    sp2_hdulist = fits.open(sp2)
    sp1_hdr = sp1_hdulist[0].header
    sp2_hdr = sp2_hdulist[0].header
    sp1_airmass = float(sp1_hdr["AIRMASS"])
    sp2_airmass = float(sp2_hdr["AIRMASS"])

    # cut sp1, sp2

    sp1x_cut = sp1x[(w_low < sp1x) & (sp1x < w_high)]
    sp2x_cut = sp2x[(w_low < sp2x) & (sp2x < w_high)]
    sp1y_cut = sp1y[(w_low < sp1x) & (sp1x < w_high)]
    sp2y_cut = sp2y[(w_low < sp2x) & (sp2x < w_high)]

    # shift sp2

    sp2y_shifted = scipy.ndimage.interpolation.shift(sp2y_cut, - relshift / b1, order=3)

    # set the parameters

    div = 100  # (1.-2./div)*100 percent of the spectrum will be used
    edgecut = 20  # the region size (in pix) to be removed in the calculation of difference of the spectrum.

    # calculate the region of sp1 to be used

    naxis1 = len(sp1x_cut)
    start_sp1 = int(naxis1 / div)
    end_sp1 = int(naxis1 / div * (div - 1))

    ax1 = fig.add_subplot(1, 1, 1)

    # ax1.step(sp1x, sp1y_norm, where="mid")
    # f = scipy.interpolate.interp1d(sp2x, sp2y_norm, kind="cubic")
    # ax1.step(sp1x, f(sp1x)-0.1, where="mid")
    ax1.step(sp1x, sp1y, where="mid", color="0.5")
    ax1.step(sp2x, sp2y - 1.0, where="mid", color="0.5")
    ax1.step(sp1x_cut, sp1y_cut, where="mid", color="r", label="Target")
    ax1.step(sp2x_cut, sp2y_shifted - 1.0, where="mid", color="b", label="Telluric")
    ax1.plot(sp2x_cut, sp1y_cut / sp2y_shifted, color="yellow", label="Target / Telluric")
    ax1.legend()

    plt.savefig(pp, format="pdf")
    plt.clf()

    scalevector = numpy.log(numpy.fabs(sp1y_cut[start_sp1 + edgecut:end_sp1 - edgecut])) / numpy.log(
        numpy.fabs(sp2y_shifted[start_sp1 + edgecut:end_sp1 - edgecut])) / sp1_airmass * sp2_airmass
    transmittance = numpy.fabs(sp2y_shifted[start_sp1 + edgecut:end_sp1 - edgecut])

    scale_av_1 = numpy.average(scalevector[transmittance < trans_thres])
    scale_low_1 = numpy.average(scalevector[transmittance < 0.5])
    scale_high_1 = numpy.average(scalevector[numpy.where((0.5 < transmittance) & (transmittance < trans_thres))])
    scale_std_1 = numpy.std(scalevector[transmittance < trans_thres])

    sig = 5.

    clip_array = numpy.where(
        (numpy.fabs(scale_av_1 - scalevector) < scale_std_1 * sig) & (transmittance < trans_thres))

    scale_av = numpy.average(scalevector[clip_array])
    scale_std = numpy.std(scalevector[clip_array])
    T_arr = numpy.array([transmittance[clip_array], numpy.ones(len(transmittance[clip_array]))])
    T_arr = T_arr.T
    scale_a, scale_b = numpy.linalg.lstsq(T_arr, scalevector[clip_array])[0]

    ax1 = fig.add_subplot(1, 1, 1)

    ax1.scatter(transmittance, scalevector, s=7, c="dimgray", edgecolors="face")
    ax1.scatter(transmittance[numpy.where(
        (numpy.fabs(scale_av_1 - scalevector) < scale_std_1 * 2.) & (transmittance < trans_thres))], scalevector[
                    numpy.where(
                        (numpy.fabs(scale_av_1 - scalevector) < scale_std_1 * 2.) & (transmittance < trans_thres))],
                s=9, c="blue", edgecolors="face")
    ax1.plot([-0.2, 1.1], [scale_av, scale_av], "k--")
    ax1.text(0.05, 1.9,
             "OBJECT: " + sp1.split("/")[-1].rstrip("fits").rstrip(".") + " (airmass=%.2f)" % sp1_airmass,
             fontsize=10)
    ax1.text(0.05, 1.8,
             "TELLURIC: " + sp2.split("/")[-1].rstrip("fits").rstrip(".") + " (airmass=%.2f)" % sp2_airmass,
             fontsize=10)
    # ax1.text(0.05, 0.3, r"RANGE: %d - %d $\AA$" % (w_low, w_high), fontsize=10)
    ax1.text(0.05, 0.2, r"SHIFT: %.3f +/- %.3f$\AA$" % (relshift, shifterr), fontsize=10)
    ax1.text(0.05, 0.1, "SCALE: %.3f +/- %.3f" % (scale_av, scale_std), fontsize=10)
    ax1.set_xlabel("Transmittance")
    ax1.set_ylabel("Scale")
    ax1.set_ylim(0, 2.)
    ax1.set_xlim(0.0, 1.1)

    plt.savefig(pp, format="pdf")
    plt.clf()

    return scale_av, scale_std, scale_a, scale_b


def telluricAstropy(inputsp, refsp, outputsp, shift=0.0, scale=1.0, thres=0.001):
    inputspec = Spectrum1D.read(inputsp)
    refspec = Spectrum1D.read(refsp)
    inputfits = fits.open(inputsp)
    reffits = fits.open(refsp)
    aminput = float(inputfits[0].header["AIRMASS"])
    amref = float(reffits[0].header["AIRMASS"])

    dlamref = refspec.wavelength.value - numpy.roll(refspec.wavelength.value, 1)
    lamperpix = numpy.median(dlamref[1:])
    shiftedspec = Spectrum1D(spectral_axis=refspec.wavelength + lamperpix * shift * u.AA, flux=refspec.flux)
    spline = SplineInterpolatedResampler(bin_edges='zero_fill')
    resampledspec = spline(shiftedspec, inputspec.wavelength)
    resampledspec.flux[resampledspec.flux < thres] = thres

    scaledspec = Spectrum1D(spectral_axis=resampledspec.wavelength, flux=resampledspec.flux ** (aminput/amref*scale))

    outputspec = inputspec / scaledspec
    savespecfits(outputspec, inputsp, outputsp)


def GetListAdvancedFits(telpath, fsr):
    dirlist = glob.glob("%s*" % telpath)
    advdircand = []
    advdir = "NA"
    print("Dir list: ")
    for i in dirlist:
        if os.path.isdir(i) and i.split(telpath)[1].find("_v") != -1:
            print(i.split(telpath)[1])
            advdircand.append(i)

    if len(advdircand) > 1:
        for i in advdircand:
            if os.path.isdir(i) and i.split(telpath)[1].find("_v") != -1:
                ans = alternativequestion("%s is the directory you want? " % i.split(telpath)[1], ["yes", "no"], "yes")
                if ans == "yes":
                    advdir = i
                    break
    else:
        print("%s was selected." % advdircand[0].split(telpath)[1])
        advdir = advdircand[0]

    if advdir == "NA":
        print("No dir found.")
        sys.exit()
    telfiles = glob.glob("%s/telluric/FITS/%s/telluric_m*.fits" % (advdir, fsr))
    if telfiles == []: telfiles = glob.glob("%s/telluric/FITS/%s/telluric_m*.fits" % (advdir, "cut?"))

    return telfiles

if __name__ == "__main__":
    # option setting

    parser = argparse.ArgumentParser()
    parser.add_argument("targetid", type=str, help="targetid")
    parser.add_argument("refid", type=str, help="refid")
    parser.add_argument("-a", "--advanced", action="store_true", help="use processed standard spectra")
    parser.add_argument("-m", "--manual", action="store_true", help="manual telluric")
    parser.add_argument("-f", "--fsr", action="store_true", help="Change to FSR1.05")

    args = parser.parse_args()

    targetid = args.targetid
    refid = args.refid
    advancedflag = args.advanced
    manualflag = args.manual
    fsrflag = args.fsr
    autoflag = not manualflag

    curdir = os.getcwd()

    print(advancedflag, manualflag, autoflag)

    advint = 1 if advancedflag else 0
    autoint = 1 if autoflag else 0

    conn, cur = openproject()
    fig = plt.figure(figsize=(20, 6))

    wsstr1 = "waveshift_measure/"
    fsr = "fsr1.05" if fsrflag else "fsr1.30"
    vacorair = "VAC"
    year = int(targetid[0:4])

    if year < 2017:
        orders_sc = [43, 44, 47, 51, 58]
        low_w = {43: 12970., 44: 12630., 47: 11800., 51: 11050., 58: 9680.}
        high_w = {43: 13150., 44: 12770., 47: 11940., 51: 11120., 58: 9780.}
        orders_o2 = [44]
        orders_h2o = [43, 47, 51, 58]
    else:
        orders_sc = [43, 44, 48, 50, 58]
        low_w = {43: 12970., 44: 12630., 48: 11540., 50: 11080., 58: 9680.}
        high_w = {43: 13150., 44: 12770., 48: 11620., 50: 11250., 58: 9780.}
        orders_o2 = [44]
        orders_h2o = [43, 48, 50, 58]


    cur.execute(
        "select autoflag, advanced, telluricNumber from telluriccorrection where pipelineIDobj='%s' and pipelineIDtel='%s';" % (
            targetid, refid))
    rows = cur.fetchall()
    if rows == []:
        telluricNumber = 1
    else:
        autoflag_prev = [i[0] for i in rows]
        advanced_prev = [i[1] for i in rows]
        telluricNumber_prev = [i[2] for i in rows]
        for i in range(len(rows)):
            if autoflag_prev[i] == autoint and advanced_prev[i] == advint:
                ans = alternativequestion("You already did the same analysis. Proceed anyway?: ", ["yes", "no"], "no")
                if ans == "no":
                    print("Bye")
                    sys.exit()
                else:
                    break

        telluricNumber = max(telluricNumber_prev) + 1

    telluricID = targetid + "--" + refid[11:] + "-T%d" % telluricNumber
    teldir = telluricID

    cur.execute(
        "select path, pipelinever, mode from datareduction where pipelineID = '%s';" % targetid)
    rows = cur.fetchall()
    targetpath = rows[0][0]
    targetpver = rows[0][1]
    targetmode = rows[0][2]
    cur.execute(
        "select z.slit from datareduction as x join reducedframe as y using(pipelineID) join observation as z on y.objectframe=z.frame where x.pipelineID='%s'" % targetid)
    rows = cur.fetchall()
    targetslit = list(set([i[0] for i in rows]))
    if len(targetslit) != 1:
        print("Multiple slit size data was included in the target data: ", targetslit)
        sys.exit()
    telluricpath = targetpath + telluricID + "/"
    targetframeglob = glob.glob(targetpath + wsstr1 + "*")
    if len(targetframeglob) == 0:
        print("waveshift_measure have not been done.")
        sys.exit()

    targetframeglob.sort()
    targetframe = [i.split("/")[-1] for i in targetframeglob]
    target_wspath = ["%s%s%s/%s/%s_%s_%s_norm_%s.txt" % (targetpath, wsstr1, i, fsr, targetid, fsr, vacorair, i) for i
                     in
                     targetframe]
    targetnormdir = [glob.glob("%s*_%s/%s_norm/%s/" % (targetpath, i, vacorair, fsr)) for i in targetframe]
    targetfluxdir = [glob.glob("%s*_%s/%s_flux/%s/" % (targetpath, i, vacorair, fsr)) for i in targetframe]
    targetspec = [glob.glob("%s*_%s/%s_norm/%s/*_norm.fits" % (targetpath, i, vacorair, fsr)) for i in targetframe]
    targetspecflux = [glob.glob("%s*_%s/%s_flux/%s/*%s.fits" % (targetpath, i, vacorair, fsr, vacorair)) for i in
                      targetframe]
    for i in targetspec: i.sort()
    for i in targetspecflux: i.sort()

    cur.execute(
        "select path, pipelinever, mode from datareduction where pipelineID = '%s';" % refid)
    rows = cur.fetchall()
    refpath = rows[0][0]
    refpver = rows[0][1]
    refmode = rows[0][2]
    cur.execute(
        "select z.slit from datareduction as x join reducedframe as y using(pipelineID) join observation as z on y.objectframe=z.frame where x.pipelineID='%s'" % targetid)
    rows = cur.fetchall()
    refslit = list(set([i[0] for i in rows]))
    if len(refslit) != 1:
        print("Multiple slit size data was included in the reference data: ", refslit)
        sys.exit()

    if targetslit != refslit:
        print("Slit size are not the same.")
        sys.exit()

    # if (targetpver, targetmode) != (refpver, refmode):
    #     print(targetpver, targetmode)
    #     print(refpver, refmode)
    #     sys.exit()

    advtext = "true" if advancedflag else "false"
    autotext = "true" if autoflag else "false"
    valuetext = "'" + telluricID + "','" + targetid + "','" + refid + "'," + autotext + "," + advtext + ",'" + targetpver + "','" + targetmode + "'," + str(
        telluricNumber) + ",'" + telluricpath + "'"

    ref_wspath = "%s%ssum/%s/%s_%s_%s_norm_sum.txt" % (refpath, wsstr1, fsr, refid, fsr, vacorair)

    if advancedflag:
        refspec = GetListAdvancedFits(refpath, "fsr?.??")
        refspecflux = [i.replace("telluric_m", "telluric_flux_m") for i in refspec]
    else:
        refspec = glob.glob("%s*_sum/%s_norm/%s/*norm.fits" % (refpath, vacorair, fsr))
        refspecflux = glob.glob("%s*_sum/%s_flux/%s/*%s.fits" % (refpath, vacorair, fsr, vacorair))

    refspec.sort()
    refspecflux.sort()

    normsp = glob.glob("%s*_sum/%s_norm/%s/*%s_norm.fits" % (refpath, vacorair, fsr, vacorair))
    fluxsp = glob.glob("%s*_sum/%s_flux/%s/*%s.fits" % (refpath, vacorair, fsr, vacorair))
    normsp.sort()
    fluxsp.sort()

    fluxdir = glob.glob("%s*_sum/%s_flux" % (refpath, vacorair))[0]
    contdir = fluxdir.replace("_flux", "_cont") + "/%s" % fsr
    refspeccont = [contdir + ("/%s" % (i.split("/")[-1].rstrip("fits").rstrip("."))) + "_cont.fits" for i in fluxsp]
    if not os.path.exists(contdir):
        os.makedirs(contdir)
        for (i,j,k) in zip(fluxsp,normsp,refspeccont):
            print(fluxsp)
            fluxspec = Spectrum1D.read(fluxsp)
            normspec = Spectrum1D.read(normsp)
            contspec = fluxspec / normspec
            savespecfits(contspec, fluxsp, refspeccont)

    if advancedflag:
        for i in range(len(refspec)):
            if not os.path.exists(refspecflux[i]):
                normspec = Spectrum1D.read(refspec[i])
                contspec = Spectrum1D.read(refspeccont[i])
                fluxspec = normspec * contspec
                savespecfits(fluxspec, refspec[i], refspecflux[i])

    if autoflag:
        print("auto")
        cur.execute(
            "INSERT IGNORE INTO telluriccorrection (telluricID,pipelineIDobj,pipelineIDtel,autoflag,advanced,pipelinever,mode,telluricNumber,telluricPath) VALUES (%s);" % valuetext)
        for i in range(len(targetframe)):
            teldir_frame = telluricpath + targetframe[i] + "/"
            os.makedirs(teldir_frame)

            pp = PdfPages(teldir_frame + telluricID + "_" + targetframe[i] + ".pdf")

            m, shift_m, shift_err_m = waveshift_dif(target_wspath[i], ref_wspath, fig, pp)
            scale_o2_arr, scale_err_o2_arr, scale_h2o_arr, scale_err_h2o_arr = [], [], [], []

            m_common = []
            targetspec_m = {}
            refspec_m = {}
            targetspecflux_m = {}
            refspecflux_m = {}
            for j in m:
                tid = -1
                rid = -1
                for l in range(len(targetspec[i])):
                    if targetspec[i][l].split("/")[-1].find("m%d" % j) != -1:
                        tid = l
                        break
                for k in range(len(refspec)):
                    if refspec[k].split("/")[-1].find("m%d" % j) != -1:
                        rid = k
                        break
                if tid != -1 and rid != -1:
                    m_common.append(j)
                    targetspec_m[j] = targetspec[i][tid]
                    refspec_m[j] = refspec[rid]
                    targetspecflux_m[j] = targetspecflux[i][tid]
                    refspecflux_m[j] = refspecflux[rid]

            for mm in m_common:
                if mm in orders_sc:
                    scale_av, scale_std, scale_a, scale_b = scale_measure(targetspec_m[mm], refspec_m[mm], shift_m[mm],
                                                                          shift_err_m[mm], low_w[mm], high_w[mm],
                                                                          fig, pp)
                    if mm in orders_o2:
                        scale_o2_arr.append(scale_av)
                        scale_err_o2_arr.append(scale_std)
                    else:
                        scale_h2o_arr.append(scale_av)
                        scale_err_h2o_arr.append(scale_std)

            scale_o2_arr = numpy.array(scale_o2_arr)
            scale_err_o2_arr = numpy.array(scale_err_o2_arr)
            scale_h2o_arr = numpy.array(scale_h2o_arr)
            scale_err_h2o_arr = numpy.array(scale_err_h2o_arr)
            scale_o2 = numpy.sum((1. / scale_err_o2_arr) ** 2. * scale_o2_arr) / numpy.sum(
                (1. / scale_err_o2_arr) ** 2.)
            scale_h2o = numpy.sum((1. / scale_err_h2o_arr) ** 2 * scale_h2o_arr) / numpy.sum(
                (1. / scale_err_h2o_arr) ** 2)
            scale_err_o2 = (1. / numpy.sum(1. / scale_err_o2_arr ** 2))**0.5
            scale_err_h2o = (1. / numpy.sum(1. / scale_err_h2o_arr ** 2))**0.5

            # os.chdir(teldir_frame)
            for j in m_common:
                scale_cur = scale_o2 if j in orders_o2 else scale_h2o
                scale_err_cur = scale_err_o2 if j in orders_o2 else scale_err_h2o
                rspec_cur = refspecflux_m[j]

                spx_t, spy_t, _, dw, _ = openspecfits(targetspecflux_m[j])
                spx_r, spy_r, _, _, _ = openspecfits(rspec_cur)

                telfile = teldir_frame + targetspecflux_m[j].split("/")[-1].rstrip("fits").rstrip(".") + "_tel.fits"
                # telfile_tmp = "work_place/" + telfile.split("/")[-1]
                stdout = open(teldir_frame + ("telluric_log_m%d.txt" % j), "w")
                sys.stdout = stdout

                # copyfile = "work_place/" + targetspecflux_m[j].split("/")[-1].replace(".fits", "-cp.fits")
                # shutil.copyfile(targetspecflux_m[j], copyfile)
                # reffile = "work_place/" + refspecflux_m[j].split("/")[-1].replace(".fits", "-cp.fits")
                # shutil.copyfile(refspecflux_m[j], reffile)

                # print("oooooooooooooooooooooooo")
                # print(teldir_frame)
                # print(targetspecflux_m[j].split("/")[-1])
                # print(copyfile, telfile_tmp, reffile)
                # iraf.telluric(copyfile, telfile_tmp, reffile, shift=shift_m[j] / dw,
                #               scale=scale_cur,
                #               interactive="no", threshold="0.001", tweakrm="no", xcorr="no")
                telluricAstropy(targetspecflux_m[j], refspecflux_m[j], telfile, shift=shift_m[j] / dw, scale=scale_cur)
                # iraf.telluric(targetspecflux_m[j], telfile, refspecflux_m[j], shift=shift_m[j] / dw,
                #               scale=scale_cur,
                #               interactive="no", threshold="0.001", tweakrm="no", xcorr="no")
                # shutil.copy(telfile_tmp, telfile)
                # os.remove(copyfile)
                # os.remove(reffile)
                # os.remove(telfile_tmp)

                spx_a, spy_a, _, _, _ = openspecfits(telfile)
                minlam = min(spx_a)
                maxlam = max(spx_a)

                ax1 = fig.add_subplot(1, 1, 1)
                max_t, max_r = numpy.amax(spy_t), numpy.amax(spy_r)
                norm_t, norm_r = numpy.median(spy_t[spy_t > max_t * 0.5]), numpy.median(spy_r[spy_r > max_r * 0.5])

                spectrum_plot(ax1, [spx_t, spx_r, spx_a],
                              [spy_t / norm_t, spy_r / norm_r, spy_a / numpy.median(spy_a)],
                              [minlam, maxlam], [0., 2.0], labels=["Target", "Telluric", "Result"],
                              yshift=[0, 0.2, 0.4], xaxis_label=r"Wavelength ($\AA$)", yaxis_label=r"Normalized flux",
                              grid_flag=True, legend_flag=True)
                plt.title("m=%d, shift=%.3f, scale=%.3f" % (j, shift_m[j] / dw, scale_cur))
                plt.savefig(pp, format="pdf")
                plt.clf()

                valuetext = "'%s',%d,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,'%s','%s'" % (
                    telluricID, j, scale_cur, scale_err_cur, shift_m[j], shift_err_m[j], minlam, maxlam,
                    targetframe[i], telfile)
                cur.execute(
                    "INSERT IGNORE INTO telluricresult (telluricID,echelleorder,scale,scaleerr,shift,shifterr,lambdamin,lambdamax,frame,telluricfilepath) VALUES (%s);" % valuetext)
                stdout.close()
            os.chdir(curdir)
            pp.close()

    else:
        cur.execute(
            "INSERT IGNORE INTO telluriccorrection (telluricID,pipelineIDobj,pipelineIDtel,autoflag,advanced,pipelinever,mode,telluricNumber,telluricPath) VALUES (%s);" % valuetext)
        for i in range(len(targetframe)):
            teldir_frame = telluricpath + targetframe[i] + "/"
            os.makedirs(teldir_frame)

            pp = PdfPages(teldir_frame + telluricID + "_" + targetframe[i] + ".pdf")

            m, shift_m, shift_err_m = waveshift_dif(target_wspath[i], ref_wspath, fig, pp)

            m_common = []
            targetspec_m = {}
            refspec_m = {}
            targetspecflux_m = {}
            refspecflux_m = {}
            for j in m:
                tid = -1
                rid = -1
                for l in range(len(targetspec[i])):
                    if targetspec[i][l].split("/")[-1].find("m%d" % j) != -1:
                        tid = l
                        break
                for k in range(len(refspec)):
                    if refspec[k].split("/")[-1].find("m%d" % j) != -1:
                        rid = k
                        break
                if tid != -1 and rid != -1:
                    m_common.append(j)
                    targetspec_m[j] = targetspec[i][tid]
                    refspec_m[j] = refspec[rid]
                    targetspecflux_m[j] = targetspecflux[i][tid]
                    refspecflux_m[j] = refspecflux[rid]

            # ln_tnorm = "target_%s_%s_norm" % (targetframe[i], vacorair)
            # ln_tflux = "target_%s_%s_flux" % (targetframe[i], vacorair)
            # ln_rnorm = "ref_sum_%s_norm" % (vacorair)
            # ln_rflux = "ref_sum_%s_flux" % (vacorair)


            wf = open(teldir_frame + "telluric_pyraf.txt", "w")

            targetflist_m = {}
            refflist_m = {}
            targetnlist_m = {}
            refnlist_m = {}
            for j in m_common:
                targetflist_m[j] = targetspecflux_m[j].split("/")[-1]
                targetnlist_m[j] = targetspec_m[j].split("/")[-1]
                refnlist_m[j] = refspec_m[j].split("/")[-1]
                refflist_m[j] = refspecflux_m[j].split("/")[-1]

            for j in m_common:
                os.symlink(targetspec_m[j], teldir_frame + targetnlist_m[j])
                os.symlink(targetspecflux_m[j], teldir_frame + targetflist_m[j])
                os.symlink(refspec_m[j], teldir_frame + refnlist_m[j])
                os.symlink(refspecflux_m[j], teldir_frame + refflist_m[j])

                spx_t, spy_t, _, dw, _ = openspecfits(targetspecflux_m[j])

                telfile = targetflist_m[j].rstrip("fits").rstrip(".") + "_tel.fits"

                wf.write("telluric %s %s %s thres=0.001 shift=%.3f tweak=no xcorr=no > telluric_log_m%d.txt\n" % (
                targetflist_m[j], telfile, refflist_m[j], shift_m[j] / dw, j))


            pp.close()
            wf.close()

        print("Go ahead!")
        print("cd %s" % telluricpath)

    conn.commit()
    conn.close()
