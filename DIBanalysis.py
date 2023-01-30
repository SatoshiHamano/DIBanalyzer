# -*- coding:utf-8 -*-

# from pyraf import iraf
import math
import sys
import numpy
import scipy.constants
import copy
import glob
import os
import time
import matplotlib.pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
from scipy import special
from astropy.modeling.polynomial import Chebyshev1D, Legendre1D
from specutils.fitting.continuum import fit_continuum
from astropy.stats import sigma_clip
from astropy import units as u
from specutils import Spectrum1D
from Spec1Dtools import openspecfits, savespecfits
from spectra_plotter import transmittance_resampling, spectrum_plot
from open_mysql_project import openproject
from combine_MySQL import openTelluricSpectra, obtainCombinePath, pycombine, saveTelluricSpectra, obtainWeight
from DIBcut import readDIBcut_vel
from add_lineDIB_mysql import GetDIBLine
from vac2air_spec import air2vac

# iraf.imred()
# iraf.echelle()


def maskRegion_line(objfits, linelam, velc, fwhm=3.):
    lightvel = scipy.constants.c * 1.e-3
    objspx, objspy, _, _, _ = openspecfits(objfits)
    minlam, maxlam = numpy.amin(objspx), numpy.amax(objspx)

    linebool = objspx > -1

    if type(linelam) == float:
        linelam = [linelam]

    for lam0 in linelam:
        if minlam < lam0 < maxlam:
            linebool *= numpy.absolute(objspx - lam0 * (1. + velc / lightvel)) > fwhm

    linebool[0] = True
    linebool[-1] = True

    return linebool


def maskRegion_tel_core(objfits, telfits, thres=0.5, thres_upside="INDEF"):
    objspx, objspy, _, _, _ = openspecfits(objfits)
    [telspx, telspy] = telfits
    telspy_tm = transmittance_resampling(objspx, objspy, telspx, telspy)
    if thres_upside == "INDEF":
        telbool = telspy_tm > thres
    else:
        telbool = (telspy_tm > thres) & (telspy_tm < thres_upside)

    telbool[0] = True
    telbool[-1] = True

    return telbool


def maskRegion_tel(objfits, telspec, eorder, mode, thres=0.5, thres_upside="INDEF"):
    maskbool_tel = maskRegion_empty(objfits)
    if mode == "obs":
        for i in telspec.keys():
            maskbool_tel *= maskRegion_tel_core(objfits, telspec[i][eorder], thres=thres, thres_upside=thres_upside)  # tel_model case is not considered here.
    elif mode == "model":
        for i in telspec.keys():
            maskbool_tel *= maskRegion_tel_core(objfits, telspec[i], thres=thres, thres_upside=thres_upside)

    return maskbool_tel


def maskRegion_empty(objfits):
    objspx, objspy, _, _, _ = openspecfits(objfits)
    return numpy.array([True for i in range(objspx.size)])


def sampleRegion(objfits, maskbool, minsep=0., upplimit=13400.):
    objspx, objspy, _, _, _ = openspecfits(objfits)
    if numpy.amax(objspx) > upplimit:
        waverange = objspx < upplimit
        objspx = objspx[waverange]
        objspy = objspy[waverange]

    normalize_range = []
    prewave = 0.
    curstart = 0.
    closed = False
    for i in range(1, len(objspx)):
        if maskbool[i]:
            if prewave != objspx[i - 1]:
                if (prewave - curstart > minsep) or curstart == 0.:
                    if prewave != 0.:
                        endwave = prewave
                        normalize_range.append((startwave * u.AA, endwave * u.AA))
                    startwave = objspx[i]
                    curstart = objspx[i]
                    closed = False
            elif i == len(objspx) - 1:
                endwave = objspx[i]
                normalize_range.append((startwave * u.AA, endwave * u.AA))
                closed = True
            prewave = objspx[i]

    return normalize_range


def IRAF_sampleRegion(objfits, maskbool, minsep=0., upplimit=13400.):
    objspx, objspy, _, _, _ = openspecfits(objfits)
    if numpy.amax(objspx) > upplimit:
        waverange = objspx < upplimit
        objspx = objspx[waverange]
        objspy = objspy[waverange]

    normalize_range = ""
    prewave = 0.
    curstart = 0.
    closed = False
    for i in range(1, len(objspx)):
        if maskbool[i]:
            if prewave != objspx[i - 1]:
                if (prewave - curstart > minsep) or curstart == 0.:
                    if prewave != 0.:
                        normalize_range += ":%.1f" % prewave
                    normalize_range += ",%.1f" % objspx[i]
                    curstart = objspx[i]
                    closed = False
            elif i == len(objspx) - 1:
                normalize_range += ":%.1f" % objspx[i]
                closed = True
            prewave = objspx[i]

    if not closed:
        normalize_range = normalize_range.rstrip("%.1f" % curstart).rstrip(",")

    return normalize_range.lstrip(",")


def readDIBspec_cut(DIBID, combineID):
    combpath = obtainCombinePath(combineID)

    fitsfiles = glob.glob("%sDIBs/DIBID%d/%s*DIBID%d.fits" % (combpath, DIBID, combineID, DIBID))
    orders = [float(i.split("/")[-1].split(combineID + "_m")[1].split("_")[0]) for i in fitsfiles]

    fitsdict = {}
    for i in range(len(fitsfiles)):
        fitsdict[orders[i]] = fitsfiles[i]

    return orders, fitsdict


def readDIBspec_norm(DIBID, combineID):
    combpath = obtainCombinePath(combineID)

    fitsfiles = glob.glob("%sDIBs/DIBID%d/%s*DIBID%d_norm*.fits" % (combpath, DIBID, combineID, DIBID))
    orders = [float(i.split("/")[-1].split(combineID + "_m")[1].split("_")[0]) for i in fitsfiles]
    orders_set = list(set(orders))
    orders_set.sort()

    fitsdict = {}
    for m in orders_set:
        fitsdict[m] = []

    for i in range(len(fitsfiles)):
        fitsdict[orders[i]].append(fitsfiles[i])

    return orders_set, fitsdict


def obtainDIBdirpath(DIBID, combineID):
    combpath = obtainCombinePath(combineID)
    DIBpath = "%sDIBs/DIBID%d/" % (combpath, DIBID)
    return DIBpath


def measureSNratio(dibfits, telfits, maskbool_tel, maskbool_line, outputfig, sigma=2., ite=10, saveresult="INDEF"):
    # spx, spy, _, _, _ = openspecfits(inputfits)

    spx, spy, _, dx, _ = openspecfits(dibfits)

    maskbool = maskbool_tel * maskbool_line

    std = numpy.std(spy[maskbool])
    stdprev = copy.copy(std)
    for i in range(ite):
        stdbool = numpy.absolute(spy - 1.) < std * sigma
        std = numpy.std(spy[maskbool * stdbool])
        if stdprev == std:
            break
        else:
            stdprev = std

    SNR = 1. / std

    spy_interp = interpolateSpec(spy, maskbool_tel)
    [telspx, telspy] = telfits
    telspy_res = transmittance_resampling(spx, spy, telspx, telspy)
    error = numpy.ones(spx.size) / SNR / (telspy_res ** 0.5)
    sigma = 1. / SNR
    error_interp = interpolateSpecErr(error, maskbool_tel)
    error_cont = dx * sigma * spy_interp / (1. - sigma ** 2.)

    if saveresult != "INDEF":
        numpy.savez(saveresult, mask=maskbool, wav=spx, fluxinterpolate=spy_interp, error=error,
                    errorinterpolate=error_interp, errorcontinuum=error_cont, SNR=SNR)

    plotForCheckSNR(dibfits, telfits, maskbool * stdbool, outputfig, error, SNR, sigma)

    print("SNR: %.1f" % SNR)

    return SNR


def openDIBresult(saveresult):
    result = numpy.load(saveresult)
    mask = result["mask"]
    wav = result["wav"]
    spy_interp = result["fluxinterpolate"]
    error = result["error"]
    error_interp = result["errorinterpolate"]
    error_cont = result["errorcontinuum"]
    sn = result["SNR"]

    return [mask, wav, spy_interp, sn, error, error_interp, error_cont]


def continuum(inputfits, outputfits, window="INDEF", order=5, high_rej=3., low_rej=2., nite=10):
    # plt.figure()
    # pp = PdfPages("continuum-test.pdf")
    inputspec = Spectrum1D.read(inputfits)
    contwindow = [(inputspec.wavelength[0], inputspec.wavelength[-1])] if window == "INDEF" else window

    sigmask = numpy.array([False for n in range(inputspec.wavelength.size)])
    for i in range(nite):
        clippedspec = Spectrum1D(spectral_axis=inputspec.wavelength[sigmask==False], flux=inputspec.flux[sigmask==False])
        fitted_continuum = fit_continuum(clippedspec, window=contwindow, model=Legendre1D(order))#, window=region)
        contspec = fitted_continuum(inputspec.wavelength)
        normspec = inputspec / Spectrum1D(spectral_axis=inputspec.wavelength, flux=contspec)
        sigmaskspec = sigma_clip(normspec.flux, sigma_lower=low_rej, sigma_upper=high_rej)
        sigmask = sigmaskspec.mask
        # plt.plot(inputspec.wavelength, inputspec.flux, label="Orig")
        # plt.scatter(inputspec.wavelength[sigmask], inputspec.flux[sigmask])
        # plt.plot(inputspec.wavelength, contspec)
        # plt.title("{}, {}".format(i, numpy.sum(sigmask)))
        # plt.savefig(pp, format="pdf")
        # plt.clf()

    savespecfits(normspec, inputfits, outputfits)
    # pp.close()


def IRAF_continuum(inputfits, outputfits, sample="*", override="yes", interactive="no", func="legendre", high_rej=3.,
              low_rej=2., order=5.):
    iraf.continuum(inputfits, outputfits, sample=sample, override=override, interactive=interactive, func=func,
                   high_rej=high_rej, low_rej=low_rej, order=order)


def indexClustering(boolarr, startid=0, falseid=-1):
    idarr = numpy.zeros(boolarr.shape) + falseid
    clid = startid - 1
    naxis = boolarr.size
    curstatus = False
    clusterids = []
    for i in range(1, naxis - 1):
        if boolarr[i] and (not curstatus):
            clid += 1
            curstatus = True
            clusterids.append(clid)
        if boolarr[i]:
            idarr[i] = clid
        if not boolarr[i]:
            curstatus = False

    return idarr, clusterids


def interpolateSpec(spy, maskarr, functype="linear"):
    x = numpy.arange(spy.size)
    f = interpolate.interp1d(x[maskarr], spy[maskarr], kind=functype)
    return f(x)


#
# The error of the interpolated pixel is accumulated to the pixel that is used for the interpolation.
#
# def interpolateSpecErr(err, maskarr):
#     errit = copy.copy(err)
#     itregion = False
#     counter = 0
#     startid = 0
#     for i in range(1, len(err) - 1):
#         if not maskarr[i]:
#             counter += 1
#             itregion = True
#             errit[i] = 0
#         else:
#             if itregion:
#                 errit[startid] += err[startid] * counter / 2.
#                 errit[i] += err[i] * counter / 2.
#             itregion = False
#             counter = 0
#             startid = i
#
#     return errit


def interpolateSpecErr(err, maskarr):
    errit = copy.copy(err)
    itregion = False
    counter = 0
    startid = 0
    for i in range(1, len(err) - 1):
        if not maskarr[i]:
            counter += 1
            itregion = True
        else:
            if itregion:
                for j in range(counter):
                    errit[startid + j + 1] = (err[startid] * (j + 1) + err[i] * (counter - j - 1)) / counter
                # errit[startid] += err[startid] * counter / 2.
                # errit[i] += err[i] * counter / 2.
            itregion = False
            counter = 0
            startid = i

    return errit


def telluricCombine(telfits1, telfits2, combinedx):
    [spx1, spy1] = telfits1
    [spx2, spy2] = telfits2
    telspx = numpy.concatenate((spx1, spx2))
    telspy = numpy.concatenate((spy1, spy2))
    f = interpolate.interp1d(telspx, telspy, kind="linear")
    combinedy = f(combinedx)

    return combinedy


def DIBrange(dibfits, maskbool, resultnpz, linelam, velc, dibsigma=2., intsigma=0.5, vlim=50.):
    lightvel = scipy.constants.c * 1.e-3
    spx, spy, _, _, _ = openspecfits(dibfits)
    results = openDIBresult(resultnpz)
    [_, _, spy_interp, sn, error, error_interp, error_cont] = results

    # spy_interp = interpolateSpec(spy, maskbool)

    spy_minus = numpy.roll(spy_interp, -1)
    spy_plus = numpy.roll(spy_interp, +1)
    spx_minus = numpy.roll(spx, -1)
    spx_plus = numpy.roll(spx, +1)
    dx = (spx_plus - spx_minus) / 2.

    absrange = (spy_interp < 1. - dibsigma * error_interp) & (spy_minus < 1.) & (spy_plus < 1.)
    negarange = spy_interp < 1. - intsigma * error_interp

    # pp = PdfPages("test_%s.pdf" % dibfits.split("/")[-1].split(".")[0])
    # plt.figure()
    # plt.plot(spx, spy_interp)
    # plt.plot(spx, 1. - dibsigma * error_interp)
    # plt.savefig(pp, format="pdf")
    # pp.close()

    absidarr, absids = indexClustering(absrange)

    if absids == []:
        print("No %d sigma absorption was detected." % dibsigma)
        return None

    abscenter = numpy.array([numpy.average(spx[absidarr == i]) for i in absids])
    absvel = (abscenter - linelam) / linelam * lightvel
    dibcandid = absids[numpy.argmin(numpy.absolute(absvel - velc))]

    rangepdf = dibfits.rstrip("fits").rstrip(".") + "_range.pdf"

    if numpy.absolute(absvel[dibcandid] - velc) > vlim:
        print("The velocity of candidate, %.2f km/s, was too apart from input velocity (%.2f km/s)." % (
            absvel[dibcandid], velc))
        plotForCheckIntrange(dibfits, maskbool, absidarr, absids, abscenter, absvel, dibcandid, rangepdf)
        return None

    negidarr, negaids = indexClustering(negarange)
    intid = numpy.median(negidarr[absidarr == dibcandid])
    intrange = numpy.array([False for i in range(spx.size)])
    intrange[negidarr == intid] = True

    plotForCheckIntrange(dibfits, maskbool, absidarr, absids, abscenter, absvel, dibcandid, rangepdf, intrange=intrange)

    return intrange


def diffProb(y1, y2, s1, s2):
    # probability that y1 is larger than y2.

    normy = (y1 - y2) / (s1 ** 2. + s2 ** 2.) ** 0.5 / 2. ** 0.5
    return 0.5 * (1. + special.erf(normy))


def diffProbDensity(y1, s1, dp):
    return 1. / math.sqrt(2 * math.pi) / s1 * numpy.exp(-((y1 - dp) / s1) ** 2. / 2.)


def DIBmeasurement(dibfits, maskbool, intrange, resultnpz, contfactor=0.5):
    spx, spy, _, dx, _ = openspecfits(dibfits)
    results = openDIBresult(resultnpz)
    [_, _, spy_interp, sn, error, error_interp, error_cont] = results
    # spy_interp = interpolateSpec(spy, maskbool)
    # [telspx, telspy] = telfits
    # telspy_res = transmittance_resampling(spx, spy, telspx, telspy)
    # error = numpy.ones(spx.size) / SNratio / (telspy_res ** 0.5)
    # sigma = 1. / SNratio
    # error_interp = interpolateSpecErr(error, maskbool)
    # error_cont = dx * sigma * spy_interp / (1. - sigma ** 2.)

    ew = numpy.sum(1. - spy_interp[intrange]) * dx * 1.e+3  # mA unit
    errorstat_sum = numpy.sum((error_interp[intrange] * dx) ** 2.) ** 0.5 * 1.e+3  # mA unit
    errorcont_sum = numpy.sum(error_cont[intrange] ** 2.) ** 0.5 * contfactor * 1.e+3  # mA unit
    ewerr = (errorstat_sum ** 2. + errorcont_sum ** 2.) ** 0.5  # mA unit
    optdepth = numpy.log(1. / spy_interp)

    lamcenter = numpy.sum(spx[intrange] * optdepth[intrange]) / numpy.sum(optdepth[intrange])
    stdev = math.sqrt(
        numpy.sum(spx[intrange] ** 2. * optdepth[intrange]) / numpy.sum(optdepth[intrange]) - lamcenter ** 2.)

    maxProb = []
    spx_reg = spx[intrange * maskbool]
    spy_reg = spy[intrange * maskbool]
    error_reg = error[intrange * maskbool]
    for i in range(len(spx_reg)):
        difP = diffProb(1. - spy_reg[i], 1. - spy_reg, error_reg[i], error_reg)
        maxProb.append(numpy.prod(difP))

    maxProb = numpy.array(maxProb)
    spx_peak = numpy.sum(spx_reg * maxProb) / numpy.sum(maxProb)
    f = interpolate.interp1d(spx_reg, spy_reg, kind="cubic")
    depth_peak = 1. - f(spx_peak)

    spx_short = spx_reg[spx_reg < spx_peak]
    spx_long = spx_reg[spx_reg >= spx_peak]
    depth_short = 1. - spy_reg[spx_reg < spx_peak]
    depth_long = 1. - spy_reg[spx_reg >= spx_peak]
    error_short = error_reg[spx_reg < spx_peak]
    error_long = error_reg[spx_reg >= spx_peak]
    half_depth = depth_peak / 2.

    probHalfMax_short = diffProbDensity(depth_short, error_short,
                                        half_depth)  # norm((depth_short - half_depth) / error_short)
    probHalfMax_long = diffProbDensity(depth_long, error_long,
                                       half_depth)  # norm((depth_long - half_depth) / error_long)

    spx_fwhm_short = numpy.sum(spx_short * probHalfMax_short) / numpy.sum(probHalfMax_short)
    spx_fwhm_long = numpy.sum(spx_long * probHalfMax_long) / numpy.sum(probHalfMax_long)
    fwhm = spx_fwhm_long - spx_fwhm_short

    print("FWHM", fwhm)
    #
    # if saveresult != "INDEF":
    #     numpy.savez(saveresult, mask=maskbool, wav=spx, fluxinterpolate=spy_interp, error=error)

    return [ew, ewerr, lamcenter, stdev, spx_peak, depth_peak, fwhm]


def DIBupperlimit(dibfits, telfits, maskbool, SNratio, DIBlam, DIBfwhm, vel, contfactor=0.5, UPLsigma=3.):
    lightvel = scipy.constants.c * 1.e-3
    spx, spy, _, dx, _ = openspecfits(dibfits)
    spy_interp = interpolateSpec(spy, maskbool)
    [telspx, telspy] = telfits
    telspy_res = transmittance_resampling(spx, spy, telspx, telspy)
    error = numpy.ones(spx.size) / SNratio / (telspy_res ** 0.5)
    sigma = 1. / SNratio
    error_interp = interpolateSpecErr(error, maskbool)
    error_cont = dx * sigma * spy_interp / (1. - sigma ** 2.)

    center = DIBlam * (1. + vel / lightvel)

    UPLrange = (spx >= center - DIBfwhm) & (spx <= center + DIBfwhm)
    errorstat_sum = numpy.sum((error_interp[UPLrange] * dx) ** 2.) ** 0.5 * 1.e+3  # mA unit
    errorcont_sum = numpy.sum(error_cont[UPLrange] ** 2.) ** 0.5 * contfactor * 1.e+3  # mA unit
    ewerr = (errorstat_sum ** 2. + errorcont_sum ** 2.) ** 0.5
    upperlimit = ewerr * UPLsigma

    return upperlimit


def registerDIBmeasurement(measurementID, measurementNum, DIBID, combineID, order, DIBspecpath, DIBimgpath, primaryflag,
                           autonormalizeflag, automeasurementflag, multiorderflag, EW, EWerr, centerlam_air,
                           centerlam_vac, helio_velocity, lampeak, depth, FWHM, FWHMerr, SNR, integration_start,
                           integration_end,
                           comment):
    conn, cur = openproject()

    cur.execute("INSERT IGNORE INTO DIBmeasurement (measurementID, measurementNum, DIBID, combineID, echelleorder, " +
                "DIBspecpath, DIBimgpath,primaryflag, autonormalizeflag, automeasurementflag, multiorderflag, EW, " +
                "EWerr, centerlam_air, centerlam_vac, helio_velocity, peaklam_air, depth, FWHM, FWHMerr, SNR, integration_start, " +
                "integration_end, comment) VALUES " +
                "('%s',%d,%d,'%s',%.1f,'%s','%s',%s,%s,%s,%s,%.3f,%.3f,%.4f,%.4f,%.3f,%.4f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,'%s');" % (
                    measurementID, measurementNum, DIBID, combineID, order, DIBspecpath, DIBimgpath, primaryflag,
                    autonormalizeflag, automeasurementflag, multiorderflag, EW, EWerr, centerlam_air, centerlam_vac,
                    helio_velocity, lampeak, depth, FWHM, FWHMerr, SNR, integration_start, integration_end, comment))

    conn.commit()
    conn.close()

    if EW != 0.:
        print("%s was registered. (detected)" % measurementID)
    else:
        print("%s was registered. (upper limit)" % measurementID)





def obtainNewMeasurementID(DIBID, combineID, order):
    conn, cur = openproject()
    cur.execute(
        "select measurementID, measurementNum from DIBmeasurement where DIBID=%d and combineID='%s' and echelleorder=%.1f" % (
            DIBID, combineID, order))
    rows = cur.fetchall()
    if rows == []:
        measureNum = 1
    else:
        measureNum = max([int(i[1]) for i in rows]) + 1
    conn.close()

    measurementID = combineID + "_DIB%d_m%.1f_M%d" % (DIBID, order, measureNum)
    return measurementID, measureNum


def obtainMeasuredFits(DIBID, combineID):
    conn, cur = openproject()
    cur.execute(
        "select DIBspecpath from DIBmeasurement where DIBID=%d and combineID='%s'" % (
            DIBID, combineID))
    rows = cur.fetchall()
    if rows == []:
        DIBspeclist = []
    else:
        DIBspeclist = [i[0] for i in rows]
    conn.close()

    return DIBspeclist


def obtainMultiIDsFromMeasurementID(measurementID):
    conn, cur = openproject()

    cur.execute("select x.DIBID, y.objectID, x.combineID, x.echelleorder from DIBmeasurement as x join combinesummary as y "
                "using(combineID) where measurementID='%s';" % measurementID)
    rows = cur.fetchall()
    DIBID = int(rows[0][0])
    objectID = int(rows[0][1])
    combineID = rows[0][2]
    echelleorder = float(rows[0][3])
    print("DIB ID: %d" % DIBID)
    print("object ID: %d" % objectID)
    print("combine ID: %s" % combineID)
    print("order: %.1f" % echelleorder)

    return DIBID, objectID, combineID, echelleorder


def normalizeDIBspec(inputfits, outputfits, telspec, maskbool):
    outputdir = os.path.dirname(outputfits) + "/"
    samplereg = sampleRegion(inputfits, maskbool)
    print("Sample region for %s:" % os.path.basename(inputfits))
    print(samplereg)
    continuum(inputfits, outputfits, window=samplereg)
    plotForCheckNorm(outputfits, telspec, maskbool, "%s.pdf" % (outputfits.rstrip("fits").rstrip(".")))


# def normalizeDIBspec_old(DIBID, combineID):
#     combinepath = obtainCombinePath(combineID)
#     tel_obs, tel_model = openTelluricSpectra(combinepath)
#     vel = readDIBcut_vel(DIBID, combineID)
#     orders, fitsdict = readDIBspec_cut(DIBID, combineID)
#     DIBdir = obtainDIBdirpath(DIBID, combineID)
#     DIBinfo = GetDIBLine(DIBID)
#     DIBlam_rest, DIBfwhm_rest = DIBinfo[1], DIBinfo[4]
#
#     SNlist = []
#     weightlist = []
#     fitsdict_norm = {}
#     for m in orders:
#         output = fitsdict[m].rstrip("fits").rstrip(".") + "_norm.fits"
#         maskbool = maskRegion_line(fitsdict[m], DIBlam_rest, vel, fwhm=DIBfwhm_rest)
#         maskbool_tel = maskRegion_empty(fitsdict[m])
#         if tel_obs != {}:
#             for i in tel_obs.keys():
#                 maskbool_tel *= maskRegion_tel(fitsdict[m], tel_obs[i][m])
#         else:
#             for i in tel_model.keys():
#                 maskbool_tel *= maskRegion_tel(fitsdict[m], tel_model[i][m])
#         maskbool *= maskbool_tel
#         samplereg = sampleRegion(fitsdict[m], maskbool)
#         print(m,samplereg)
#         continuum(fitsdict[m], output, sample=samplereg, func="legendre")
#         fitsdict_norm[m] = output
#         SN = measureSNratio(output, maskbool)
#         SNlist.append(SN)
#         weightlist.append(1. / SN ** 2)
#
#         plotForCheck(fitsdict_norm[m], tel_obs[i][m], maskbool, "%snormalize_m%d.pdf" % (DIBdir, m))
#
#     if len(orders) > 1:
#         combfits = DIBdir + combineID + "_m%.1f_DIBID%d_norm.fits" % (numpy.average(orders), DIBID)
#         pycombine([fitsdict_norm[m] for m in orders], combfits, weightlist)
#         SNcomb = math.sqrt(SNlist[0] ** 2 + SNlist[1] ** 2)
#         SNlist.append(SNcomb)
#         newm = numpy.average(orders)
#         fitsdict_norm[newm] = combfits
#         orders.append(newm)
#         combinedx, _, _, _, _ = openspecfits(combfits)
#         if tel_obs != {}:
#             for i in tel_obs.keys():
#                 tel_obs[i][orders[-1]] = [combinedx, telluricCombine(tel_obs[i][orders[0]], tel_obs[i][orders[1]], combinedx)]
#
#     saveTelluricSpectra(DIBdir, tel_obs, tel_model)


def DIBanalysis(DIBID, combineID, autonorm=True, combine=True, detectionsigma=5., setprimary=False, telthres=0.5):
    c = scipy.constants.c * 1.e-3
    combinepath = obtainCombinePath(combineID)
    vel = readDIBcut_vel(DIBID, combineID)
    DIBdir = obtainDIBdirpath(DIBID, combineID)
    DIBinfo = GetDIBLine(DIBID)
    DIBlam_rest, DIBfwhm = DIBinfo[1], DIBinfo[4]
    datasetID, weight = obtainWeight(combineID)
    maxid = datasetID[weight.index(max(weight))]
    orders_cut, fitsdict = readDIBspec_cut(DIBID, combineID)
    measuredFits = obtainMeasuredFits(DIBID, combineID)

    if os.path.exists(DIBdir + "telluric_spectra.pickle"):
        tel_obs, tel_model = openTelluricSpectra(DIBdir)
    else:
        tel_obs, tel_model = openTelluricSpectra(combinepath)

    if tel_obs != {}:
        telspec = tel_obs
        telmode = "obs"
    else:
        telspec = tel_model
        telmode = "model"

    if autonorm:
        maskbool_dib_fornorm = {}
        maskbool_tel_fornorm = {}
        maskbool_ful_fornorm = {}
        for m in orders_cut:
            output = fitsdict[m].rstrip("fits").rstrip(".") + "_norm.fits"
            if not os.path.exists(output):
                maskbool_dib_fornorm[m] = maskRegion_line(fitsdict[m], DIBlam_rest, vel, fwhm=DIBfwhm)
                maskbool_tel_fornorm[m] = maskRegion_tel(fitsdict[m], telspec, m, telmode)
                # maskbool_tel_fornorm[m] = maskRegion_empty(fitsdict[m])
                # if telmode == "obs":
                #     for i in telspec.keys():
                #         maskbool_tel_fornorm[m] *= maskRegion_tel(fitsdict[m],
                #                                                   telspec[i][
                #                                                       m])  # tel_model case is not considered here.
                # elif telmode == "model":
                #     for i in telspec.keys():
                #         maskbool_tel_fornorm[m] *= maskRegion_tel(fitsdict[m], telspec[i])
                maskbool_ful_fornorm[m] = maskbool_dib_fornorm[m] * maskbool_tel_fornorm[m]

                if telmode == "obs":
                    normalizeDIBspec(fitsdict[m], output, telspec[datasetID[0]][m], maskbool_ful_fornorm[m])
                elif telmode == "model":
                    normalizeDIBspec(fitsdict[m], output, telspec[datasetID[0]], maskbool_ful_fornorm[m])


    orders, fitsdict_norm = readDIBspec_norm(DIBID, combineID)

    maskbool_dib = {}
    maskbool_tel = {}
    maskbool_ful = {}
    resultnpz = {}
    SNlist = {}
    SNfig = {}
    for m in orders:
        maskbool_tel[m] = []
        maskbool_dib[m] = []
        maskbool_ful[m] = []
        resultnpz[m] = []
        SNlist[m] = []
        SNfig[m] = []
        for j in range(len(fitsdict_norm[m])):
            maskbool_dib[m].append(maskRegion_line(fitsdict_norm[m][j], DIBlam_rest, vel, fwhm=DIBfwhm))
            maskbool_tel[m].append(maskRegion_tel(fitsdict_norm[m][j], telspec, m, telmode, thres=telthres))
            # maskbool_tel[m].append(maskRegion_empty(fitsdict_norm[m][j]))
            # for i in telspec.keys():
            #     if telmode == "obs":
            #         maskbool_tel[m][j] *= maskRegion_tel(fitsdict_norm[m][j], telspec[i][m])
            #     elif telmode == "model":
            #         maskbool_tel[m][j] *= maskRegion_tel(fitsdict_norm[m][j], telspec[i])
            maskbool_ful[m].append(maskbool_dib[m][j] * maskbool_tel[m][j])
            resultnpz[m].append(fitsdict_norm[m][j].rstrip("fits") + "npz")
            SNfig[m].append(fitsdict_norm[m][j].rstrip("fits").rstrip(".") + "_snr.pdf")

            if os.path.exists(resultnpz[m][j]):
                snresults = openDIBresult(resultnpz[m][j])
                [_, _, _, sn, _, _, _] = snresults
                SNlist[m].append(sn)
            else:
                if telmode == "obs":
                    SNlist[m].append(
                        measureSNratio(fitsdict_norm[m][j], telspec[maxid][m], maskbool_tel[m][j], maskbool_dib[m][j],
                                       SNfig[m][j], saveresult=resultnpz[m][j]))
                elif telmode == "model":
                    SNlist[m].append(
                        measureSNratio(fitsdict_norm[m][j], telspec[maxid], maskbool_tel[m][j], maskbool_dib[m][j],
                                       SNfig[m][j], saveresult=resultnpz[m][j]))

    if combine and len(orders) == 2:
        combfits = DIBdir + combineID + "_m%.1f_DIBID%d_norm.fits" % (numpy.average(orders), DIBID)
        pycombine([fitsdict_norm[m][0] for m in orders], combfits, [SNlist[m][0] ** 2. for m in orders])
        combinedx, _, _, _, _ = openspecfits(combfits)
        newm = numpy.average(orders)
        fitsdict_norm[newm] = []
        maskbool_dib[newm] = []
        maskbool_tel[newm] = []
        maskbool_ful[newm] = []
        resultnpz[newm] = []
        SNlist[newm] = []
        SNfig[newm] = []
        fitsdict_norm[newm].append(combfits)
        orders.append(newm)
        maskbool_dib[newm].append(maskRegion_line(fitsdict_norm[newm][0], DIBlam_rest, vel, fwhm=DIBfwhm))
        # maskbool_tel[newm].append(maskRegion_empty(fitsdict_norm[newm][0]))
        if telmode == "obs":
            for i in telspec.keys():
                telspec[i][newm] = [combinedx, telluricCombine(telspec[i][orders[0]], telspec[i][orders[1]], combinedx)]
                tel_obs[i][newm] = telspec[i][newm]
                # maskbool_tel[newm][0] *= maskRegion_tel(fitsdict_norm[newm][0], telspec[i][newm])
        # elif telmode == "model":
        #     for i in telspec.keys():
        #         maskbool_tel[newm][0] *= maskRegion_tel(fitsdict_norm[newm][0], telspec[i])
        maskbool_tel[newm].append(maskRegion_tel(fitsdict_norm[newm][0], telspec, newm, telmode, thres=telthres))

        maskbool_ful[newm].append(maskbool_dib[newm][0] * maskbool_tel[newm][0])
        resultnpz[newm].append(DIBdir + combineID + "_m%.1f_DIBID%d_norm.npz" % (newm, DIBID))
        SNfig[newm].append(fitsdict_norm[newm][0].rstrip("fits").rstrip(".") + "_snr.pdf")

        if telmode == "obs":
            SNlist[newm].append(
                measureSNratio(fitsdict_norm[newm][0], telspec[maxid][newm], maskbool_tel[newm][0],
                               maskbool_dib[newm][0], SNfig[newm][0],
                               saveresult=resultnpz[newm][0]))
        elif telmode == "model":
            SNlist[newm].append(
                measureSNratio(fitsdict_norm[newm][0], telspec[maxid], maskbool_tel[newm][0], maskbool_dib[newm][0],
                               SNfig[newm][0], saveresult=resultnpz[newm][0]))

    if not os.path.exists(DIBdir + "telluric_spectra.pickle"):
        saveTelluricSpectra(DIBdir, tel_obs, tel_model)

    SNall = []
    for i in orders:
        for j in range(len(fitsdict_norm[i])):
            SNall.append(SNlist[i][j])
            if SNlist[i][j] >= max(SNall):
                SNhighest = [i, j]

    primaryFlag = {}
    for i in orders:
        primaryFlag[i] = []
        for j in range(len(fitsdict_norm[i])):
            if [i, j] == SNhighest and setprimary:
                primaryFlag[i].append(True)
            else:
                primaryFlag[i].append(False)

    for i in orders:
        print(i, "is being measured.")
        if i.is_integer():
            multiorder = False
        else:
            multiorder = True
        for j in range(len(fitsdict_norm[i])):
            if not fitsdict_norm[i][j] in measuredFits:
                measureID, measureNum = obtainNewMeasurementID(DIBID, combineID, i)
                spx, _, _, _, _ = openspecfits(fitsdict_norm[i][j])

                intrange = DIBrange(fitsdict_norm[i][j], maskbool_tel[i][j], resultnpz[i][j], DIBlam_rest, vel)
                # resultnpz = DIBdir + combineID + "_m%.1f_DIBID%d_norm_M%d.npz" % (i, DIBID, measureNum)
                resultfig = DIBdir + combineID + "_m%.1f_DIBID%d_norm_M%d.pdf" % (i, DIBID, measureNum)

                if intrange is not None and numpy.sum(intrange) >= 5:
                    DIBparam = DIBmeasurement(fitsdict_norm[i][j], maskbool_tel[i][j], intrange, resultnpz[i][j])
                    [ew, ewerr, lamcenter, stdev, spx_peak, depth_peak, fwhm] = DIBparam
                    v_helio = (lamcenter - DIBlam_rest) / DIBlam_rest * c
                    lammin = min(spx[intrange])
                    lammax = max(spx[intrange])
                    registerDIBmeasurement(measureID, measureNum, DIBID, combineID, i, fitsdict_norm[i][j],
                                           "INDEF", primaryFlag[i][j], autonorm, True, multiorder, ew, ewerr, lamcenter,
                                           air2vac(lamcenter),
                                           v_helio, spx_peak, depth_peak, fwhm, 0., SNlist[i][j], lammin, lammax, "")

                    DIBtext = r"$\lambda%d$ (ID=%d)" % (DIBlam_rest, DIBID) + \
                              "\n SNR = %.1f\n EW = %.1f " % (SNlist[i][j], ew) + \
                              r"$\pm$ %.1f m$\AA$" % (ewerr) + "\n" + \
                              r" FWHM = %.2f $\AA$" % (fwhm) + \
                              "\n depth = %.3f\n" % (depth_peak) + \
                              r" v = %.2f km s$^{-1}$" % (v_helio) + \
                              "\n" + r" range = %.2f-%.2f $\AA$" % (lammin, lammax)

                    loadresult = openDIBresult(resultnpz[i][j])
                    [_, _, _, _, _, errorspec, _] = loadresult
                    # errorspec = loadresult["error"]
                    if telmode == "obs":
                        DIBmeasureplot(fitsdict_norm[i][j], telspec[maxid][i], errorspec, maskbool_tel[i][j], intrange,
                                       DIBID, resultfig, measureID, SNlist[i][j], DIBparam=DIBparam, detection=True,
                                       DIBtext=DIBtext)
                    elif telmode == "model":
                        DIBmeasureplot(fitsdict_norm[i][j], telspec[maxid], errorspec, maskbool_tel[i][j], intrange,
                                       DIBID, resultfig, measureID, SNlist[i][j], DIBparam=DIBparam, detection=True,
                                       DIBtext=DIBtext)

                else:
                    lamcenter = DIBlam_rest * (1. + vel / c)
                    if telmode == "obs":
                        upperlimit = DIBupperlimit(fitsdict_norm[i][j], telspec[maxid][i], maskbool_tel[i][j],
                                                   SNlist[i][j], DIBlam_rest, DIBfwhm, vel)
                    elif telmode == "model":
                        upperlimit = DIBupperlimit(fitsdict_norm[i][j], telspec[maxid], maskbool_tel[i][j],
                                                   SNlist[i][j], DIBlam_rest, DIBfwhm, vel)

                    registerDIBmeasurement(measureID, measureNum, DIBID, combineID, i, fitsdict_norm[i][j],
                                           "INDEF", primaryFlag[i][j], autonorm, True, multiorder, 0., upperlimit, lamcenter,
                                           air2vac(lamcenter), vel, 0., 0., 0., 0., SNlist[i][j], lamcenter - DIBfwhm,
                                           lamcenter + DIBfwhm, "")
                    DIBtext = r"$\lambda%d$ (ID=%d)" % (DIBlam_rest, DIBID) + \
                              "\n SNR = %.1f\n EW " % (SNlist[i][j]) + \
                              r"$<$ %.1f m$\AA$" % (upperlimit) + "\n" + \
                              "\n" + r" range = %.2f-%.2f $\AA$" % (lamcenter - DIBfwhm, lamcenter + DIBfwhm)

                    loadresult = openDIBresult(resultnpz[i][j])
                    [_, _, _, _, _, errorspec, _] = loadresult

                    if telmode == "obs":
                        DIBmeasureplot(fitsdict_norm[i][j], telspec[maxid][i], errorspec, maskbool_tel[i][j], intrange,
                                       DIBID, resultfig, measureID, SNlist[i][j], detection=False, DIBtext=DIBtext)
                    elif telmode == "model":
                        DIBmeasureplot(fitsdict_norm[i][j], telspec[maxid], errorspec, maskbool_tel[i][j], intrange,
                                       DIBID, resultfig, measureID, SNlist[i][j], detection=False, DIBtext=DIBtext)

            else:
                print("%s was already measured." % fitsdict_norm[i][j].split("/")[-1])

def DIBmeasureplot(dibfits, telfits, error, maskbool, intrange, DIBID, figfile, title, SN, DIBparam="INDEF",
                   DIBtext="INDEF", detection=True, velocity=False):
    lightvel = scipy.constants.c * 1.e-3
    spx, spy, _, dx, _ = openspecfits(dibfits)
    spy_interp = interpolateSpec(spy, maskbool)
    onespec = numpy.ones(spy.shape)
    [telspx, telspy] = telfits
    telspy_res = transmittance_resampling(spx, spy, telspx, telspy)
    DIBinfo = GetDIBLine(DIBID)
    DIBlam_rest, DIBfwhm = DIBinfo[1], DIBinfo[4]

    if detection and DIBparam != "INDEF":
        [ew, ewerr, lamcenter, stdev, spx_peak, depth_peak, fwhm] = DIBparam

    if velocity:
        spx = (spx - DIBlam_rest) / DIBlam_rest * lightvel

    if figfile.find(".pdf") != -1:
        pp = PdfPages(figfile)
    else:
        pp = figfile

    xrange = [numpy.amin(spx), numpy.amax(spx)]
    uppsig = 5.
    lowsig = 15.

    if detection and DIBparam != "INDEF":
        yrange = [min(1. - depth_peak * 1.2, 1. - lowsig / SN), 1. + uppsig / SN]
    else:
        yrange = [1. - lowsig / SN, 1. + uppsig / SN]

    figw = 0.87
    figh_tel = 0.16
    figh_obj = 0.6
    orig_w = 0.08
    orig_h = 0.1
    gap_telobj = 0.05

    plt.figure(figsize=(12, 6))
    ax1 = plt.axes([orig_w, orig_h, figw, figh_tel])
    spectrum_plot(ax1, [telspx], [telspy], xrange, [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                  yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

    ax2 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj, figw, figh_obj])
    spectrum_plot(ax2, [spx, spx, spx, spx, spx], [onespec, spy, spy_interp, spy_interp + error, spy_interp - error],
                  xrange, yrange, yaxis_label="Normalized flux", colors=["0.5", "0.5", "k", "b", "b"],
                  lines=["--", "-", "-", "-", "-"], linew=[1., 1., 2., 1., 1.])

    if intrange is not None:
        ax2.fill_between(spx[intrange], onespec[intrange], spy_interp[intrange], facecolor="y", alpha=0.5)

    if detection and DIBparam != "INDEF":
        lineleng = (yrange[1] - yrange[0]) / 15.
        if velocity:
            spx_peak = (spx_peak - DIBlam_rest) / DIBlam_rest * lightvel
            fwhm = fwhm / DIBlam_rest * lightvel
        ax2.plot([lamcenter - fwhm / 2., lamcenter + fwhm / 2.], [1. - depth_peak / 2., 1. - depth_peak / 2.],
                 color="orange")
        ax2.plot([lamcenter, lamcenter], [1. - depth_peak - lineleng, 1. - depth_peak - lineleng * 2.], "r", lw=1.5)

    plt.title(title)

    if DIBtext != "INDEF":
        plt.figtext(0.1, 0.37, DIBtext, fontsize=8, bbox=dict(facecolor='white', alpha=0.8))

    if figfile.find(".pdf") != -1:
        plt.savefig(pp, format="pdf")
    else:
        plt.savefig(pp)

    plt.clf()

    if figfile.find(".pdf") != -1:
        pp.close()


def plotForCheckNorm(dibfits, telfits, maskbool, output):
    spx, spy, _, dx, _ = openspecfits(dibfits)
    spy_mask = numpy.ma.masked_where(numpy.logical_not(maskbool), spy)
    [telspx, telspy] = telfits

    pp = PdfPages(output)

    xrange = [numpy.amin(spx), numpy.amax(spx)]

    figw = 0.87
    figh_tel = 0.16
    figh_obj = 0.6
    orig_w = 0.08
    orig_h = 0.1
    gap_telobj = 0.05

    plt.figure(figsize=(12, 6))
    ax1 = plt.axes([orig_w, orig_h, figw, figh_tel])
    spectrum_plot(ax1, [telspx], [telspy], xrange, [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                  yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

    ax2 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj, figw, figh_obj])
    spectrum_plot(ax2, [spx, spx], [spy, spy_mask], xrange, [min(spy[spy > 0.]) - 0.1, max(spy) + 0.3],
                  yaxis_label="Normalized flux", colors=["0.5", "b"], lines=["-", "-"])

    plt.savefig(pp, format="pdf")
    plt.clf()

    pp.close()


def plotForCheckSNR(dibfits, telfits, maskbool, output, errorspec, SNR, sigma):
    spx, spy, _, dx, _ = openspecfits(dibfits)
    spy_mask = numpy.ma.masked_where(numpy.logical_not(maskbool), spy)
    spy_plus = spy + errorspec
    spy_minus = spy - errorspec
    spy_plus_mask = numpy.ma.masked_where(numpy.logical_not(maskbool), spy_plus)
    spy_minus_mask = numpy.ma.masked_where(numpy.logical_not(maskbool), spy_minus)
    snr_plus = numpy.zeros(spx.size) + 1. / SNR
    snr_minus = numpy.zeros(spx.size) - 1. / SNR

    SNRtext = "%s\nSNR=%.1f\nclipping sigma=%.1f" % (dibfits.split("/")[-1], SNR, sigma)

    [telspx, telspy] = telfits

    pp = PdfPages(output)

    xrange = [numpy.amin(spx), numpy.amax(spx)]

    figw = 0.87
    figh_tel = 0.16
    figh_obj = 0.6
    orig_w = 0.08
    orig_h = 0.1
    gap_telobj = 0.05
    inset_orig_w = 0.6
    inset_orig_h = 0.1
    inset_w = 0.25
    inset_h = 0.2

    plt.figure(figsize=(12, 6))
    ax1 = plt.axes([orig_w, orig_h, figw, figh_tel])
    spectrum_plot(ax1, [telspx], [telspy], xrange, [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                  yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

    ax2 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj, figw, figh_obj])
    spectrum_plot(ax2, [spx, spx, spx, spx, spx, spx, spx, spx],
                  [spy, spy_mask, spy_plus, spy_minus, spy_plus_mask, spy_minus_mask, snr_plus, snr_minus], xrange,
                  [min(spy[spy > 0.]) - 0.3, max(spy) + 0.1],
                  yaxis_label="Normalized flux", colors=["0.5", "b", "0.5", "0.5", "cyan", "cyan", "red", "red"],
                  lines=["-", "-", "-", "-", "-", "-", "--", "--"])

    ax3 = plt.axes([orig_w + inset_orig_w, orig_h + inset_orig_h + figh_tel + gap_telobj, inset_w, inset_h])
    ax3.hist([(spy[maskbool] - 1), (spy - 1)], bins=40, histtype="step", range=(-5. * 1. / SNR, 5. * 1. / SNR))
    ax3.plot([-3. / SNR, -3. / SNR], [0, 200], "k--")
    ax3.plot([3. / SNR, 3. / SNR], [0, 200], "k--")
    ax3.plot([-1. / SNR, -1. / SNR], [0, 200], "k--")
    ax3.plot([1. / SNR, 1. / SNR], [0, 200], "k--")
    ax3.set_xlabel("Normalized flux - 1.")
    ax3.set_ylabel(r"$N$")
    ax3.set_xlim(-5. * 1. / SNR, 5. * 1. / SNR)
    ax3.set_ylim(0, numpy.sum(maskbool) * 1.1 / 8.)

    plt.figtext(0.1, 0.37, SNRtext, fontsize=8, bbox=dict(facecolor='white', alpha=0.8))

    plt.savefig(pp, format="pdf")
    plt.clf()

    pp.close()


def plotForCheckIntrange(dibfits, maskbool, absidarr, absids, abscenter, absvel, dibcandid, output, intrange="INDEF"):
    spx, spy, _, dx, _ = openspecfits(dibfits)
    spy_interp = interpolateSpec(spy, maskbool)
    onearr = numpy.ones(spx.shape)

    pp = PdfPages(output)

    xrange = [numpy.amin(spx), numpy.amax(spx)]

    figw = 0.87
    figh_obj = 0.81
    orig_w = 0.08
    orig_h = 0.1

    plt.figure(figsize=(12, 6))
    ax1 = plt.axes([orig_w, orig_h, figw, figh_obj])
    spectrum_plot(ax1, [spx, spx, spx], [spy, spy_interp, onearr], xrange, [min(spy[spy > 0.]) - 0.1, max(spy) + 0.3],
                  yaxis_label="Normalized flux", colors=["0.5", "k", "0.5"], lines=["-", "-", "--"],
                  labels=["Original", "Interpolated", "Continuum level"], legend_flag=True)

    for i in range(len(absids)):
        if dibcandid == absids[i]:
            ax1.plot(spx[absidarr == absids[i]], onearr[absidarr == absids[i]] + 0.1, "r")
            ax1.text(abscenter[i], 1.15, "ID=%d\nv=%.1f" % (absids[i], absvel[i]), color="r", ha="center", va="bottom")
        else:
            ax1.plot(spx[absidarr == absids[i]], onearr[absidarr == absids[i]] + 0.1, "b")
            ax1.text(abscenter[i], 1.15, "ID=%d\nv=%.1f" % (absids[i], absvel[i]), color="b", ha="center", va="bottom")

    if intrange != "INDEF":
        ax1.fill_between(spx[intrange], onearr[intrange], spy_interp[intrange], facecolor="y", alpha=0.5)

    plt.savefig(pp, format="pdf")
    plt.clf()

    pp.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dibid", type=int, help="dibid", nargs='*')
    parser.add_argument("-c", "--combineid", type=str, help="combineid", nargs='*')
    parser.add_argument("-p", "--primary", action="store_true", help="Primary flag")
    parser.add_argument("-t", "--telthres", type=float, default=0.5, help="Transmittance threshold")

    args = parser.parse_args()

    DIBID = args.dibid
    combineID = args.combineid
    primary = args.primary
    telthres = args.telthres

    for comb in combineID:
        print("[combine ID: %s]" % comb)
        for dib in DIBID:
            print("   DIBID=%d" % dib)
            DIBanalysis(dib, comb, setprimary=primary, telthres=telthres, combine=True)
            if len(combineID) * len(DIBID) > 1:
                time.sleep(2)

    # iraf.continuum(filename[0], filename[1], sample=filename[2], override="yes", interactive="yes", func="legendre",
    #               high_rej="3", order="3")
