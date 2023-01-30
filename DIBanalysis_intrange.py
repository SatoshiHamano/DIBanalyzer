# -*- coding:utf-8 -*-

import sys, os
import numpy
import scipy.constants
import argparse
from DIBanalysis import obtainCombinePath, readDIBcut_vel, obtainDIBdirpath, GetDIBLine, obtainWeight, readDIBspec_cut, \
    openTelluricSpectra, readDIBspec_norm, obtainNewMeasurementID, maskRegion_tel, DIBmeasurement, \
    registerDIBmeasurement, openDIBresult, DIBmeasureplot, obtainMultiIDsFromMeasurementID
from Spec1Dtools import openspecfits
from vac2air_spec import air2vac

if __name__ == "__main__":
    c = scipy.constants.c * 1.e-3

    parser = argparse.ArgumentParser()
    parser.add_argument("measurementID", type=str, help="Measurement ID")
    parser.add_argument("intmin", type=float, help="Integration start wavelength")
    parser.add_argument("intmax", type=float, help="Integration end wavelength")
    parser.add_argument("-p", "--primary", action="store_true", help="Primary flag")
    parser.add_argument("-m", "--manual", action="store_true", help="Primary flag")

    args = parser.parse_args()

    measurementID = args.measurementID
    intmin = args.intmin
    intmax = args.intmax
    primary = args.primary
    autonorm = not args.manual

    # filename = sys.argv[1:]
    # measurementID = filename[0]
    DIBID, _, combineID, eorder = obtainMultiIDsFromMeasurementID(measurementID)
    # DIBID = int(filename[0])
    # combineID = filename[1]
    # eorder = float(filename[2])
    # intmin = float(filename[1])
    # intmax = float(filename[2])
    # autonorm = True
    # if len(filename) == 6:
    #     if filename[5] == "manual":
    #         autonorm = False

    combinepath = obtainCombinePath(combineID)
    vel = readDIBcut_vel(DIBID, combineID)
    DIBdir = obtainDIBdirpath(DIBID, combineID)
    DIBinfo = GetDIBLine(DIBID)
    DIBlam_rest, DIBfwhm = DIBinfo[1], DIBinfo[4]
    datasetID, weight = obtainWeight(combineID)
    maxid = datasetID[weight.index(max(weight))]
    orders_cut, fitsdict = readDIBspec_cut(DIBID, combineID)
    orders, fitsdict_norm = readDIBspec_norm(DIBID, combineID)

    if len(fitsdict_norm[eorder]) == 1:
        mainid = 0
    else:
        for i in range(len(fitsdict_norm[eorder])):
            print("ID=%d: " % i, fitsdict_norm[eorder][i].split("/")[-1])
        mainid = int(input("Enter ID: "))

    resultnpz = fitsdict_norm[eorder][mainid].rstrip("fits") + "npz"
    snresults = openDIBresult(resultnpz)
    [_, _, _, sn, _, errorspec, _] = snresults

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

    spx, spy, _, _, _ = openspecfits(fitsdict_norm[eorder][mainid])
    maskbool_tel = maskRegion_tel(fitsdict_norm[eorder][mainid], telspec, eorder, telmode)

    intrange = numpy.array([False for i in range(spx.size)])
    intrange[(spx >= intmin) & (spx <= intmax)] = True

    measureID, measureNum = obtainNewMeasurementID(DIBID, combineID, eorder)
    resultfig = DIBdir + combineID + "_m%.1f_DIBID%d_norm_M%d.pdf" % (eorder, DIBID, measureNum)

    if eorder.is_integer():
        multiorder = False
    else:
        multiorder = True

    DIBparam = DIBmeasurement(fitsdict_norm[eorder][mainid], maskbool_tel, intrange, resultnpz)
    [ew, ewerr, lamcenter, stdev, spx_peak, depth_peak, fwhm] = DIBparam
    v_helio = (lamcenter - DIBlam_rest) / DIBlam_rest * c
    lammin = min(spx[intrange])
    lammax = max(spx[intrange])
    registerDIBmeasurement(measureID, measureNum, DIBID, combineID, eorder, fitsdict_norm[eorder][mainid],
                           "INDEF", primary, autonorm, True, multiorder, ew, ewerr, lamcenter,
                           air2vac(lamcenter),
                           v_helio, spx_peak, depth_peak, fwhm, 0., sn, lammin, lammax, "")

    DIBtext = r"$\lambda%d$ (ID=%d)" % (DIBlam_rest, DIBID) + \
              "\n SNR = %.1f\n EW = %.1f " % (sn, ew) + \
              r"$\pm$ %.1f m$\AA$" % (ewerr) + "\n" + \
              r" FWHM = %.2f $\AA$" % (fwhm) + \
              "\n depth = %.3f\n" % (depth_peak) + \
              r" v = %.2f km s$^{-1}$" % (v_helio) + \
              "\n" + r" range = %.2f-%.2f $\AA$" % (lammin, lammax)

    DIBmeasureplot(fitsdict_norm[eorder][mainid], telspec[maxid][eorder], errorspec, maskbool_tel, intrange,
                   DIBID, resultfig, measureID, sn, DIBparam=DIBparam, detection=True,
                   DIBtext=DIBtext)
