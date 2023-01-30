import sys
import mysql.connector
from urllib.parse import urlparse
import numpy
from waveshift_measure import read_waveshift
import astropy.io.fits as fits
import os
import glob
import scipy.constants
import time
from vac2air_spec import vac2air_spec
from dopcor_vcorr import rvcorrect, dopcor_rvcorrect, read_rvfile
from telluric_auto import GetListAdvancedFits
from open_mysql_project import openproject
from atran_model_npz import modelnpzopen


def fitslambdacorrection(inputsp, outputsp, lamc, a, shiftm, offset):
    spfits = fits.open(inputsp)
    crval1 = float(spfits[0].header["CRVAL1"])
    cdelt1 = float(spfits[0].header["CDELT1"])

    spfits[0].header["CRVAL1"] = crval1 * (1. + a) - a * lamc + shiftm + offset
    spfits[0].header["CDELT1"] = cdelt1 * (1. + a)
    spfits[0].header["CD1_1"] = cdelt1 * (1. + a)

    spfits.writeto(outputsp)

    spfits.close()


if __name__ == "__main__":
    conn, cur = openproject()

    fsr = "fsr1.30"
    vacorair = "VAC"
    rvfile = "rvcorrect.txt"
    nite = 5
    lowsig = 2
    telmodelhelio = "telluric_model_AIR_helio.npz"

    cur.execute(
        "select pipelineID, FrameNum, totalSNR, mode, obsdate, path, pipelinever from datareduction where pipelineID = '%s';" %
        sys.argv[1])
    rows = cur.fetchall()

    telluricflag = False

    if len(rows) == 1:
        ppid = rows[0][0]
        fnum = int(rows[0][1])
        path = rows[0][5]
        pver = rows[0][6]
    elif rows == []:
        cur.execute(
            "select y.pipelineID, y.FrameNum, x.telluricPath, y.pipelinever, y.path from telluriccorrection as x join datareduction as y on x.pipelineIDobj=y.pipelineID where x.telluricID='%s';" %
            sys.argv[1])
        rows = cur.fetchall()
        if rows == []:
            print(sys.argv[1] + " could not be found.")
            sys.exit()
        elif len(rows) == 1:
            telluricflag = True
            ppid = rows[0][0]
            fnum = int(rows[0][1])
            telpath = rows[0][2]
            pver = rows[0][3]
            path = rows[0][4]
        cur.execute(
            "select y.pipelineID, y.pipelinever, y.path, x.advanced from telluriccorrection as x join datareduction as y on x.pipelineIDtel=y.pipelineID where x.telluricID='%s';" %
            sys.argv[1])
        rows = cur.fetchall()
        if rows == []:
            print(sys.argv[1] + " could not be found.")
            sys.exit()
        elif len(rows) == 1:
            ppid_tel = rows[0][0]
            pver_tel = rows[0][1]
            path_tel = rows[0][2]
            advanced = int(rows[0][3])
        cur.execute(
            "select echelleorder, shift from telluricresult where telluricID='%s' and frame='sum';" % sys.argv[1])
        rows = cur.fetchall()
        m_tel = [int(i[0]) for i in rows]
        shift_tel = {}
        for i in range(len(rows)):
            shift_tel[m_tel[i]] = float(rows[i][1])


    else:
        print("Error: multiple records hit for pipelineID=" + sys.argv[1])
        sys.exit()

    # if not "3.5" in pver:
    #     print("Check Pipeline ver.")
    #     sys.exit()

    rf = open(sys.argv[2], "r")
    rl = rf.readlines()
    rf.close()

    am = {}
    shiftpa = {}
    for i in rl:
        orderpa = int(i.split()[0])
        am[orderpa] = float(i.split()[1])
        shiftpa[orderpa] = float(i.split()[3])

    if not os.path.exists(path + rvfile):
        fits_ex = glob.glob(path + "*_sum/VAC_norm/fsr*/*fits")
        rvcorrect(fits_ex[0], path + rvfile)
    if not os.path.exists(path + telmodelhelio):
        vhelio = read_rvfile(path + rvfile)
        wav, flux, id = modelnpzopen(vacorair="AIR")
        wav_helio = wav * (1. - vhelio / (scipy.constants.c * 1.e-3))
        numpy.savez(path + telmodelhelio, wav=wav_helio, flux=flux, id=id)

    frame = "sum"  # ["NO%d" % (j+1) for j in range(fnum[i])] + ["sum"]
    wsfile = "%s%s%s/%s/%s_%s_%s_norm_%s.txt" % (
        path, "waveshift_measure/", frame, fsr, ppid, fsr, vacorair, frame)
    order, split, wavc, shift, r_edge, r_min, r_dif, r_ratio = read_waveshift(wsfile)
    m = list(set(order))
    m.sort()

    wavc_m = {}
    for k in m:
        wavc_m[k] = numpy.average(wavc[order == k])
    req1 = numpy.array([r_edge[k] > 0.5 for k in range(len(order))])

    clip1 = numpy.array([True for k in range(len(order))])
    for k in range(nite):
        r_ratio_av = numpy.average(r_ratio[numpy.logical_not(req1) & clip1])
        r_ratio_std = numpy.std(r_ratio[numpy.logical_not(req1) & clip1])
        clip1[(r_ratio - r_ratio_av) < - lowsig * r_ratio_std] = False

    req = numpy.logical_not(numpy.logical_or(req1, numpy.logical_not(clip1)))

    wavc_cor = numpy.array([])
    shift_cor = numpy.array([])

    for k in m:
        wavc_cor = numpy.append(wavc_cor, wavc[numpy.logical_and(req, order == k)])
        shift_cor = numpy.append(shift_cor, shift[numpy.logical_and(req, order == k)] - am[k] * (
                wavc[numpy.logical_and(req, order == k)] - wavc_m[k]) - shiftpa[k])

    a, b = numpy.polyfit(wavc_cor, shift_cor, 1)

    if telluricflag:
        frames = glob.glob(telpath + "*")
        frames.sort()
        for fr in frames:
            # os.chdir(fr)
            spfiles = glob.glob(fr + "/*_tel.fits")
            spfiles.sort()
            m = [int(spf.split("_fsr")[-2].split("_m")[-1]) for spf in spfiles]
            for spf in range(len(spfiles)):
                outputf = spfiles[spf].rstrip("fits").rstrip(".") + "_wscor.fits"
                outputf_air = outputf.replace("VAC_", "AIR_")
                outputf_helio = outputf.rstrip("fits").rstrip(".") + "_helio.fits"
                outputf_air_helio = outputf_air.rstrip("fits").rstrip(".") + "_helio.fits"
                if not os.path.exists(outputf): fitslambdacorrection(spfiles[spf], outputf, wavc_m[m[spf]], am[m[spf]],
                                                                     shiftpa[m[spf]], a * wavc_m[m[spf]] + b)
                if not os.path.exists(outputf_air): vac2air_spec(outputf, outputf_air)
                if not os.path.exists(outputf_helio): dopcor_rvcorrect(outputf, outputf_helio, path + rvfile)
                if not os.path.exists(outputf_air_helio): dopcor_rvcorrect(outputf_air, outputf_air_helio,
                                                                           path + rvfile)

        cur.execute(
            "select echelleorder,shift,shifterr,scale,scaleerr,frame,telluricfilepath from telluricresult where telluricID='%s';" %
            sys.argv[1])
        rows = cur.fetchall()
        orders = [int(i[0]) for i in rows]
        shift = [float(i[1]) for i in rows]
        shifterr = [float(i[2]) for i in rows]
        scale = [float(i[3]) for i in rows]
        scaleerr = [float(i[4]) for i in rows]
        frame = [i[5] for i in rows]
        telfile = [i[6] for i in rows]

        if advanced == 0:
            wscor = path + ppid_tel + "_wscor"
            if not os.path.exists(wscor):
                os.makedirs(wscor)
                spfiles = glob.glob(path_tel + "*_sum/VAC_norm/fsr1.30/*_norm.fits")
                spfiles.sort()
                m = [int(spf.split("_fsr")[-2].split("_m")[-1]) for spf in spfiles]
                for spf in range(len(spfiles)):
                    outputf = wscor + "/" + spfiles[spf].split("/")[-1].rstrip("fits").rstrip(".") + "_wscor.fits"
                    outputf_air = outputf.replace("VAC_", "AIR_")
                    outputf_helio = outputf.rstrip("fits").rstrip(".") + "_helio.fits"
                    outputf_air_helio = outputf_air.rstrip("fits").rstrip(".") + "_helio.fits"
                    if not os.path.exists(outputf): fitslambdacorrection(spfiles[spf], outputf, wavc_m[m[spf]],
                                                                         am[m[spf]], shiftpa[m[spf]],
                                                                         a * wavc_m[m[spf]] + b - shift_tel[m[spf]])
                    if not os.path.exists(outputf_air): vac2air_spec(outputf, outputf_air)
                    if not os.path.exists(outputf_helio): dopcor_rvcorrect(outputf, outputf_helio, path + rvfile)
                    if not os.path.exists(outputf_air_helio): dopcor_rvcorrect(outputf_air, outputf_air_helio,
                                                                               path + rvfile)
        else:
            wscor = path + ppid_tel + "_adv_wscor"
            if not os.path.exists(wscor):
                os.makedirs(wscor)
                spfiles = GetListAdvancedFits(path_tel, "fsr?.??")
                spfiles.sort()
                m = [int(spf.split("_m")[-1].rstrip("fits").rstrip(".")) for spf in spfiles]

                for spf in range(len(spfiles)):
                    outputf = wscor + "/" + spfiles[spf].split("/")[-1].rstrip("fits").rstrip(".") + "_VAC_wscor.fits"
                    outputf_air = outputf.replace("VAC_", "AIR_")
                    outputf_helio = outputf.rstrip("fits").rstrip(".") + "_helio.fits"
                    outputf_air_helio = outputf_air.rstrip("fits").rstrip(".") + "_helio.fits"
                    if not os.path.exists(outputf): fitslambdacorrection(spfiles[spf], outputf, wavc_m[m[spf]],
                                                                         am[m[spf]], shiftpa[m[spf]],
                                                                         a * wavc_m[m[spf]] + b - shift_tel[m[spf]])
                    if not os.path.exists(outputf_air): vac2air_spec(outputf, outputf_air)
                    if not os.path.exists(outputf_helio): dopcor_rvcorrect(outputf, outputf_helio, path + rvfile)
                    if not os.path.exists(outputf_air_helio): dopcor_rvcorrect(outputf_air, outputf_air_helio,
                                                                               path + rvfile)


    else:
        wcdir = path + "waveshift_correct/"
        if os.path.exists(wcdir):
            print("Already done previously.")
            sys.exit()
        frames = ["NO%d" % (f + 1) for f in range(fnum)] + ["sum"]

        os.makedirs(wcdir)
        for fr in frames:
            os.makedirs(wcdir + fr)

        for fr in frames:
            fsrdir_norm = glob.glob(path + "*_%s/VAC_norm/fsr*" % fr)
            fsrdir_norm.sort()
            fsrdir_flux = glob.glob(path + "*_%s/VAC_flux/fsr*" % fr)
            fsrdir_flux.sort()
            for fsr in fsrdir_norm:
                os.makedirs(wcdir + fr + "/" + fsr.split("/")[-1])
                spfiles = glob.glob(fsr + "/*_norm.fits")
                spfiles.sort()
                m = [int(spf.rstrip("_VAC_norm.fits").rstrip(fsr.split("/")[-1]).rstrip("_").split("_m")[-1]) for spf in
                     spfiles]

                for spf in range(len(spfiles)):
                    outputf = wcdir + fr + "/" + fsr.split("/")[-1] + "/" + spfiles[spf].split("/")[-1].rstrip(
                        "fits").rstrip(".") + "_wscor.fits"
                    outputf_air = outputf.replace("VAC_", "AIR_")
                    outputf_helio = outputf.rstrip("fits").rstrip(".") + "_helio.fits"
                    outputf_air_helio = outputf_air.rstrip("fits").rstrip(".") + "_helio.fits"
                    if not os.path.exists(outputf): fitslambdacorrection(spfiles[spf], outputf, wavc_m[m[spf]],
                                                                         am[m[spf]], shiftpa[m[spf]],
                                                                         a * wavc_m[m[spf]] + b)
                    if not os.path.exists(outputf_air): vac2air_spec(outputf, outputf_air)
                    if not os.path.exists(outputf_helio): dopcor_rvcorrect(outputf, outputf_helio, path + rvfile)
                    if not os.path.exists(outputf_air_helio): dopcor_rvcorrect(outputf_air, outputf_air_helio,
                                                                               path + rvfile)

            for fsr in fsrdir_flux:
                spfiles = glob.glob(fsr + "/*VAC.fits")
                spfiles.sort()
                m = [int(spf.rstrip("_VAC.fits").rstrip(fsr.split("/")[-1]).rstrip("_").split("_m")[-1]) for spf in
                     spfiles]

                for spf in range(len(spfiles)):
                    outputf = wcdir + fr + "/" + fsr.split("/")[-1] + "/" + spfiles[spf].split("/")[-1].rstrip(
                        "fits").rstrip(".") + "_wscor.fits"
                    outputf_air = outputf.replace("VAC_", "AIR_")
                    outputf_helio = outputf.rstrip("fits").rstrip(".") + "_helio.fits"
                    outputf_air_helio = outputf_air.rstrip("fits").rstrip(".") + "_helio.fits"
                    if not os.path.exists(outputf): fitslambdacorrection(spfiles[spf], outputf, wavc_m[m[spf]],
                                                                         am[m[spf]], shiftpa[m[spf]],
                                                                         a * wavc_m[m[spf]] + b)
                    if not os.path.exists(outputf_air): vac2air_spec(outputf, outputf_air)
                    if not os.path.exists(outputf_helio): dopcor_rvcorrect(outputf, outputf_helio, path + rvfile)
                    if not os.path.exists(outputf_air_helio): dopcor_rvcorrect(outputf_air, outputf_air_helio,
                                                                               path + rvfile)

    conn.close()
