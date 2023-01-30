#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sys
import matplotlib.pyplot as plt
import glob
import pickle
from Spec1Dtools import openspecfits, savespecfits
from specutils import Spectrum1D
import numpy
import copy
import os
import shutil
# from pyraf import iraf
from astropy import units as u
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler
from spectra_plotter import MultiSpecPlotter
from open_mysql_project import openproject
from obtain_filelist import *
from matplotlib.backends.backend_pdf import PdfPages
from add_lineDIB_mysql import *
from atran_model_npz import modelnpzopenHelio, modelnpzopen
from dopcor_vcorr import read_rvfile_ID


def pycombine(splist, output, weightlist):
    if len(splist) != len(weightlist):
        print("Error: the lengths of spectrum list ({}) and weight list ({}) must be equal.".format(len(splist),len(weightlist)))
        return None

    speclist = []
    minlam, maxlam, deltalam, aveflux = [], [], [], []
    for i in splist:
        spec = Spectrum1D.read(i)
        if spec.flux[0] == 0.:
            spec = spec[spec.wavelength[1] - 1e-7 * u.AA:spec.wavelength[-1] + 1e-7 * u.AA]
        if spec.flux[-1] == 0.:
            spec = spec[spec.wavelength[0] - 1e-7 * u.AA:spec.wavelength[-2] + 1e-7 * u.AA]
        speclist.append(spec)
        minlam.append(numpy.amin(spec.wavelength.value))
        maxlam.append(numpy.amax(spec.wavelength.value))
        deltaspec = spec.wavelength.value - numpy.roll(spec.wavelength.value, 1)
        deltalam.append(numpy.median(deltaspec[1:]))
        aveflux.append(numpy.average(spec.flux.value))

    minnewsp = numpy.amin(minlam)
    maxnewsp = numpy.amax(maxlam)
    deltanewsp = numpy.average(deltalam)
    newspecgrid = numpy.arange(minnewsp, maxnewsp, deltanewsp) * speclist[0].wavelength.unit
    spline = SplineInterpolatedResampler(bin_edges='zero_fill')
    resampledSpeclist = []
    for i in range(len(speclist)):
        resampledSpec = spline(speclist[i], newspecgrid)
        resampledSpec.flux.value[resampledSpec.flux.value == 0] += aveflux[i]
        resampledSpeclist.append(resampledSpec * weightlist[i])

    combinedSpec = numpy.sum(resampledSpeclist) / numpy.sum(weightlist)
    savespecfits(combinedSpec, splist[0], output)


def IRAF_pycombine(splist, output, weightlist):
    curdir = os.getcwd()
    # print(curdir)
    outputdir = os.path.dirname(output)
    outputfile = os.path.basename(output)
    # print(outputdir)

    os.chdir(outputdir)

    scombfile = "scombine_file.list"
    weightfile = "scombine_weight.list"
    sf = open(scombfile, "w")
    counter = 1
    copylist = []
    for i in splist:
        inputfile = os.path.basename(i)
        copyfile = inputfile.replace(".fits", "-cp%d.fits" % counter)
        copylist.append(copyfile)
        if not os.path.exists(outputdir + "/" + copyfile):
            shutil.copy(i, outputdir + "/" + copyfile)
        sf.write("%s\n" % copyfile)
        counter += 1
    sf.close()
    wf = open(weightfile, "w")
    for i in weightlist:
        wf.write("%.3e\n" % i)
    wf.close()

    iraf.scombine("@%s" % scombfile, outputfile, combine="average", weight="@%s" % weightfile, group="all",
                  reject="avsigclip", lsigma=3., hsigma=3.)

    for i in copylist:
        os.remove("./%s" % os.path.basename(i))

    os.chdir(curdir)


def obtainWeight(combineID):
    conn, cur = openproject()
    cur.execute("SELECT datasetID,weight from combinedataset where combineID='%s';" % combineID)
    rows = cur.fetchall()
    datasetID = [i[0] for i in rows]
    weight = [float(i[1]) for i in rows]

    return datasetID, weight


def openTelluricSpectra(picklepath):
    telspecpcl = open(picklepath + "telluric_spectra.pickle", "rb")
    telspnpz = pickle.load(telspecpcl)
    tel_obs = telspnpz["obs"]
    tel_model = telspnpz["model"]

    # cur.execute("SELECT datasetID,weight from combinedataset where combineID='%s';" % combineID)
    # rows = cur.fetchall()
    # datasetID = [i[0] for i in rows]
    # weight = [float(i[1]) for i in rows]
    #
    telspecpcl.close()
    # conn.close()

    return tel_obs, tel_model


def saveTelluricSpectra(combinepath, obs, model):
    telspecpcl = open(combinepath + "telluric_spectra.pickle", "wb")
    spec = {}
    spec["obs"] = obs
    spec["model"] = model
    pickle.dump(spec, telspecpcl)
    telspecpcl.close()


def obtainCombinePath(combineID):
    conn, cur = openproject()
    cur.execute("SELECT combinepath from combinesummary where combineID='%s';" % combineID)
    rows = cur.fetchall()

    if rows == []:
        print("combineID='%s' was not found." % combineID)
        return None
    else:
        combpath = rows[0][0]

    conn.close()

    return combpath


def obtainCombineFileList(combineID, normflag=False):
    conn, cur = openproject()

    cur.execute("SELECT telluricflag, combineflag from combinesummary where combineID='%s';" % combineID)
    rows = cur.fetchall()
    telflag = rows[0][0]
    combineflag = rows[0][1]

    cur.execute("SELECT combineID,echelleorder,combinefilepath, lambdamin, lambdamax from combinedspectrum where combineID='%s';" % combineID)
    rows = cur.fetchall()

    filepath = {}
    lammin = {}
    lammax = {}
    if rows == []:
        print("combineID='%s' was not found." % combineID)
        return None
    else:
        m = [int(i[1]) for i in rows]
        if normflag:
            for i in range(len(m)):
                if telflag == 0 and combineflag == 0:
                    filepath[m[i]] = rows[i][2].replace("AIR_wscor", "AIR_norm_wscor")
                else:
                    filepath[m[i]] = rows[i][2].rstrip("fits").rstrip(".") + "_norm.fits"
                lammin[m[i]] = float(rows[i][3])
                lammax[m[i]] = float(rows[i][4])
        else:
            for i in range(len(m)):
                filepath[m[i]] = rows[i][2]
                lammin[m[i]] = float(rows[i][3])
                lammax[m[i]] = float(rows[i][4])

    conn.close()
    m.sort()

    return [m, filepath, lammin, lammax]


if __name__ == "__main__":
    pathcore = "/Users/hamano/DIB_analysis/DIB_pipeline_dir/"
    plotsig = 5
    ite = 20
    vacorair = "AIR"

    conn, cur = openproject()

    idsinput = sys.argv[1:]
    combnum = len(idsinput)

    path, objectid_list, objname_list, mode_list, obsdate_list, telluricflags = [], [], [], [], [], []
    for i in idsinput:
        cur.execute(
            "select x.telluricPath, y.objectID, z.objectname, y.mode, y.obsdate from telluriccorrection as x JOIN datareduction as y on x.pipelineIDobj=y.pipelineID JOIN object as z on y.objectid=z.objectid where x.telluricID='%s';" % i)
        rows = cur.fetchall()
        if rows == []:
            cur.execute(
                "select y.path, y.objectID, z.objectname, y.mode, y.obsdate from datareduction as y JOIN object as z on y.objectid=z.objectid where y.pipelineID='%s';" % i)
            rows = cur.fetchall()
            if rows == []:
                print("%s is not registered." % i)
                sys.exit()
            elif len(rows) == 1:
                path.append(rows[0][0])
                objectid_list.append(rows[0][1])
                objname_list.append(rows[0][2])
                mode_list.append(rows[0][3])
                obsdate_list.append(rows[0][4])
                telluricflags.append(False)
        else:
            path.append(rows[0][0])
            objectid_list.append(rows[0][1])
            objname_list.append(rows[0][2])
            mode_list.append(rows[0][3])
            obsdate_list.append(rows[0][4])
            telluricflags.append(True)

    telconverter = {True: "T", False: "R"}
    boolconv = {True: "True", False: "False"}

    v_helio = [read_rvfile_ID(i) for i in idsinput]

    if all(telluricflags):
        print("Data type: Telluric corrected data")
        telflag = True
    elif not any(telluricflags):
        print("Data type: Telluric not-corrected data")
        telflag = False
    else:
        for i in range(combnum):
            print(idsinput[i], ":", telconverter[telluricflags[i]])
            print("The data type are not the same.")
            sys.exit()

    objectid = objectid_list[0]
    mode = mode_list[0]
    objname = objname_list[0].replace("*", "").lstrip(" ").replace("   ", " ").replace("  ", " ").replace(" ", "_")
    for i in range(combnum):
        if objectid_list[i] != objectid:
            print("The object ID error (1: %s, %d: %s)" % (objname_list[0], i + 1, objname_list[i]))
            sys.exit()
        if mode_list[i] != mode:
            print("The mode error (1: %s, %d: %s)" % (mode_list[0], i + 1, mode_list[i]))
            sys.exit()

    cur.execute("select combineNumber from combinesummary where objectID=%d;" % objectid)
    rows = cur.fetchall()

    if rows == []:
        curID = 1
    else:
        combIDs = [int(i[0]) for i in rows]
        curID = max(combIDs) + 1

    combineID = "%s_o%d_c%d%s" % (objname, objectid, curID, telconverter[telflag])

    frameinput = []
    print("No input: frame is set as 'sum'.")
    print("Frame number (int): frame is set as 'No#'.")
    for i in range(combnum):
        ans = input("Frame for %s: " % idsinput[i])
        if ans == "" or ans == "sum":
            frameinput.append("sum")
        else:
            frameinput.append("NO%d" % int(ans))

    if combnum > 1:
        print("\nSet weights.")
        weight = []
        for i in range(combnum):
            weight.append(float(input("Weight for %s: " % idsinput[i])))
        wmaxid = weight.index(max(weight))

        combinepath = pathcore + objname + "/" + combineID + "/"
        outputdir = combinepath + "spectrum_files/"
        cur.execute(
            "INSERT IGNORE INTO combinesummary (combineID,objectID,combinepath,mode,telluricflag,combineflag,combineNumber) VALUES ('%s',%d,'%s','%s',%s,%s,%d);" % (
                combineID, objectid, combinepath, mode, boolconv[telflag], "True", curID))

        spfiles_list = []
        for i in range(combnum):
            if telflag:
                orders, frameset, spfiles = obtainObjTelFileList(path[i], fluxornorm="flux", vacorair=vacorair)
            else:
                orders, frameset, spfiles = obtainWscorFileList(path[i], fluxornorm="flux", vacorair=vacorair)
            spfiles_list.append(spfiles)
            if not frameinput[i] in frameset:
                print("Frame '%s' does not exist for ID=%s." % (frameinput[i], idsinput[i]))
                sys.exit()
            cur.execute(
                "INSERT IGNORE INTO combinedataset (combineID,datasetID,weight,frame) VALUES ('%s','%s',%.3f,'%s');" % (
                    combineID, idsinput[i], weight[i], frameinput[i]))

        tel_obs = {}
        if telflag:
            telfiles = {}
            for i in idsinput:
                cur.execute(
                    "select y.path, x.pipelineIDtel, x.advanced from telluriccorrection as x join datareduction as y on x.pipelineIDobj=y.pipelineID where x.telluricID='%s';" %
                    i)
                rows = cur.fetchall()
                plpath = rows[0][0]
                telID = rows[0][1]
                advflag = int(rows[0][2])
                telfiles[i] = obtainWscorTelluricFileList(plpath, telID, advflag, vacorair=vacorair)
                tel_obs[i] = {}
                for m in orders:
                    telspx, telspy, _, _, _ = openspecfits(telfiles[i][orders.index(m)])
                    tel_obs[i][m] = [telspx, telspy]

        tel_model = {}
        for i in range(len(idsinput)):
            telwav_model1, telflux_model1, _ = modelnpzopenHelio(v_helio[i], vacorair="AIR")
            tel_model[idsinput[i]] = [telwav_model1, telflux_model1]

        os.makedirs(outputdir)
        pp = PdfPages(outputdir + combineID + ".pdf")
        saveTelluricSpectra(combinepath, tel_obs, tel_model)

        fig = plt.figure(figsize=(20, 4 + 2 * combnum))
        colorsp = ["k" for i in range(combnum + 1)]
        colorsp[-1] = "b"
        conn.commit()

        for m in orders:
            spxlist, spylist, spylist_norm = [], [], []
            lammin, lammax, ymin, ymax, yrange, ymed, ystd = [], [], [], [], [], [], []
            combinesplist = []
            for i in range(combnum):
                combinesplist.append(spfiles_list[i][frameinput[i]][m])
                spx, spy, _, _, _ = openspecfits(combinesplist[i])
                spxlist.append(spx)
                spylist.append(spy)

            outputfile = outputdir + combineID + "_m%d.fits" % m
            pycombine(combinesplist, outputfile, weight)
            spx, spy, _, _, _ = openspecfits(outputfile)
            spxlist.append(spx)
            spylist.append(spy)
            lambdamin = min(spx)
            lambdamax = max(spx)

            texts = []
            for i in range(combnum + 1):
                if i != combnum:
                    texts.append(" %s\n  (weight=%.3f)" % (idsinput[i].replace("--", "\n"), weight[i]))
                else:
                    texts.append(" %s\n" % (combineID))

            if telflag:
                MultiSpecPlotter(spxlist, spylist, pp, colorsp, texts, "%s (m=%d)" % (combineID, m), fs=6.,
                                 obsdates=obsdate_list, order=m, telx=tel_obs[idsinput[wmaxid]][m][0], v_helio=v_helio,
                                 tely=tel_obs[idsinput[wmaxid]][m][1], telflag=True)
            else:
                MultiSpecPlotter(spxlist, spylist, pp, colorsp, texts, "%s (m=%d)" % (combineID, m), fs=6.,
                                 obsdates=obsdate_list, order=m, telx=tel_model[idsinput[wmaxid]][0], v_helio=v_helio,
                                 tely=tel_model[idsinput[wmaxid]][1], telflag=True)

            cur.execute(
                "INSERT IGNORE INTO combinedspectrum (combineID, echelleorder, combinefilepath, lambdamin, lambdamax) VALUES ('%s', '%s', '%s', %.4f, %.4f);"
                % (combineID, m, outputfile, lambdamin, lambdamax)
            )
        pp.close()

        conn.commit()

    else:
        tel_obs = {}
        telfiles = {}
        if telflag:
            orders, frameset, spfiles = obtainObjTelFileList(path[0], fluxornorm="flux", vacorair=vacorair)
            cur.execute(
                "select y.path, x.pipelineIDtel, x.advanced from telluriccorrection as x join datareduction as y on x.pipelineIDobj=y.pipelineID where x.telluricID='%s';" %
                idsinput[0])
            rows = cur.fetchall()
            plpath = rows[0][0]
            telID = rows[0][1]
            advflag = int(rows[0][2])
            telfiles[idsinput[0]] = obtainWscorTelluricFileList(plpath, telID, advflag, vacorair=vacorair)
            tel_obs[idsinput[0]] = {}
            for m in orders:
                telspx, telspy, _, _, _ = openspecfits(telfiles[idsinput[0]][orders.index(m)])
                tel_obs[idsinput[0]][m] = [telspx, telspy]
        else:
            orders, frameset, spfiles = obtainWscorFileList(path[0], fluxornorm="flux", vacorair=vacorair)

        combinepath = os.path.dirname(spfiles[frameinput[0]][orders[0]]).rstrip(spfiles[frameinput[0]][orders[0]].split("/")[-2])

        tel_model = {}
        telwav_model1, telflux_model1, _ = modelnpzopenHelio(v_helio[0], vacorair="AIR")
        tel_model[idsinput[0]] = [telwav_model1, telflux_model1]
        saveTelluricSpectra(combinepath, tel_obs, tel_model)

        if not frameinput[0] in frameset:
            print("Frame '%s' does not exist for ID=%s." % (frameinput[0], idsinput[0]))
            sys.exit()
        cur.execute(
            "INSERT IGNORE INTO combinesummary (combineID,objectID,combinepath,mode,telluricflag,combineflag,combineNumber) VALUES ('%s',%d,'%s','%s',%s,%s,%d);" % (
                combineID, objectid, combinepath, mode, boolconv[telflag], "False", curID))
        cur.execute(
            "INSERT IGNORE INTO combinedataset (combineID,datasetID,weight,frame) VALUES ('%s','%s',%.3f,'%s');" % (
                combineID, idsinput[0], 1.0, frameinput[0]))
        for m in orders:
            spx, _, _, _, _ = openspecfits(spfiles[frameinput[0]][m])
            lambdamin = min(spx)
            lambdamax = max(spx)
            cur.execute(
                "INSERT IGNORE INTO combinedspectrum (combineID, echelleorder, combinefilepath, lambdamin, lambdamax) VALUES ('%s', '%s', '%s', %.4f, %.4f);"
                % (combineID, m, spfiles[frameinput[0]][m], lambdamin, lambdamax)
            )
        conn.commit()

    print("combine ID: %s" % combineID)

    conn.close()
