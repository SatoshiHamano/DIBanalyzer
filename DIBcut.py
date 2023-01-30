#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sys
import glob
import numpy
import time
import copy
import os
import scipy.constants
from obtain_filelist import *
from add_lineDIB_mysql import *
import argparse
from open_mysql_project import openproject
from combine_MySQL import obtainCombinePath, obtainCombineFileList
from Spec1Dtools import cutSpectrum

def DIBcut(spfile, output, DIBlam, vel_c, vel_hw=1000.):
    # inputcopy = "./work_place/" + spfile.split("/")[-1].replace(".fits", "-cp.fits")
    # outputcopy = "./work_place/" + output.split("/")[-1].replace(".fits", "-cp.fits")
    # shutil.copyfile(spfile, inputcopy)
    #
    # v_light = scipy.constants.c * 1.e-3
    # minlam = DIBlam * (1. + (vel_c - vel_hw) / v_light)
    # maxlam = DIBlam * (1. + (vel_c + vel_hw) / v_light)
    # iraf.scopy(inputcopy, outputcopy, w1=minlam, w2=maxlam)
    #
    # shutil.copyfile(outputcopy, output)
    # os.remove(inputcopy)
    # os.remove(outputcopy)
    v_light = scipy.constants.c * 1.e-3
    minlam = DIBlam * (1. + (vel_c - vel_hw) / v_light)
    maxlam = DIBlam * (1. + (vel_c + vel_hw) / v_light)
    cutSpectrum(spfile, output, minlam, maxlam)


def writeDIBcut_info(wfname, text):
    wf = open(wfname, "a")
    wf.write(time.ctime())
    wf.write(text)
    wf.close()


def readDIBcut_vel(DIBID, combineID):
    combpath = obtainCombinePath(combineID)
    textfile = glob.glob(combpath + "DIBs/DIBID%d/DIBID%d_lambda*.txt" % (DIBID, DIBID))

    rf = open(textfile[0], "r")
    rl = rf.readlines()
    rf.close()
    vel = 0.
    for i in rl:
        if i.find("velocity: ") != -1:
            vel = float(i.split("velocity: ")[1].split()[0])

    return vel


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("combineid", type=str, help="combineid")
    parser.add_argument("-d", "--dibid", type=int, default=0, help="DIB ID")
    parser.add_argument("-v", "--velocity", type=float, default=0., help="central velocity")
    parser.add_argument("-c", "--category", type=str, default="all", help="DIB category")

    args = parser.parse_args()

    combineid = args.combineid
    dibid = args.dibid
    velocity = args.velocity
    category = args.category

    atranmodel = "atran.smo.11513_R28000_0ft.npz"
    vacorair = "AIR"

    conn, cur = openproject()
    cur.execute(
        "SELECT echelleorder,combinefilepath,lambdamin,lambdamax from combinedspectrum where combineID='%s';" % combineid)
    rows = cur.fetchall()
    if rows == []:
        print("combineID=%s was not found." % combineid)
        sys.exit()
    else:
        orders = [int(i[0]) for i in rows]
        spfilepath = [i[1] for i in rows]
        lammin = [float(i[2]) for i in rows]
        lammax = [float(i[3]) for i in rows]

    cur.execute(
        "SELECT x.telluricflag, y.datasetID, y.weight from combinesummary as x JOIN combinedataset as y using (combineID) where x.combineID='%s';" % combineid)
    rows = cur.fetchall()
    weight = [float(i[2]) for i in rows]
    maxid = weight.index(max(weight))
    telflag = rows[maxid][0]
    datasetID = rows[maxid][1]

    if telflag:
        cur.execute(
            "select y.path, x.pipelineIDtel, x.advanced from telluriccorrection as x join datareduction as y on x.pipelineIDobj=y.pipelineID where x.telluricID='%s';" %
            datasetID)
        rows = cur.fetchall()
        plpath = rows[0][0]
        telID = rows[0][1]
        advflag = int(rows[0][2])
        telfiles = obtainWscorTelluricFileList(plpath, telID, advflag, vacorair=vacorair)
    else:
        atran_npz = numpy.load(atranmodel)
        telwav = atran_npz["wav"]
        telflux = atran_npz["flux"]
        id = atran_npz["id"]

    conn.close()

    for i in range(len(orders)):
        if dibid == 0:
            dibinfo = GetDIBList(lammin[i], lammax[i], category=category)
            if dibinfo == None:
                print("No DIBs were found in the database.")
                continue
            else:
                [DIBIDs, wav_air, reference, DIBcate, fwhm, wavenumber, comment] = dibinfo
        else:
            dibinfo = GetDIBLine(dibid)
            if dibinfo == None:
                print("No DIBs were found in the database.")
                continue
            else:
                [DIBIDs, wav_air, reference, DIBcate, fwhm, wavenumber, comment] = dibinfo
                DIBIDs, wav_air, reference, DIBcate, fwhm, wavenumber, comment = [DIBIDs], [wav_air], [reference], [
                    DIBcate], [fwhm], [wavenumber], [comment]
        outputdir = os.path.dirname(spfilepath[i]).rstrip(spfilepath[i].split("/")[-2]) + "DIBs/"
        # sprelpath = "./%s/%s" % (spfilepath[i].split("/")[-2], spfilepath[i].split("/")[-1])
        for j in range(len(DIBIDs)):
            if lammin[i] < wav_air[j] < lammax[i]:
                # outputfile = sprelpath.split("/")[-1].rstrip("fits").rstrip(".") + "_DIBID%d.fits" % (DIBIDs[j])
                # outputfile = spfilepath[i].split("/")[-1].rstrip("fits").rstrip(".") + "_DIBID%d.fits" % (DIBIDs[j])
                outputfile = combineid + "_m%d_DIBID%d.fits" % (orders[i], DIBIDs[j])
                DIBdir = outputdir + "DIBID%d/" % (DIBIDs[j])
                # print(sprelpath, outputfile)
                if not os.path.exists(DIBdir):
                    os.makedirs(DIBdir)
                if not os.path.exists(DIBdir + outputfile):
                    DIBcut(spfilepath[i], DIBdir + outputfile, wav_air[j], velocity)
                    # DIBcut(sprelpath, DIBdir + outputfile, wav_air[j], velocity)
                    # print("DIBcut('%s', '%s', %.5f, %.2f)" % (spfilepath[i], DIBdir + outputfile, wav_air[j], velocity))
                    writeDIBcut_info("%sDIBID%d_lambda%.1f.txt" % (DIBdir, DIBIDs[j], wav_air[j]),
                                     "\n%s\nvelocity: %.1f\n\n" % (DIBdir + outputfile, velocity))
