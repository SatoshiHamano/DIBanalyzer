# -*- coding:utf-8 -*-

import sys
import mysql.connector
from urllib.parse import urlparse
import glob
import os

def obtainObjTelFileList(telpath, frame="all", fluxornorm="norm", vacorair="VAC", wscor=True, helio=True):
    alllist = glob.glob(telpath + "*")
    framelist = []
    for i in alllist:
        if os.path.isdir(i):
            if i.find("pycache") == -1 and i.find("DIBs") == -1:
                framelist.append(i)
    framelist.sort()
    frameset = [framelist[i].split("/")[-1] for i in range(len(framelist))]
    filelist = {}

    filename = vacorair + "_norm_tel" if fluxornorm == "norm" else vacorair + "_tel"
    if wscor:
        filename += "_wscor"
    if helio:
        filename += "_helio"
    filename += ".fits"

    orderlist = []

    if frame == "all":
        for i in range(len(framelist)):
            filelist[frameset[i]] = {}
            files = glob.glob(framelist[i] + "/*%s" % filename)
            orders = [int(files[j].split("/")[-1].split("%s_m" % frameset[i])[1].split("_")[0]) for j in range(len(files))]
            for j in range(len(orders)):
                filelist[frameset[i]][orders[j]] = files[j]
                orderlist.append(orders[j])
    else:
        for i in range(len(framelist)):
            if framelist[i].split("/")[-1] == "sum":
                filelist[frameset[i]] = {}
                files = glob.glob(framelist[i] + "/*%s" % filename)
                orders = [int(files[j].split("/")[-1].split("%s_m" % frameset[i])[1].split("_")[0]) for j in range(len(files))]
                for j in range(len(orders)):
                    filelist[frameset[i]][orders[j]] = files[j]
                    orderlist.append(orders[j])

    orderlist = list(set(orderlist))
    orderlist.sort()

    return orderlist, frameset, filelist

def obtainWscorTelluricFileList(plpath, telID, advflag, vacorair="VAC", helio=True):
    telpath = plpath + telID + "_adv_wscor/" if advflag else plpath + telID + "_wscor/"
    if vacorair == "VAC":
        if helio:
            files = glob.glob("%s*VAC*wscor_helio.fits" % telpath)
        else:
            files = glob.glob("%s*VAC*wscor.fits" % telpath)
    else:
        if helio:
            files = glob.glob("%s*AIR*wscor_helio.fits" % telpath)
        else:
            files = glob.glob("%s*AIR*wscor.fits" % telpath)

    files.sort()

    return files


def obtainWscorFileList(plpath, fluxornorm="norm", frame="all", fsr="1.30", vacorair="VAC", helio=True):
    alllist = glob.glob(plpath + "waveshift_correct/*")
    framelist = []
    for i in alllist:
        if os.path.isdir(i):
            framelist.append(i)
    framelist.sort()
    frameset = [framelist[i].split("/")[-1] for i in range(len(framelist))]
    filelist = {}

    filename = "fsr" + fsr +  "_" + vacorair + "_norm_wscor" if fluxornorm=="norm" else "fsr" + fsr +  "_" + vacorair + "_wscor"
    if helio:
        filename += "_helio"
    filename += ".fits"

    orderlist = []

    if frame == "all":
        for i in range(len(framelist)):
            filelist[frameset[i]] = {}
            files = glob.glob(framelist[i] + "/fsr%s/*%s" % (fsr, filename))
            orders = [int(files[j].split("/")[-1].split("%s_m" % frameset[i])[1].split("_")[0]) for j in range(len(files))]
            for j in range(len(orders)):
                filelist[frameset[i]][orders[j]] = files[j]
                orderlist.append(orders[j])
    else:
        for i in range(len(framelist)):
            if framelist[i].split("/")[-1] == "sum":
                filelist[frameset[i]] = {}
                files = glob.glob(framelist[i] + "/fsr%s/*%s" % (fsr, filename))
                orders = [int(files[j].split("/")[-1].split("%s_m" % frameset[i])[1].split("_")[0]) for j in range(len(files))]
                for j in range(len(orders)):
                    filelist[frameset[i]][orders[j]] = files[j]
                    orderlist.append(orders[j])

    orderlist = list(set(orderlist))
    orderlist.sort()

    return orderlist, frameset, filelist


#def obtain