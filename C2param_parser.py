# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import configparser
from C2_parameters import *

C2params = moleclinelist("C2")
C2settingFile = "C2setting.txt"


def alternativequestion(question, ans1, ans2):
    flagans = 0
    while flagans == 0:
        flag = input(question)
        if flag == ans1 or flag == ans2:
            flagans += 1
        else:
            print("please answer to the question with \"%s\" or \"%s\" !!!" % (ans1, ans2))

    return flag


def WriteSetting(settingfile):
    setf = open(C2settingFile, "a")
    setf.write(settingfile)
    setf.write("\n")
    setf.close()


def OverwriteConfig(config):
    setf = open(C2settingFile, "r")
    setl = setf.readlines()
    setf.close()

    for i in range(len(setl)):
        configpath = setl[i].split("#")[0].split()[0]

    with open(configpath, 'w') as configfile:
        config.write(configfile)


def ReadConfig(band="AX00"):
    if not os.path.exists(C2settingFile):
        pfile = input("Parameter setting file: ")
        WriteSetting(pfile)
        InitializeParameter(pfile, band)

    setf = open(C2settingFile, "r")
    setl = setf.readlines()
    setf.close()

    for i in range(len(setl)):
        configpath = setl[i].split("#")[0].split()[0]

    if not os.path.exists(configpath):
        InitializeParameter(configpath, band)

    config = configparser.ConfigParser()
    config.read(configpath)

    return config


def InitializeParameter(pfile, band):
    C2lines = C2params.bandlist[band].lines

    config = configparser.ConfigParser()

    flagfloat = False
    while not flagfloat:
        resolution = input("Resolution: ")
        try:
            resolution = float(resolution)
            flagfloat = True
        except:
            print("%s is invalid input." % resolution)

    config["Spectrum"] = {"Target": "", "Telluric": "", "Resolution": resolution}

    config["Cloud parameters"] = {"n_comp": 0}

    config["C2 parameters"] = {"band": "", "detection": ""}
    for i in C2lines:
        config["C2 parameters"][i.transition] = "on"

    with open(pfile, 'w') as configfile:
        config.write(configfile)


def VelocityReader(config):
    n_comp = int(config["Cloud parameters"]["n_comp"])
    vel_comp = []
    for i in range(n_comp):
        vel = float(config["Cloud parameters"]["vel%d" % (i + 1)])
        vel_comp.append(vel)

    return n_comp, vel_comp


def CloudParameterReader(config):
    n_comp = int(config["Cloud parameters"]["n_comp"])
    vel_comp = []
    velerr_comp = []
    b_comp = []
    berr_comp = []
    for i in range(n_comp):
        vel = float(config["Cloud parameters"]["vel%d" % (i + 1)])
        velerr = float(config["Cloud parameters"]["vel%d_err" % (i + 1)])
        b = float(config["Cloud parameters"]["b%d" % (i + 1)])
        berr = float(config["Cloud parameters"]["b%d_err" % (i + 1)])
        vel_comp.append(vel)
        velerr_comp.append(velerr)
        b_comp.append(b)
        berr_comp.append(berr)

    return n_comp, vel_comp, velerr_comp, b_comp, berr_comp


def RegionReader(config):
    n_region = int(config["Spectrum"]["region number"])
    region_ranges = []
    target_region = []
    telluric_region = []
    target_region_normfits = []
    target_region_normtxt = []
    region_snr = []
    for i in range(n_region):
        reg_range = config["Spectrum"]["region%d range" % (i + 1)]
        regmin, regmax = float(reg_range.split(",")[0]), float(reg_range.split(",")[1])
        region_ranges.append([regmin, regmax])
        if config.has_option("Spectrum", "region%d target" % (i + 1)):
            target = config["Spectrum"]["region%d target" % (i + 1)]
            target_region.append(target)
        if config.has_option("Spectrum", "region%d telluric" % (i + 1)):
            telluric = config["Spectrum"]["region%d telluric" % (i + 1)]
            telluric_region.append(telluric)
        if config.has_option("Spectrum", "region%d target normfits" % (i + 1)):
            target_normfits = config["Spectrum"]["region%d target normfits" % (i + 1)]
            target_region_normfits.append(target_normfits)
        if config.has_option("Spectrum", "region%d target normtxt" % (i + 1)):
            target_normtxt = config["Spectrum"]["region%d target normtxt" % (i + 1)]
            target_region_normtxt.append(target_normtxt)
        if config.has_option("Spectrum", "region%d snr" % (i + 1)):
            snr = config["Spectrum"]["region%d snr" % (i + 1)]
            try:
                snr = float(snr)
                region_snr.append(snr)
            except:
                pass

    return n_region, region_ranges, target_region, telluric_region, target_region_normfits, target_region_normtxt, region_snr


def ColumnDensityReader(config):
    J = []
    logN = []
    logN_err = []
    minJ = int(config["C2 parameters"]["min J"])
    maxJ = int(config["C2 parameters"]["max J"])
    n_comp = int(config["Cloud parameters"]["n_comp"])

    exflag = False
    for j in range(minJ, maxJ + 1, 2):
        N = []
        Nerrone = []
        for k in range(n_comp):
            if config.has_option("Cloud parameters", "logn_j%d_comp%d" % (j, (k + 1))):
                N.append(float(config["Cloud parameters"]["logn_j%d_comp%d" % (j, (k + 1))]))
                Nerrone.append(float(config["Cloud parameters"]["logn_j%d_comp%d_err" % (j, (k + 1))]))
                exflag = True
        if exflag:
            J.append(j)
            logN.append(N)
            logN_err.append(Nerrone)
        exflag = False

    return J, logN, logN_err


def EWReader(config):
    n_comp = int(config["Cloud parameters"]["n_comp"])
    band = config["C2 parameters"]["band"]
    C2lines = C2params.bandlist[band].lines
    C2bool = [config.getboolean("C2 parameters", c2line.transition) for c2line in C2lines]

    lines = []
    ew = []
    ewerr = []

    for i in range(len(C2lines)):
        if C2bool[i]:
            lines.append(C2lines[i])
            ew.append([])
            ewerr.append([])
            for j in range(len(n_comp)):
                if config.has_option("Cloud parameters", "ew_%s_comp%d" % (C2lines[i].transition, j + 1)):
                    ew[i].append(float(config["Cloud parameters"]["ew_%s_comp%d" % (C2lines[i].transition, j + 1)]))
                    ewerr[i].append(
                        float(config["Cloud parameters"]["ew_%s_comp%d_err" % (C2lines[i].transition, j + 1)]))
                else:
                    ew[i].append(0.)
                    ewerr[i].append(0.)

    return lines, ew, ewerr


if __name__ == "__main__":
    InitializeParameter("test.txt", "AX00")
