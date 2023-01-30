# -*- coding:utf-8 -*-

import sys, os, shutil
from DIBanalysis import sampleRegion
# from pyraf import iraf
from combine_MySQL import openTelluricSpectra, obtainCombinePath, pycombine, saveTelluricSpectra, obtainWeight, \
    obtainCombineFileList
from open_mysql_project import openproject
from DIBanalysis import maskRegion_empty, maskRegion_tel, continuum

def normalizeCombinedspec(inputfits, outputfits, maskbool):
    samplereg = sampleRegion(inputfits, maskbool, minsep=5.)
    # print("Sample region for %s:" % os.path.basename(inputfits))
    # print("continuum %s %s sample='%s' func='legendre' order=7" % (inputfits.split("/")[-1], outputfits.split("/")[-1], samplereg))
    curdir = os.getcwd()
    outputdir = os.path.dirname(outputfits)

    os.chdir(outputdir)
    continuum(inputfits.split("/")[-1], "temporaly_norm.fits", window=samplereg, order=7.)
    shutil.move("temporaly_norm.fits", outputfits.split("/")[-1])
    print(outputfits.split("/")[-1])
    os.chdir(curdir)

if __name__ == "__main__":

    combineID = sys.argv[1]
    combinepath = obtainCombinePath(combineID)
    tel_obs, tel_model = openTelluricSpectra(combinepath)
    [orders, spfilepath, lammin, lammax] = obtainCombineFileList(combineID)
    datasetID, weight = obtainWeight(combineID)
    curdir = os.getcwd()

    print(combineID)
    print("cd %s" % combinepath)


    if tel_obs != {}:
        telspec = tel_obs
        telmode = "obs"
    else:
        telspec = tel_model
        telmode = "model"

    maskbool_tel_fornorm = {}
    for m in orders:
        output = spfilepath[m].rstrip("fits").rstrip(".") + "_norm.fits"
        if os.path.exists(output):
            os.remove(output)
        if not os.path.exists(output):
            maskbool_tel_fornorm[m] = maskRegion_tel(spfilepath[m], telspec, m, telmode, thres=0.1, thres_upside=1.5)
            # maskbool_tel_fornorm[m] = maskRegion_empty(spfilepath[m])
            # if telmode == "obs":
            #     for i in telspec.keys():
            #         maskbool_tel_fornorm[m] *= maskRegion_tel(spfilepath[m], telspec[i][m], thres=0.1, thres_upside=1.5)
            # elif telmode == "model":
            #     for i in telspec.keys():
            #         maskbool_tel_fornorm[m] *= maskRegion_tel(spfilepath[m], telspec[i], thres=0.1)

            normalizeCombinedspec(spfilepath[m], output, maskbool_tel_fornorm[m])