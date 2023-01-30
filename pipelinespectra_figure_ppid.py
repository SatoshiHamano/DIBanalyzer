#!/usr/bin/env python
# -*- coding:utf-8 -*-

from combine_MySQL import obtainCombineFileList, obtainCombinePath, openTelluricSpectra, obtainWeight
from spectra_plotter import MultiSpecPlotter
from Spec1Dtools import openspecfits
import sys
import numpy
from dopcor_vcorr import read_rvfile_ID
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse
from open_mysql_project import openproject


def pipelineSpecFig(outputpdf, pipelinepath, xmax=13400.):
    splist, lammin, lammax = [], [], []
    orders = set([])
    rtlist = []
    for i in range(len(pipelineID)):
        sqlresult = obtainCombineFileList(combineID[i], normflag=normalizedflag)
        if sqlresult != None:
            orders = orders | set(sqlresult[0])
            splist.append(sqlresult[1])
            lammin.append(sqlresult[2])
            lammax.append(sqlresult[3])
            rtlist.append(righttexts[i])
        else:
            print(i, " is not found.")

    orderslist = list(orders)
    orderslist.sort()

    pp = PdfPages(outputpdf)
    plt.figure(figsize=(15, 20))

    for m in orderslist:
        spxlist, spylist = [], []
        if telmode == "obs":
            telx, tely = telspec[datasetID[maxid]][m]
        elif telmode == "model":
            telx, tely = telspec[datasetID[maxid]]
        for i in range(len(splist)):
            spx, spy, _, _, _ = openspecfits(splist[i][m])
            spxlist.append(spx[spx<xmax])
            spylist.append(spy[spx<xmax] / numpy.median(spy[spx<xmax]))
        colors = ["k" for j in spxlist]
        MultiSpecPlotter(spxlist, spylist, pp, colors, righttexts, "m=%d" % m, order=m, obsdates=obsdate, telx=telx,
                         tely=tely, v_helio=vhelio, DIBgrid=True, plotsig=5.)

    pp.close()


if __name__ == "__main__":
    # file mode
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("readfile", type=str, help="readfile")
    parser.add_argument("outputpdf", type=str, help="outputpdf")
    parser.add_argument("-i", "--repid", type=int, default=0, help="Representative dataset id")
    parser.add_argument("-n", "--normalized", action="store_true", help="normalized spec")

    args = parser.parse_args()

    readfile = args.readfile
    outputpdf = args.outputpdf
    rep_id = args.repid
    normalizedflag = args.normalized

    rf = open(readfile, "r")
    rl = rf.readlines()
    rf.close()

    combineID = [i.split("|")[0] for i in rl]
    righttexts = [i.split("|")[1].replace("\\n", "\n") for i in rl]

    ### object id mode
    # parser = argparse.ArgumentParser()
    # parser.add_argument("objid", type=int, help="object ID")
    # parser.add_argument("outputpdf", type=str, help="outputpdf")
    # parser.add_argument("-i", "--repid", type=int, default=0, help="Representative dataset id")
    # parser.add_argument("-n", "--normalized", action="store_true", help="normalized spec")
    #
    # args = parser.parse_args()
    #
    # objid = args.objid
    # outputpdf = args.outputpdf
    # rep_id = args.repid
    # normalizedflag = args.normalized
    #
    # conn, cur = openproject()
    # cur.execute("SELECT combineID,mode,telluricflag,combineflag from combinesummary where objectID=%d;" % objid)
    # rows = cur.fetchall()
    # if rows == []:
    #     print("objectID='%d' was not found." % objid)
    #     sys.exit()
    # else:
    #     combineID = [i[0] for i in rows]
    #     mode = [i[1] for i in rows]
    #     telflag = [i[2] for i in rows]
    #     combineflag = [i[3] for i in rows]
    #
    # righttexts = ["%s\n%s\n(T,C)=(%d,%d)\n" % (combineID[i], mode[i], telflag[i], combineflag[i]) for i in range(len(combineID))]

    # common part
    combineSpecFig(outputpdf, combineID, righttexts, rep_id)
