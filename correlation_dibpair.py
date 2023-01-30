#!/usr/bin/env python
# -*- coding:utf-8 -*-
import math

from open_mysql_project import openproject
import argparse
import numpy
import scipy.stats
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
from correlation_ebv import styleTxtParsar, confband
from add_lineDIB_mysql import GetDIBLine


def DIBcorrelation(ax, ew1, ewerr1, objid1, DIBwav1, ew2, ewerr2, objid2, DIBwav2, stylefile="INDEF",
                   colordef="k", fmtdef="o", msdef=3., alpha=0.95):
    objidcom = numpy.intersect1d(objid1, objid2)
    ew1com = numpy.array([ew1[i] for i in objidcom])
    ewerr1com = numpy.array([ewerr1[i] for i in objidcom])
    ew2com = numpy.array([ew2[i] for i in objidcom])
    ewerr2com = numpy.array([ewerr2[i] for i in objidcom])

    detection12 = (ew1com != 0.) & (ew2com != 0.)
    detection1 = (ew1com != 0.) & (ew2com == 0.)
    detection2 = (ew1com == 0.) & (ew2com != 0.)
    uplimit = (ew1com == 0.) & (ew2com == 0.)

    if stylefile != "INDEF":
        objids, ps, mark, color, legend = styleTxtParsar(stylefile)
        objidsall = []
        for i in objids:
            objidsall += i
        objidsallset = list(set(objidsall))
        if len(objidsallset) != len(objidsall):
            print("Error: overlapped object IDs.")
            sys.exit()
        for i in range(len(objids)):
            req = numpy.array([j in objids[i] for j in objidcom])
            if any(req):
                det12req = numpy.logical_and(req, detection12)
                det1req = numpy.logical_and(req, detection1)
                det2req = numpy.logical_and(req, detection2)
                upreq = numpy.logical_and(req, uplimit)
                ax.errorbar(ew1com[det12req], ew2com[det12req], xerr=ewerr1com[det12req], yerr=ewerr2com[det12req],
                            capsize=0, ecolor=color[i], fmt=mark[i], color=color[i], label=legend[i])
                ax.quiver(ew1com[det1req], ewerr2com[det1req], 0, -1, width=0.003, scale=40, color=color[i])
                ax.errorbar(ew1com[det1req], ewerr2com[det1req], xerr=ewerr1com[det1req], fmt=" ", capsize=0,
                            ecolor=color[i])
                ax.quiver(ewerr1com[det2req], ew2com[det2req], -1, 0, width=0.003, scale=40, color=color[i])
                ax.errorbar(ewerr1com[det2req], ew2com[det2req], yerr=ewerr2com[det2req], fmt=" ", capsize=0,
                            ecolor=color[i])
                ax.quiver(ewerr1com[upreq], ewerr2com[upreq], -0.707, -0.707, width=0.003, scale=40, color=color[i])

        req = numpy.array([j not in objidsall for j in objidcom])
        det12req = numpy.logical_and(req, detection12)
        det1req = numpy.logical_and(req, detection1)
        det2req = numpy.logical_and(req, detection2)
        upreq = numpy.logical_and(req, uplimit)
        ax.errorbar(ew1com[det12req], ew2com[det12req], xerr=ewerr1com[det12req], yerr=ewerr2com[det12req],
                    capsize=0, ecolor=colordef, fmt=fmtdef, color=colordef)
        ax.quiver(ew1com[det1req], ewerr2com[det1req], 0, -0.5, scale=10, color=colordef)
        ax.errorbar(ew1com[det1req], ewerr2com[det1req], xerr=ewerr1com[det1req], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[det2req], ew2com[det2req], -0.5, 0, scale=10, color=colordef)
        ax.errorbar(ewerr1com[det2req], ew2com[det2req], yerr=ewerr2com[det2req], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[upreq], ewerr2com[upreq], -0.354, -0.354, scale=10, color=colordef)
        # ax.legend()
    else:
        ax.errorbar(ew1com[detection12], ew2com[detection12], xerr=ewerr1com[detection12], yerr=ewerr2com[detection12],
                    capsize=0, ecolor=colordef, fmt=fmtdef, color=colordef)
        ax.quiver(ew1com[detection1], ewerr2com[detection1], 0, -0.5, scale=10, color=colordef)
        ax.errorbar(ew1com[detection1], ewerr2com[detection1], xerr=ewerr1com[detection1], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[detection2], ew2com[detection2], -0.5, 0, scale=10, color=colordef)
        ax.errorbar(ewerr1com[detection2], ew2com[detection2], yerr=ewerr2com[detection2], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[uplimit], ewerr2com[uplimit], -0.354, -0.354, scale=10, color=colordef)

    z = numpy.polyfit(ew1com[detection12], ew2com[detection12], 1)
    p = numpy.poly1d(z)
    ax.plot(ew1com[detection12], p(ew1com[detection12]), "k")
    # ymin95, ymax95, x95 = confband(ew1com[detection12], ew2com[detection12], z[0], z[1])
    # ax.fill_between(x95, ymin95, ymax95, facecolor="gray", alpha=0.3)

    ax.set_xlabel(r"$\lambda$ %d EW (m$\AA$)" % DIBwav1)
    ax.set_ylabel(r"$\lambda$ %d EW (m$\AA$)" % DIBwav2)

    r, p = scipy.stats.pearsonr(ew1com[detection12], ew2com[detection12])
    # n = numpy.sum(detection12)
    # z = 0.5*numpy.log((1+r)/(1-r))
    # za = scipy.stats.norm.ppf(0.5 + 0.5 * alpha)
    # zl = z - za * math.sqrt(1/(n-3))
    # zu = z + za * math.sqrt(1/(n-3))
    # rhol = (math.exp(2 * zl) - 1) / (math.exp(2 * zl) + 1)
    # rhou = (math.exp(2 * zu) - 1) / (math.exp(2 * zu) + 1)
    rhol, rhou = 0, 0

    return r, p, numpy.sum(detection12), rhol, rhou


def DIBcorrelation2(ax, ew1, ewerr1, objid1, objname1, DIBwav1, ew2, ewerr2, objid2, objname2, DIBwav2, xlabeltext, ylabeltext,
                    stylefile="INDEF", colordef="k", fmtdef="o", msdef=3., alpha=0.95):
    objidcom = numpy.intersect1d(objid1, objid2)
    ew1com = numpy.array([ew1[i] for i in objidcom])
    ewerr1com = numpy.array([ewerr1[i] for i in objidcom])
    objname1com = numpy.array(["{}({})".format(objname1[i], i) for i in objidcom])
    ew2com = numpy.array([ew2[i] for i in objidcom])
    ewerr2com = numpy.array([ewerr2[i] for i in objidcom])
    objname2com = numpy.array(["{}({})".format(objname2[i], i) for i in objidcom])

    detection12 = (ew1com != 0.) & (ew2com != 0.)
    detection1 = (ew1com != 0.) & (ew2com == 0.)
    detection2 = (ew1com == 0.) & (ew2com != 0.)
    uplimit = (ew1com == 0.) & (ew2com == 0.)

    if stylefile != "INDEF":
        objids, ps, mark, color, legend = styleTxtParsar(stylefile)
        objidsall = []
        for i in objids:
            objidsall += i
        objidsallset = list(set(objidsall))
        if len(objidsallset) != len(objidsall):
            print("Error: overlapped object IDs.")
            sys.exit()
        for i in range(len(objids)):
            req = numpy.array([j in objids[i] for j in objidcom])
            if any(req):
                det12req = numpy.logical_and(req, detection12)
                det1req = numpy.logical_and(req, detection1)
                det2req = numpy.logical_and(req, detection2)
                upreq = numpy.logical_and(req, uplimit)
                ax.errorbar(ew1com[det12req], ew2com[det12req], xerr=ewerr1com[det12req], yerr=ewerr2com[det12req],
                            capsize=0, ecolor=color[i], fmt=mark[i], color=color[i], label=legend[i])
                ax.quiver(ew1com[det1req], ewerr2com[det1req], 0, -1, width=0.003, scale=40, color=color[i])
                ax.errorbar(ew1com[det1req], ewerr2com[det1req], xerr=ewerr1com[det1req], fmt=" ", capsize=0,
                            ecolor=color[i])
                ax.quiver(ewerr1com[det2req], ew2com[det2req], -1, 0, width=0.003, scale=40, color=color[i])
                ax.errorbar(ewerr1com[det2req], ew2com[det2req], yerr=ewerr2com[det2req], fmt=" ", capsize=0,
                            ecolor=color[i])
                ax.quiver(ewerr1com[upreq], ewerr2com[upreq], -0.707, -0.707, width=0.003, scale=40, color=color[i])
                textx = numpy.array([])
                texty = numpy.array([])
                textobj = numpy.array([])
                textx = numpy.append(textx, ew1com[det12req])
                texty = numpy.append(texty, ew2com[det12req])
                textobj = numpy.append(textobj, objname1com[det12req])
                textx = numpy.append(textx, ew1com[det1req])
                texty = numpy.append(texty, ewerr2com[det1req])
                textobj = numpy.append(textobj, objname1com[det1req])
                textx = numpy.append(textx, ewerr1com[det2req])
                texty = numpy.append(texty, ew2com[det2req])
                textobj = numpy.append(textobj, objname1com[det2req])
                textx = numpy.append(textx, ewerr1com[upreq])
                texty = numpy.append(texty, ewerr2com[upreq])
                textobj = numpy.append(textobj, objname1com[upreq])
                for k in range(textx.size):
                    plt.text(textx[k], texty[k], textobj[k], fontsize=5)

        req = numpy.array([j not in objidsall for j in objidcom])
        det12req = numpy.logical_and(req, detection12)
        det1req = numpy.logical_and(req, detection1)
        det2req = numpy.logical_and(req, detection2)
        upreq = numpy.logical_and(req, uplimit)
        ax.errorbar(ew1com[det12req], ew2com[det12req], xerr=ewerr1com[det12req], yerr=ewerr2com[det12req],
                    capsize=0, ecolor=colordef, fmt=fmtdef, color=colordef)
        ax.quiver(ew1com[det1req], ewerr2com[det1req], 0, -0.5, scale=10, color=colordef)
        ax.errorbar(ew1com[det1req], ewerr2com[det1req], xerr=ewerr1com[det1req], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[det2req], ew2com[det2req], -0.5, 0, scale=10, color=colordef)
        ax.errorbar(ewerr1com[det2req], ew2com[det2req], yerr=ewerr2com[det2req], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[upreq], ewerr2com[upreq], -0.354, -0.354, scale=10, color=colordef)
        textx = numpy.array([])
        texty = numpy.array([])
        textobj = numpy.array([])
        textx = numpy.append(textx, ew1com[det12req])
        texty = numpy.append(texty, ew2com[det12req])
        textobj = numpy.append(textobj, objname1com[det12req])
        textx = numpy.append(textx, ew1com[det1req])
        texty = numpy.append(texty, ewerr2com[det1req])
        textobj = numpy.append(textobj, objname1com[det1req])
        textx = numpy.append(textx, ewerr1com[det2req])
        texty = numpy.append(texty, ew2com[det2req])
        textobj = numpy.append(textobj, objname1com[det2req])
        textx = numpy.append(textx, ewerr1com[upreq])
        texty = numpy.append(texty, ewerr2com[upreq])
        textobj = numpy.append(textobj, objname1com[upreq])
        for k in range(textx.size):
            plt.text(textx[k], texty[k], textobj[k], fontsize=5)

        # ax.legend()
    else:
        ax.errorbar(ew1com[detection12], ew2com[detection12], xerr=ewerr1com[detection12], yerr=ewerr2com[detection12],
                    capsize=0, ecolor=colordef, fmt=fmtdef, color=colordef)
        ax.quiver(ew1com[detection1], ewerr2com[detection1], 0, -0.5, scale=10, color=colordef)
        ax.errorbar(ew1com[detection1], ewerr2com[detection1], xerr=ewerr1com[detection1], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[detection2], ew2com[detection2], -0.5, 0, scale=10, color=colordef)
        ax.errorbar(ewerr1com[detection2], ew2com[detection2], yerr=ewerr2com[detection2], fmt=" ", capsize=0,
                    ecolor=colordef)
        ax.quiver(ewerr1com[uplimit], ewerr2com[uplimit], -0.354, -0.354, scale=10, color=colordef)
        textx = numpy.array([])
        texty = numpy.array([])
        textobj = numpy.array([])
        textx = numpy.append(textx, ew1com[detection12])
        texty = numpy.append(texty, ew2com[detection12])
        textobj = numpy.append(textobj, objname1com[detection12])
        textx = numpy.append(textx, ew1com[detection1])
        texty = numpy.append(texty, ewerr2com[detection1])
        textobj = numpy.append(textobj, objname1com[detection1])
        textx = numpy.append(textx, ewerr1com[detection2])
        texty = numpy.append(texty, ew2com[detection2])
        textobj = numpy.append(textobj, objname1com[detection2])
        textx = numpy.append(textx, ewerr1com[uplimit])
        texty = numpy.append(texty, ewerr2com[uplimit])
        textobj = numpy.append(textobj, objname1com[uplimit])
        for k in range(textx.size):
            plt.text(textx[k], texty[k], textobj[k], fontsize=5)

    z = numpy.polyfit(ew1com[detection12], ew2com[detection12], 1)
    p = numpy.poly1d(z)
    ewfine = numpy.arange(numpy.amin(ew1com[detection12]), numpy.amax(ew1com[detection12]), 0.1)
    ax.plot(ewfine, p(ewfine), "k")
    # ymin95, ymax95, x95 = confband(ew1com[detection12], ew2com[detection12], z[0], z[1])
    # ax.fill_between(x95, ymin95, ymax95, facecolor="gray", alpha=0.3)

    ax.set_xlabel(xlabeltext)
    ax.set_ylabel(ylabeltext)

    r, p = scipy.stats.pearsonr(ew1com[detection12], ew2com[detection12])
    # n = numpy.sum(detection12)
    # z = 0.5*numpy.log((1+r)/(1-r))
    # za = scipy.stats.norm.ppf(0.5 + 0.5 * alpha)
    # zl = z - za * math.sqrt(1/(n-3))
    # zu = z + za * math.sqrt(1/(n-3))
    # rhol = (math.exp(2 * zl) - 1) / (math.exp(2 * zl) + 1)
    # rhou = (math.exp(2 * zu) - 1) / (math.exp(2 * zu) + 1)
    rhol, rhou = 0, 0

    return r, p, numpy.sum(detection12), rhol, rhou


def Vcorrelation(ax, ew1, objid1, DIBwav1, ew2, objid2, DIBwav2, stylefile="INDEF", colordef="k",
                 fmtdef="o", msdef=3.):
    objidcom = numpy.intersect1d(objid1, objid2)
    ew1com = numpy.array([ew1[i] for i in objidcom])
    ew2com = numpy.array([ew2[i] for i in objidcom])

    if stylefile != "INDEF":
        objids, ps, mark, color, legend = styleTxtParsar(stylefile)
        objidsall = []
        for i in objids:
            objidsall += i
        objidsallset = list(set(objidsall))
        if len(objidsallset) != len(objidsall):
            print("Error: overlapped object IDs.")
            sys.exit()
        for i in range(len(objids)):
            req = numpy.array([j in objids[i] for j in objidcom])
            if any(req):
                ax.scatter(ew1com[req], ew2com[req], marker=mark[i], c=color[i], label=legend[i])

        req = numpy.array([j not in objidsall for j in objidcom])
        ax.scatter(ew1com[req], ew2com[req], marker=fmtdef, c=colordef)  # , label='FWHM$_{10780}>1.5\AA$')
        # ax.legend()
    else:
        ax.scatter(ew1com, ew2com, marker=fmtdef, c=colordef)

    ax.plot([min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com))],
            [min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com))], "k")
    ax.plot([min(min(ew1com), min(ew2com)) + 10, max(max(ew1com), max(ew2com))],
            [min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com)) - 10], "k--")
    ax.plot([min(min(ew1com), min(ew2com)) + 20, max(max(ew1com), max(ew2com))],
            [min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com)) - 20], "k--")
    ax.plot([min(min(ew1com), min(ew2com)) + 30, max(max(ew1com), max(ew2com))],
            [min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com)) - 30], "k--")
    ax.plot([min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com)) - 30],
            [min(min(ew1com), min(ew2com)) + 30, max(max(ew1com), max(ew2com))], "k--")
    ax.plot([min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com)) - 20],
            [min(min(ew1com), min(ew2com)) + 20, max(max(ew1com), max(ew2com))], "k--")
    ax.plot([min(min(ew1com), min(ew2com)), max(max(ew1com), max(ew2com)) - 10],
            [min(min(ew1com), min(ew2com)) + 10, max(max(ew1com), max(ew2com))], "k--")

    ax.set_xlabel(r"$v$ ($\lambda$%d) (km/s)" % DIBwav1)
    ax.set_ylabel(r"$v$ ($\lambda$%d) (km/s)" % DIBwav2)

    return ew2com - ew1com


if __name__ == '__main__':
    conn, cur = openproject()

    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=str, help="Output pdf")
    # parser.add_argument("dibid1", type=int, help="DIB ID")
    # parser.add_argument("dibid2", type=int, help="DIB ID")
    parser.add_argument("-s", "--stylefile", type=str, default="INDEF", help="Style file")
    parser.add_argument("-n", "--norm", action="store_true")
    parser.add_argument("-l", "--log", action="store_true")

    args = parser.parse_args()
    dibid1 = [33, 53, 39, 54, 55, 58, 40, 41, 42, 43, 60, 61, 44, 64, 66, 35, 36, 67, 68, 69, 70, 72, 73,
                74, 28, 75, 76, 37, 679, 30, 80, 81, 82, 45, 680, 46, 47, 48, 85, 49, 86, 87, 50, 80, 51, 90,
                52, 91, 38]  # [33, 34, 35, 37, 38]  # args.dibid1
    dibid2 = [33, 53, 39, 54, 55, 58, 40, 41, 42, 43, 60, 61, 44, 64, 66, 35, 36, 67, 68, 69, 70, 72, 73,
                74, 28, 75, 76, 37, 679, 30, 80, 81, 82, 45, 680, 46, 47, 48, 85, 49, 86, 87, 50, 80, 51, 90,
                52, 91, 38]  # args.dibid2
    dibidselect2 = [33, 53, 39, 54, 55, 58, 40, 41, 42, 43, 60, 61, 44, 64, 66, 35, 36, 67, 68, 69, 70, 72, 73,
                74, 28, 75, 76, 37, 679, 30, 80, 81, 82, 45, 680, 46, 47, 48, 85, 49, 86, 87, 50, 80, 51, 90,
                52, 91, 38]
    dibidselect = [33, 53, 39, 58, 40, 41, 42, 61, 44, 64, 66, 35, 36, 68, 69, 70, 73, 74, 75, 76, 37, 679, 80, 81, 82,
                   45, 680, 46, 47, 48, 85, 49, 86, 87, 50, 80, 51, 90, 52, 91, 38]
    dibidall = [33, 34, 53, 39, 54, 55, 56, 58, 40, 41, 42, 43, 60, 61, 44, 64, 66, 35, 36, 67, 68, 69, 70, 71, 72, 73,
                74, 28, 75, 76, 37, 679, 30, 79, 80, 81, 82, 45, 680, 46, 47, 48, 85, 49, 86, 87, 50, 80, 51, 89, 90,
                52, 91, 38]
    output = args.output
    stylefile = args.stylefile
    norm = args.norm
    log = args.log
    objectIDinput = [50, 75, 23, 22, 24, 33, 43, 30, 123, 31, 42, 124, 90, 246, 32, 27, 63, 45, 59, 41, 66, 40, 39,
                     14, 10, 12, 13, 15, 9, 132, 144]
    print(len(objectIDinput))
    objstr = ''
    for i in objectIDinput:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    pp = PdfPages(output)
    pp2 = PdfPages("Correlation/pair_select2_hist.pdf")
    wf = open("Correlation/corrnorm_select2_tabletext_r.csv", "w")
    wf2 = open("Correlation/corrnorm_select2_tabletext_n.csv", "w")

    plt.figure(figsize=(8, 8))
    headtext = "\\colhead{{DIB}} & \\colhead{{N}} "
    ctext = "cc"
    for n in range(len(dibid1)):
        pairnum = []
        cur.execute(
            "select z.objectID,x.EW,x.EWerr,z.E_BV,z.objectname from DIBmeasurement as x join combinesummary as y using(combineID) "
            "join object as z using(objectID) where x.DIBID=%d and x.primaryflag = 1 and z.objectID IN (%s);" % (
                dibid1[n], objstr))
        rows = cur.fetchall()
        objectid1 = numpy.array([int(i[0]) for i in rows])
        ew1, ewerr1, ebv1, objname1 = {}, {}, {}, {}
        numdet = 0
        for i in rows:
            if norm:
                ew1[i[0]] = float(i[1])/float(i[3])
                ewerr1[i[0]] = float(i[2])/float(i[3])
            else:
                ew1[i[0]] = float(i[1])#/float(i[3])
                ewerr1[i[0]] = float(i[2])#/float(i[3])

            objname1[i[0]] = i[4]
            if float(i[1]) > 0:
                numdet += 1
        DIBinfo = GetDIBLine(dibid1[n])
        DIBwav1 = DIBinfo[1]
        tabletext = "{:.1f} \t {} ".format(DIBwav1, numdet)
        tabletext2 = "{:.1f} \t {} ".format(DIBwav1, numdet)
        headtext += " & \\colhead{{{:.1f}}}".format(DIBwav1)
        ctext += "c"
        for m in range(n+1, len(dibid2)):
            # if dibid1[n] != dibid2[m]:
            cur.execute(
                "select z.objectID,x.EW,x.EWerr,z.E_BV,z.objectname from DIBmeasurement as x join combinesummary as y using(combineID) "
                "join object as z using(objectID) where x.DIBID=%d and x.primaryflag = 1 and z.objectID IN (%s);" % (
                    dibid2[m], objstr))
            rows = cur.fetchall()
            objectid2 = numpy.array([int(i[0]) for i in rows])
            ew2, ewerr2, objname2 = {}, {}, {}
            for i in rows:
                if norm:
                    ew2[i[0]] = float(i[1])/float(i[3])
                    ewerr2[i[0]] = float(i[2])/float(i[3])
                else:
                    ew2[i[0]] = float(i[1])
                    ewerr2[i[0]] = float(i[2])
                objname2[i[0]] = i[4]
            DIBinfo = GetDIBLine(dibid2[m])
            DIBwav2 = DIBinfo[1]

            fig, ax = plt.subplots(1, 1)

            if norm:
                xl = r"$\lambda$ %d EW / E(B-V) (m$\AA$/mag)" % DIBwav1
                yl = r"$\lambda$ %d EW / E(B-V) (m$\AA$/mag)" % DIBwav2
            else:
                xl = r"$\lambda$ %d EW (m$\AA$)" % DIBwav1
                yl = r"$\lambda$ %d EW (m$\AA$)" % DIBwav2

            r, p, num, rhol, rhou = DIBcorrelation2(ax, ew1, ewerr1, objectid1, objname1, DIBwav1, ew2, ewerr2,
                                                    objectid2, objname2,
                                                    DIBwav2, xl, yl,  stylefile=stylefile)
            if log:
                plt.xscale("log")
                plt.yscale("log")
            plt.title("{} - {} (R={:.3f} (N={}))".format(DIBwav1, DIBwav2, r, num))
            plt.savefig(pp, format="pdf")
            plt.clf()
            pairnum.append(num)
            tabletext += "\t {:.2f} ".format(r)
            tabletext2 += "\t {} ".format(num)

            # if r > 0.9:
            #     # tabletext += "& \\textcolor{{red}}{{{:.2f}({})}} ".format(r, num)
            #     tabletext += "& \\textcolor{{red}}{{{:.2f}}} ".format(r)
            #     tabletext2 += "& \\textcolor{{red}}{{{}}} ".format(num)
            # elif r < 0.6:
            #     # tabletext += "& \\textcolor{{blue}}{{{:.2f}({})}} ".format(r, num)
            #     tabletext += "& \\textcolor{{blue}}{{{:.2f}}} ".format(r)
            #     tabletext2 += "& \\textcolor{{blue}}{{{}}} ".format(num)
            # else:
            #     # tabletext += "& \\textcolor{{black}}{{{:.2f}({})}} ".format(r, num)
            #     tabletext += "& \\textcolor{{black}}{{{:.2f}}} ".format(r)
            #     tabletext2 += "& \\textcolor{{black}}{{{}}} ".format(r)

        plt.hist(pairnum)
        plt.title("{}".format(DIBwav1))
        plt.savefig(pp2, format="pdf")
        plt.clf()
        tabletext += "\n"
        tabletext2 += "\n"
        wf.write(tabletext)
        wf2.write(tabletext2)

    wf.write(headtext)
    wf.write("\n")
    wf.write(ctext)
    wf.write("\n")

    conn.close()
    pp.close()
    pp2.close()
    wf.close()
    wf2.close()
