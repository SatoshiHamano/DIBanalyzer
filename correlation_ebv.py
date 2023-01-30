#!/usr/bin/env python
# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import argparse
import numpy
import scipy.stats
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from add_lineDIB_mysql import GetDIBLine
import math


def confband(xd, yd, a, b, x=None, conf=0.95):
    """
        Calculates the confidence band of the linear regression model at the desired confidence
        level, using analytical methods. The 2sigma confidence interval is 95% sure to contain
        the best-fit regression line. This is not the same as saying it will contain 95% of
        the data points.
        Arguments:
        - conf: desired confidence level, by default 0.95 (2 sigma)
        - xd,yd: data arrays
        - a,b: linear fit parameters as in y=ax+b
        - x: (optional) array with x values to calculate the confidence band. If none is provided, will
        by default generate 100 points in the original x-range of the data.

        Returns:
        Sequence (lcb,ucb,x) with the arrays holding the lower and upper confidence bands
        corresponding to the [input] x array.
        Usage:
        >>> lcb,ucb,x=nemmen.confband(all.kp,all.lg,a,b,conf=0.95)
        calculates the confidence bands for the given input arrays
        >>> pylab.fill_between(x, lcb, ucb, alpha=0.3, facecolor='gray')
        plots a shaded area containing the confidence band
        References:
        1. http://en.wikipedia.org/wiki/Simple_linear_regression, see Section Confidence intervals
        2. http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
        Author: Rodrigo Nemmen
        v1 Dec. 2011
        v2 Jun. 2012: corrected bug in computing dy
        """
    alpha = 1. - conf  # significance
    n = xd.size  # data sample size

    if x == None: x = numpy.linspace(xd.min(), xd.max(), 100)#; print("none")

    # Predicted values (best-fit model)
    y = a * x + b

    # Auxiliary definitions
    # sd=scatterfit(xd,yd,a,b)	# Scatter of data about the model
    sd = 1. / (n - 2.) * numpy.sum((yd - a * xd - b) ** 2);
    sd = numpy.sqrt(sd)
    sxd = numpy.sum((xd - xd.mean()) ** 2)
    sx = (x - xd.mean()) ** 2  # array
    #print(xd.mean())

    # Quantile of Student's t distribution for p=1-alpha/2
    q = scipy.stats.t.ppf(1. - alpha / 2., n - 2)

    # Confidence band
    # print q, sd, n, sx, sxd
    dy = q * sd * numpy.sqrt(1. / n + sx / sxd)
    # print x
    ucb = y + dy  # Upper confidence band
    lcb = y - dy  # Lower confidence band

    return lcb, ucb, x



def ebv_correlation(ax, ebv, ew, ewerr, objectid, DIBwav, stylefile="INDEF", colordef="k", fmtdef="o", msdef=3., alpha=0.95):
    # pp = PdfPages(outputpdf)

    detection = ew != 0.
    uplimit = numpy.logical_not(detection)

    # plt.figure()
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
            req = numpy.array([j in objids[i] for j in objectid])
            if any(req):
                detreq = numpy.logical_and(req, detection)
                upreq = numpy.logical_and(req, uplimit)
                ax.errorbar(ebv[detreq], ew[detreq], yerr=ewerr[detreq], capsize=0, ecolor=color[i], fmt=mark[i],
                             color=color[i], label=legend[i])
                ax.quiver(ebv[upreq], ewerr[upreq], 0, -1, width=0.003, scale=40, color=color[i])

        req  = numpy.array([j not in objidsall for j in objectid])
        detreq = numpy.logical_and(req, detection)
        upreq = numpy.logical_and(req, uplimit)
        ax.errorbar(ebv[detreq], ew[detreq], yerr=ewerr[detreq], capsize=0, ecolor=colordef, fmt=fmtdef,
                     color=colordef)
        ax.quiver(ebv[upreq], ewerr[upreq], 0, -0.5, scale=10, color=colordef)
        ax.legend()
    else:
        ax.errorbar(ebv[detection], ew[detection], yerr=ewerr[detection], capsize=0, ecolor=colordef, fmt=fmtdef,
                     color=colordef)
        ax.quiver(ebv[uplimit], ewerr[uplimit], 0, -0.5, scale=10, color=colordef)

    z = numpy.polyfit(ebv[detection], ew[detection], 1)
    p = numpy.poly1d(z)
    ax.plot(ebv[detection], p(ebv[detection]), "k")
    ymin95, ymax95, x95 = confband(ebv[detection], ew[detection], z[0], z[1])
    ax.fill_between(x95,ymin95,ymax95,facecolor="gray",alpha=0.3)


    ax.set_xlabel("E(B-V) (mag)")
    ax.set_ylabel(r"EW (m$\AA$)")
    # ax.set_title(r"$\lambda$ %d" % DIBwav)
    # ax.set_yscale("log")

    # pp.close()

    r, p = scipy.stats.pearsonr(ebv[detection], ew[detection])
    n = numpy.sum(detection)
    z = 0.5*numpy.log((1+r)/(1-r))
    za = scipy.stats.norm.ppf(0.5 + 0.5 * alpha)
    zl = z - za * math.sqrt(1/(n-3))
    zu = z + za * math.sqrt(1/(n-3))
    rhol = (math.exp(2 * zl) - 1) / (math.exp(2 * zl) + 1)
    rhou = (math.exp(2 * zu) - 1) / (math.exp(2 * zu) + 1)

    return r, p, rhol, rhou


def styleTxtParsar(styletxt):
    rf = open(styletxt)
    rl = rf.readlines()
    rf.close()

    objids = [[int(j) for j in i.split("|")[0].split(",")] for i in rl]
    ps = [float(i.split("|")[1].split()[0]) for i in rl]
    mark = [i.split("|")[2].split()[0] for i in rl]
    color = [i.split("|")[3].split()[0] for i in rl]
    legend = [i.split("|")[4].rstrip("\n") for i in rl]

    return objids, ps, mark, color, legend


if __name__ == '__main__':
    conn, cur = openproject()

    parser = argparse.ArgumentParser()
    # parser.add_argument("output", type=str, help="Output pdf")
    # parser.add_argument("-o", "--objectID", type=int, help="object ID", nargs='*')
    # parser.add_argument("-d", "--dibid", type=int, default=0, help="DIB ID")
    parser.add_argument("-s", "--stylefile", type=str, default="INDEF", help="Style file")

    args = parser.parse_args()
    # objectID = args.objectID
    # dibid = args.dibid
    stylefile = args.stylefile
    # output = args.output

    dibidall = [33, 34, 53, 39, 54, 55, 56, 58, 40, 41, 42, 43, 60, 61, 44, 64, 66, 35, 36, 67, 68, 69, 70, 71, 72, 73,
                74, 28, 75, 76, 37, 679, 30, 79, 80, 81, 82, 45, 680, 46, 47, 48, 85, 49, 86, 87, 50, 80, 51, 89, 90,
                52, 91, 38]
    objectIDinput = [138, 50, 75, 23, 22, 24, 33, 43, 30, 123, 31, 42, 124, 90, 246, 32, 27, 63, 45, 59, 41, 66, 40, 39,
                     14, 10, 12, 13, 15, 9, 132, 144]
    print(len(objectIDinput))
    objstr = ''
    for i in objectIDinput:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    for dibid in dibidall:
        cur.execute(
            "select z.objectID,x.EW,x.EWerr,z.E_BV from DIBmeasurement as x join combinesummary as y using(combineID) "
            "join object as z using(objectID) where x.DIBID=%d and x.primaryflag = 1 and z.E_BV IS NOT NULL and objectID IN (%s);" % (dibid,objstr))
        rows = cur.fetchall()
        objectid = numpy.array([int(i[0]) for i in rows])
        ew = numpy.array([float(i[1]) for i in rows])
        ewerr = numpy.array([float(i[2]) for i in rows])
        ebv = numpy.array([float(i[3]) for i in rows])

        DIBinfo = GetDIBLine(dibid)
        DIBwav = DIBinfo[1]

        pp = PdfPages("Correlation/DIB{:.0f}_ebv.pdf".format(DIBwav))

        plt.figure(figsize=(3,2))
        fig, ax = plt.subplots(1,1)
        ebv_correlation(ax, ebv, ew, ewerr, objectid, DIBwav, stylefile=stylefile)
        plt.xlim(-0.2, 5.0)
        plt.title(DIBwav)
        plt.xlabel("E(B-V)")
        plt.ylabel("EW (mA)")
        plt.savefig(pp, format="pdf")
        plt.clf()

        pp.close()
        print("\\includegraphics[width=8cm,clip]{{/Users/hamano/PycharmProjects/DIBproject/Correlation/DIB{:.0f}_ebv.pdf}}".format(DIBwav))

    conn.close()
