#!/usr/bin/env python
# -*- coding:utf-8 -*-
import sys,os,datetime,glob
import matplotlib.pyplot as plt
import numpy,math
import scipy.optimize
import scipy.stats

from open_mysql_project import openproject

#python DIBcorrelations.py paramlist0 paramlist1 output2 ebvtable3 colortable4 -pair(*) -ebv(*) -conf(*) -reg(*) confidence=0.95 xmin(*)= xmax(*)= ymin(*)= ymax(*)= xstep(*)= ystep(*)=

#   Start developping: 2016-03-29 9:43(JST) in the flight from SFO to KIX.
#   Description: plot the correlations of DIBs
#   Necessary inputs:
#       paramlist0: The ascii list of the directories containing the DIB parameter files for each observed stars.
#       paramlist1: The ascii list of the directories containing the DIB parameter files for each observed stars.
#       output2: The letters attached to the top of the output files.
#       ebvtable3: The ascii file containing table of E(B-V) of observed stars
#       colortable4: The ascii file containing table of color, point types, sizes of observed stars
#
#   Optional inputs:
#       confidence: indicate the value of confidnece lebel
#       -pair:  cancel to create the correlation plot of DIB pairs
#       -ebv:   cancel to create the correlation between DIBs and E(B-V)
#       -conf:  cancel to draw the confidence bands in the plots
#       -reg:   cancel to draw the regression lines in the plots
#       xmin:   indicate the value of the xmin value in all plots
#       ymin:   indicate the value of the ymin value in all plots
#       xmax:   indicate the value of the xmax value in all plots
#       ymax:   indicate the value of the ymax value in all plots
#       xstep:  indicate the value of the xstep value in all plots
#       ystep:  indicate the value of the ystep value in all plots

from matplotlib import cm

from matplotlib.backends.backend_pdf import PdfPages

#definition of functions

def confband(xd,yd,a,b,x,conf=0.95):
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
    alpha=1.-conf	# significance
    n=xd.size	# data sample size
    
    if x==None: x=numpy.linspace(xd.min(),xd.max(),100); print("none")
    
    # Predicted values (best-fit model)
    y=a*x+b
    
    # Auxiliary definitions
    #sd=scatterfit(xd,yd,a,b)	# Scatter of data about the model
    sd = 1./(n-2.) * numpy.sum((yd-a*xd-b)**2); sd = numpy.sqrt(sd)
    sxd=numpy.sum((xd-xd.mean())**2)
    sx=(x-xd.mean())**2	# array
    print(xd.mean())
    
    # Quantile of Student's t distribution for p=1-alpha/2
    q=scipy.stats.t.ppf(1.-alpha/2.,n-2)
    
    # Confidence band
    #print q, sd, n, sx, sxd
    dy=q*sd*numpy.sqrt( 1./n + sx/sxd )
    #print x
    ucb=y+dy	# Upper confidence band
    lcb=y-dy	# Lower confidence band
    
    return lcb,ucb,x

def regressionline(x,y,pstatus,pregression):
    
    ccx = []
    ccy = []
    for i in range(len(x)):
        if pstatus[i] == ["Detection" for n in range(len(pstatus[i]))] and pregression[i] == "in":
            ccx.append(x[i])
            ccy.append(y[i])

    ccxarray = numpy.array(ccx)
    ccyarray = numpy.array(ccy)

    sumx = numpy.sum(ccxarray)
    sumy = numpy.sum(ccyarray)
    sumxy = numpy.sum(ccxarray*ccyarray)
    sumxx = numpy.sum(ccxarray**2)
    sumyy = numpy.sum(ccyarray**2)
    narray = len(ccxarray)
        
    ssxx = sumxx - narray * (sumx / narray)**2
    ssyy = sumyy - narray * (sumy / narray)**2
    ssxy = sumxy - narray * sumx / narray * sumy / narray
    paa = ssxy / ssxx
    pab = sumy / narray - paa * sumx / narray

    return paa, pab, ccxarray, ccyarray

def linearfunc1(x,p):
    func = []
    for i in range(len(x)):
        if p[1] * x[i] + p[0] < 0:
            #print 0
            func.append(0.)
        else:
            func.append(p[1] * x[i] + p[0])
    afunc = numpy.array(func)
    return (afunc)

def linearfunc2(x,p):
    func = p[0] * x
    return (func)

def residue1(p,y,x,err):
    res = ((y-linearfunc1(x,p))/err)
    #res = (y-linearfunc1(x,p))
    return (res)

def residue2(p,y,x,err):
    res = ((y-linearfunc2(x,p))/err)
    return (res)

#read parameters of DIBs

filename = sys.argv[1:]

flagpair = 1
flagebv = 1
flagconf = 1
flagreg = 1
flagxmin = 0
flagxmax = 0
flagymin = 0
flagymax = 0
flagxstep = 0
flagystep = 0
flagconfidence = 0
xmin = 0.
xmax = 0.
ymin = 0.
ymax = 0.
xstep = 0.
ystep = 0.
confidence = 0.95
ll = 1000
cfactor = 1.0e+3

for i in range(4,len(filename)):
    if filename[i].find("-pair") != -1:
        flagpair = 0
    if filename[i].find("-ebv") != -1:
        flagebv = 0
    if filename[i].find("-conf") != -1:
        flagconf = 0
    if filename[i].find("-reg") != -1:
        flagreg = 0
    if filename[i].find("xmin=") != -1:
        flagxmin = 1
        xmin = float(filename[i].split("=")[1])
    if filename[i].find("xmax=") != -1:
        flagxmax = 1
        xmax = float(filename[i].split("=")[1])
    if filename[i].find("ymin=") != -1:
        flagymin = 1
        ymin = float(filename[i].split("=")[1])
    if filename[i].find("ymax=") != -1:
        flagymax = 1
        ymax = float(filename[i].split("=")[1])
    if filename[i].find("xstep=") != -1:
        flagxstep = 1
        xstep = float(filename[i].split("=")[1])
    if filename[i].find("ystep=") != -1:
        flagystep = 1
        ystep = float(filename[i].split("=")[1])
    if filename[i].find("confidence=") != -1:
        flagconfidence = 1
        confidence = float(filename[i].split("=")[1])

dibdirlist = open(filename[0],"r")
dibdirlines = dibdirlist.readlines()
dibdirlist.close()



dir1_DIBname = []
param1_DIBname = []
param1_SN = []
param1_inverseSN = []
param1_From = []
param1_To = []
param1_EW = []
param1_FWZI = []
param1_Depth = []
param1_Status = []
param1_EWerr = []
param1_Thres = []
param1_object = []

for i in range(len(dibdirlines)):
    
    dibpalines = glob.glob("%s/*param" % dibdirlines[i].split()[0].rstrip("/"))
    
    dir1_DIBname.append(dibdirlines[i].split("/")[-1].split()[0])
    
    param1_DIBname.append([])
    param1_SN.append([])
    param1_inverseSN.append([])
    param1_From.append([])
    param1_To.append([])
    param1_EW.append([])
    param1_FWZI.append([])
    param1_Depth.append([])
    param1_Status.append([])
    param1_EWerr.append([])
    param1_Thres.append([])
    param1_object.append([])
    
    for j in range(len(dibpalines)):
        
        paramfile = open(dibpalines[j],"r")
        paramlines = paramfile.readlines()
        paramfile.close()
                
        param1_DIBname[i].append(float(paramlines[0].split()[1]))
        param1_SN[i].append(float(paramlines[1].split()[1]))
        param1_inverseSN[i].append(float(paramlines[2].split()[1]))
        param1_From[i].append(float(paramlines[3].split()[1]))
        param1_To[i].append(float(paramlines[4].split()[1]))
        param1_EW[i].append(float(paramlines[5].split()[1])*cfactor)
        param1_FWZI[i].append(float(paramlines[6].split()[1]))
        param1_Depth[i].append(float(paramlines[7].split()[1]))
        param1_Status[i].append(paramlines[8].split()[1])
        param1_EWerr[i].append(float(paramlines[9].split()[1])*cfactor)
        param1_Thres[i].append(float(paramlines[10].split()[1]))
        param1_object[i].append(paramlines[11].split()[1])

dibdirlist = open(filename[1],"r")
dibdirlines = dibdirlist.readlines()
dibdirlist.close()



dir2_DIBname = []
param2_DIBname = []
param2_SN = []
param2_inverseSN = []
param2_From = []
param2_To = []
param2_EW = []
param2_FWZI = []
param2_Depth = []
param2_Status = []
param2_EWerr = []
param2_Thres = []
param2_object = []

for i in range(len(dibdirlines)):
    
    dibpalines = glob.glob("%s/*param" % dibdirlines[i].split()[0].rstrip("/"))
    
    dir2_DIBname.append(dibdirlines[i].split("/")[-1].split()[0])
    
    param2_DIBname.append([])
    param2_SN.append([])
    param2_inverseSN.append([])
    param2_From.append([])
    param2_To.append([])
    param2_EW.append([])
    param2_FWZI.append([])
    param2_Depth.append([])
    param2_Status.append([])
    param2_EWerr.append([])
    param2_Thres.append([])
    param2_object.append([])
    
    for j in range(len(dibpalines)):
        
        paramfile = open(dibpalines[j],"r")
        paramlines = paramfile.readlines()
        paramfile.close()
        
        param2_DIBname[i].append(float(paramlines[0].split()[1]))
        param2_SN[i].append(float(paramlines[1].split()[1]))
        param2_inverseSN[i].append(float(paramlines[2].split()[1]))
        param2_From[i].append(float(paramlines[3].split()[1]))
        param2_To[i].append(float(paramlines[4].split()[1]))
        param2_EW[i].append(float(paramlines[5].split()[1])*cfactor)
        param2_FWZI[i].append(float(paramlines[6].split()[1]))
        param2_Depth[i].append(float(paramlines[7].split()[1]))
        param2_Status[i].append(paramlines[8].split()[1])
        param2_EWerr[i].append(float(paramlines[9].split()[1])*cfactor)
        param2_Thres[i].append(float(paramlines[10].split()[1]))
        param2_object[i].append(paramlines[11].split()[1])

outputroot = filename[2]

ebvtable = open(filename[3],"r")
ebvlines = ebvtable.readlines()
ebvtable.close()

ebv_object = []
ebv_ebv = []
ebv_ebverror = []

for i in range(len(ebvlines)):
    ebv_object.append(ebvlines[i].split()[0])
    ebv_ebv.append(float(ebvlines[i].split()[1]))
    ebv_ebverror.append(float(ebvlines[i].split()[2]))

styletable = open(filename[4],"r")
stylelines = styletable.readlines()
styletable.close()

style_object = []
style_type = []
style_size = []
style_color = []
style_regression = []

for i in range(len(stylelines)):
    style_object.append(stylelines[i].split()[0])
    style_type.append(stylelines[i].split()[1])
    style_size.append(stylelines[i].split()[2])
    style_color.append(stylelines[i].split()[3])
    style_regression.append(stylelines[i].split()[4])

fig, ax = plt.subplots()
plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
[i.set_linewidth(0.7) for i in ax.spines.itervalues()]

if flagpair == 1:
    pp = PdfPages("%s_DIBpairs.pdf" % outputroot)
    for i in range(len(dir1_DIBname)):
        for j in range(len(dir2_DIBname)):
            if dir1_DIBname[i] != dir2_DIBname[j]:
                print "Now working on %s - %s correlation plot..." % (dir1_DIBname[i],dir2_DIBname[j])
                x = []
                y = []
                xerr = []
                yerr = []
                
                pobj = []
                pstatus = []
                ptype = []
                psize = []
                pcolor = []
                pregression = []
                for ii in range(len(param1_object[i])):
                    for jj in range(len(param2_object[j])):
                        if param1_object[i][ii] == param2_object[j][jj]:
                            y.append(param1_EW[i][ii])
                            x.append(param2_EW[j][jj])
                            yerr.append(param1_EWerr[i][ii])
                            xerr.append(param2_EWerr[j][jj])
                            pobj.append(param1_object[i][ii])
                            pstatus.append([param2_Status[j][jj],param1_Status[i][ii]])
                for n in range(len(pobj)):
                    for k in range(len(style_object)):
                        if pobj[n] == style_object[k]:
                            ptype.append(style_type[k])
                            psize.append(style_size[k])
                            pcolor.append(style_color[k])
                            pregression.append(style_regression[k])

                if flagreg == 1:
                    paa, pab, ccxarray, ccyarray = regressionline(x,y,pstatus,pregression)
                    plt.plot([0.0,max(x)*2], [pab, paa*max(x)*2 + pab], "k", label="w/ intercept")

                if flagconf == 1:
                    if flagreg != 1:
                        paa, pab, ccxarray, ccyarray = regressionline(x,y,pstatus,pregression)
                    x95 = numpy.array([float(k)/float(ll) * max(x)*1.1 for k in range(ll+1)])
                    ymin95, ymax95, x95 = confband(ccxarray, ccyarray, paa, pab, x95)
                    plt.fill_between(x95,ymin95,ymax95,facecolor="black",alpha=0.3)
    

                for n in range(len(pobj)):
                    if pstatus[n] == ["Detection","Detection"]:
                        plt.errorbar(x[n], y[n], xerr=xerr[n], yerr=yerr[n],fmt=ptype[n], mfc=pcolor[n], markersize=psize[n], mec=pcolor[n], capsize=0, ecolor=pcolor[n])
                    elif pstatus[n] == ["Detection","UpperLimit"]:
                        plt.quiver(x[n], yerr[n], 0, -1, width=0.003, scale=40, color=pcolor[n])
                        plt.errorbar(x[n], yerr[n], xerr=xerr[n], fmt=" ", capsize=0, ecolor=pcolor[n])
                    elif pstatus[n] == ["UpperLimit","Detection"]:
                        plt.quiver(xerr[n], y[n], -1, 0, width=0.003, scale=40, color=pcolor[n])
                        plt.errorbar(xerr[n], y[n], yerr=yerr[n], fmt=" ", capsize=0, ecolor=pcolor[n])
                    elif pstatus[n] == ["UpperLimit","UpperLimit"]:
                        plt.quiver(xerr[n], yerr[n], -0.707, -0.707, width=0.003, scale=40, color=pcolor[n])

                if flagxmin != 1:
                    xmin = 0.
                if flagxmax != 1:
                    xmax = max(x)*1.1
                if flagymin != 1:
                    ymin = 0.
                if flagymax != 1:
                    ymax = max(y)*1.1
                
                plt.xlim(xmin=xmin, xmax=xmax)
                plt.ylim(ymin=ymin, ymax=ymax)

                if flagxstep == 1:
                    plt.xticks(xmin, xmax, xstep)

                if flagystep == 1:
                    plt.yticks(ymin, ymax, ystep)

                plt.xlabel("%s" % dir2_DIBname[j], fontsize=15)
                plt.ylabel("%s" % dir1_DIBname[i], fontsize=15)

                plt.savefig("%s_%s_%s.eps" % (outputroot, dir2_DIBname[j], dir1_DIBname[i]))
                plt.savefig(pp, format="pdf", dpi=400)
                plt.clf()

    pp.close()


fig, ax = plt.subplots()
plt.figure(figsize=(9,5))
ax = fig.add_subplot(111)
[i.set_linewidth(0.7) for i in ax.spines.itervalues()]

if flagebv == 1:
    pp = PdfPages("%s_DIB_ebv.pdf" % outputroot)
    for i in range(len(dir1_DIBname)):
        print "Now working on E(B-V) - %s correlation plot..." % (dir1_DIBname[i])

        x = []
        y = []
        xerr = []
        yerr = []
        
        pobj = []
        pstatus = []
        ptype = []
        psize = []
        pcolor = []
        pregression = []
        for ii in range(len(param1_object[i])):
            for jj in range(len(ebv_object)):
                if param1_object[i][ii] == ebv_object[jj]:
                    x.append(ebv_ebv[jj])
                    y.append(param1_EW[i][ii])
                    xerr.append(ebv_ebverror[jj])
                    yerr.append(param1_EWerr[i][ii])
                    pobj.append(param1_object[i][ii])
                    pstatus.append([param1_Status[i][ii]])
        for n in range(len(pobj)):
            for k in range(len(style_object)):
                if pobj[n] == style_object[k]:
                    ptype.append(style_type[k])
                    psize.append(style_size[k])
                    pcolor.append(style_color[k])
                    pregression.append(style_regression[k])
            
        if flagreg == 1:
            paa, pab, ccxarray, ccyarray = regressionline(x,y,pstatus,pregression)
            plt.plot([0.0,max(x)*2], [pab, paa*max(x)*2 + pab], "k", label="w/ intercept")

        if flagconf == 1:
            if flagreg != 1:
                paa, pab, ccxarray, ccyarray = regressionline(x,y,pstatus,pregression)
            x95 = numpy.array([float(k)/float(ll) * max(x)*1.1 for k in range(ll+1)])
            ymin95, ymax95, x95 = confband(ccxarray, ccyarray, paa, pab, x95)
            plt.fill_between(x95,ymin95,ymax95,facecolor="black",alpha=0.3)
                
                
        for n in range(len(pobj)):
            if pstatus[n] == ["Detection"]:
                plt.errorbar(x[n], y[n], xerr=xerr[n], yerr=yerr[n],fmt=ptype[n], mfc=pcolor[n], markersize=psize[n], mec=pcolor[n], capsize=0, ecolor=pcolor[n])
            elif pstatus[n] == ["UpperLimit"]:
                plt.quiver(x[n], yerr[n], 0, -1, width=0.003, scale=40, color=pcolor[n])
                plt.errorbar(x[n], yerr[n], xerr=xerr[n], fmt=" ", capsize=0, ecolor=pcolor[n])
            
        if flagxmin != 1:
            xmin = 0.
        if flagxmax != 1:
            xmax = max(x)*1.1
        if flagymin != 1:
            ymin = 0.
        if flagymax != 1:
            ymax = max(y)*1.1
                
        plt.xlim(xmin=xmin, xmax=xmax)
        plt.ylim(ymin=ymin, ymax=ymax)
                
        if flagxstep == 1:
            plt.xticks(xmin, xmax, xstep)

        if flagystep == 1:
            plt.yticks(ymin, ymax, ystep)
        
        #plt.xlabel("N(C2) (10^13 cm^-2)")
        plt.xlabel("E(B-V)")
        plt.ylabel("%s" % dir1_DIBname[i])
                
        #plt.savefig("%s_c2_%s.eps" % (outputroot, dir1_DIBname[i]))
        plt.savefig("%s_ebv_%s.eps" % (outputroot, dir1_DIBname[i]))
        plt.savefig(pp, format="pdf", dpi=400)
        plt.clf()

    pp.close()

