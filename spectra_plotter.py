# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import math
import numpy as np
import copy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from Spec1Dtools import openspecfits, model_opener
import scipy.constants
from scipy import signal, interpolate
import glob
from add_lineDIB_mysql import *
from C2plotter import C2marks
from C2_parameters import *
from vac2air_spec import air2vac


def spectrum_plot(ax, spx, spy, xrange, yrange, spx_func="NA", spy_func="NA", colors=["k", "b", "r", "g"], linew=[1.],
                  lines=["-"], labels=[""], colors_func=["k", "b", "r", "g"], linew_func=[1.],
                  lines_func=["-"], labels_func=[""], yshift=[0], yfactor=[1.], xlshift=[0], xvshift=[0],
                  xticks_val="NA", yticks_val="NA", xticks_label="NA", yticks_label="NA", xaxis_label="NA",
                  yaxis_label="NA", grid_flag=False, legend_flag=False, legend_loc=0):
    # spx, spy: list of numpy.array
    # xrange, yrange: [xmin, xmax], [ymin, ymax]

    # colors, linew, lines, labels: list
    # yshift, yfactor, xlinear, xvshift: list of floats
    # xticks_val, yticks_val, xticks_label, yticks_label: list
    # grid_flag, legend_flag: boolean
    # legend_loc: string or int

    if type(spx) is not list:
        print("Parameter \"spx\" is not list.")
        sys.exit()
    if type(spy) is not list:
        print("Parameter \"spy\" is not list.")
        sys.exit()

    Nspec = len(spx)

    colors = colors * math.ceil(Nspec / len(colors)) if len(colors) < Nspec else colors
    linew = linew * math.ceil(Nspec / len(linew)) if len(linew) < Nspec else linew
    lines = lines * math.ceil(Nspec / len(lines)) if len(lines) < Nspec else lines
    labels = labels * math.ceil(Nspec / len(labels)) if len(labels) < Nspec else labels
    yshift = yshift + [0] * (Nspec - len(yshift)) if len(yshift) < Nspec else yshift
    yfactor = yfactor + [1] * (Nspec - len(yfactor)) if len(yfactor) < Nspec else yfactor
    xlshift = xlshift + [0] * (Nspec - len(xlshift)) if len(xlshift) < Nspec else xlshift
    xvshift = xvshift + [0] * (Nspec - len(xvshift)) if len(xvshift) < Nspec else xvshift

    for i in range(Nspec):
        spx_plot = (spx[i] + xlshift[i]) * (1. + xvshift[i] / (scipy.constants.c * 1.e-3))
        spy_plot = spy[i] * yfactor[i] + yshift[i]
        ax.step(spx_plot, spy_plot, where="mid", c=colors[i], label=labels[i], lw=linew[i], linestyle=lines[i])

    if spx_func != "NA" and spy_func != "NA":
        if type(spx_func) is not list:
            print("Parameter \"spx_func\" is not list.")
            sys.exit()
        if type(spy_func) is not list:
            print("Parameter \"spy_func\" is not list.")
            sys.exit()

        Nspec_func = len(spx_func)

        colors_func = colors_func * math.ceil(Nspec_func / len(colors_func)) if len(
            colors_func) < Nspec_func else colors_func
        linew_func = linew_func * math.ceil(Nspec_func / len(linew_func)) if len(
            linew_func) < Nspec_func else linew_func
        lines_func = lines_func * math.ceil(Nspec_func / len(lines_func)) if len(
            lines_func) < Nspec_func else lines_func
        labels_func = labels_func * math.ceil(Nspec_func / len(labels_func)) if len(
            labels_func) < Nspec_func else labels_func

        for i in range(Nspec_func):
            ax.plot(spx_func[i], spy_func[i], c=colors_func[i], label=labels_func[i], lw=linew_func[i],
                    linestyle=lines_func[i])

    if xticks_val != "NA":
        ax.set_xticks(xticks_val)
        if xticks_label != "NA":
            ax.set_xticklabels(xticks_label)
    if yticks_val != "NA":
        ax.set_yticks(yticks_val)
        if yticks_label != "NA":
            ax.set_yticklabels(yticks_label)

    if xaxis_label != "NA":
        ax.set_xlabel(xaxis_label)
    if yaxis_label != "NA":
        ax.set_ylabel(yaxis_label)

    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])

    if grid_flag:
        ax.grid()
    if legend_flag:
        ax.legend(loc=legend_loc)


def transmittance_resampling(spx, spy, telx, tely):
    telx_subp = np.linspace(min(telx), max(telx), 10 * telx.size)

    it = interpolate.interp1d(telx, tely, kind="cubic")
    tely_subp = it(telx_subp)

    spy_tm = np.copy(spy)
    for i in range(len(spx)):
        spy_tm[i] = tely_subp[np.argmin(np.absolute(telx_subp - spx[i]))]

    return spy_tm


def linelist_marker(ax, linelam, linetags, vshift, lammin, lammax, ystart, yend1, yend2, midfactor1=0.3, midfactor2=0.7,
                    fs=9., linew=0.5, linec="k", lettershift=0., itemax=100, crflag="INDEF"):
    linezip = zip(linelam, linetags)
    linezip_sorted = sorted(linezip)
    linelam, linetags = zip(*linezip_sorted)
    linelam = np.array(linelam)
    linetags = np.array(linetags)

    figw = ax.get_window_extent().width
    textwidth = 2. * (lammax - lammin) / 100. * fs / 9. * 900. / figw

    linelam_shifted = linelam * (1. + (vshift * 1.e+3 / scipy.constants.c))

    linelamIN = linelam_shifted[(linelam_shifted > lammin) & (linelam_shifted < lammax)]
    lineID = np.array([i for i in range(linelamIN.size)])
    linetagsIN = linetags[(linelam_shifted > lammin) & (linelam_shifted < lammax)]
    linenum = linelamIN.size

    if linenum == 0:
        return True

    crowded_flag = False
    lam_split = [lammin + (lammax - lammin) / 2. * i for i in range(3)]
    for i in range(2):
        linesplit = linelamIN[(linelamIN > lammin) & (linelamIN < lammax)]
        if linesplit.size * textwidth > 1.3 * (lam_split[i + 1] - lam_split[i]):
            crowded_flag = True

    if type(crflag) == bool:
        crowded_flag = crflag

    # crowded_flag = False
    # if linenum * textwidth > (lammax - lammin):
    #     crowded_flag = True

    LG_lam = [[linelamIN[0]]]
    LG_id = [[lineID[0]]]
    cur_groupid = 0

    for i in range(1, linenum):
        if math.fabs(linelamIN[i] - linelamIN[i - 1]) < textwidth:
            LG_lam[cur_groupid].append(linelamIN[i])
            LG_id[cur_groupid].append(lineID[i])
        else:
            LG_lam.append([linelamIN[i]])
            LG_id.append([lineID[i]])
            cur_groupid += 1

    LG_lamave = np.array([np.average(LG_lam[i]) for i in range(len(LG_lam))])
    LG_size = np.array([len(LG_lam[i]) for i in range(len(LG_lam))])
    LG_size_ave = np.average(LG_size)

    if crowded_flag:
        LG_centerid = []
        LG_width = []
        for i in range(len(LG_lam)):
            if LG_size[i] < LG_size_ave:
                LG_centerid.append((LG_size[i] - 1.) / 2.)
                LG_width.append(LG_size[i])
            else:
                LG_centerid.append((LG_size[i] - 1.) / 4.)
                LG_width.append((1 + LG_size[i]) / 2.)
        LG_centerid = np.array(LG_centerid)
        LG_width = np.array(LG_width)
    else:
        LG_centerid = np.array([(len(LG_id[i]) - 1.) / 2. for i in range(len(LG_lam))])
        LG_width = LG_size

    LG_shift = np.array([0. for i in range(len(LG_lam))])

    LG_left = LG_lamave + LG_shift - LG_width / 2. * textwidth
    LG_right = LG_lamave + LG_shift + LG_width / 2. * textwidth

    LG_left_ol = [max(0., (lammin + textwidth / 2. - LG_left[0]) * 2.)]
    for i in range(1, len(LG_lam)):
        LG_left_ol.append(max(0., LG_right[i - 1] - LG_left[i]))

    LG_right_ol = []
    for i in range(len(LG_lam) - 1):
        LG_right_ol.append(max(0., LG_right[i] - LG_left[i + 1]))
    LG_right_ol.append(max(0., (LG_right[len(LG_lam) - 1] - lammax + textwidth / 2.) * 2.))

    LG_left_ol = np.array(LG_left_ol)
    LG_right_ol = np.array(LG_right_ol)

    # for i in range(len(LG_lam)):
    #     print(LG_lam[i], LG_left_ol[i], LG_right_ol[i])

    nite = 1
    while ((0. in LG_left_ol) or (0. in LG_right_ol)) and (
            (np.sum(LG_left_ol) > textwidth * 1.e-1) or (np.sum(LG_right_ol) > textwidth * 1.e-1)):
        LG_shift += LG_left_ol / 2. - LG_right_ol / 2.
        for i in range(len(LG_lam)):
            if i != 0 and i != len(LG_lam) - 1:
                if LG_left_ol[i] > 0. and LG_right_ol[i] > 0:
                    # print(i,LG_left_ol[i],LG_right_ol[i])
                    LG_shift[i - 1] -= LG_left_ol[i] / 2.
                    LG_shift[i + 1] += LG_right_ol[i] / 2.

        LG_left = LG_lamave + LG_shift - LG_width / 2. * textwidth
        LG_right = LG_lamave + LG_shift + LG_width / 2. * textwidth

        LG_left_ol = [max(0., (lammin + textwidth / 2. - LG_left[0]) * 2.)]
        for i in range(1, len(LG_lam)):
            LG_left_ol.append(max(0., LG_right[i - 1] - LG_left[i]))

        LG_right_ol = []
        for i in range(len(LG_lam) - 1):
            LG_right_ol.append(max(0., LG_right[i] - LG_left[i + 1]))
        LG_right_ol.append(max(0., (LG_right[len(LG_lam) - 1] - lammax + textwidth / 2.) * 2.))

        LG_left_ol = np.array(LG_left_ol)
        LG_right_ol = np.array(LG_right_ol)

        nite += 1
        if nite > itemax:
            break

    LG_plotx = [[] for i in range(len(LG_lam))]
    for i in range(len(LG_lam)):
        if crowded_flag:
            if LG_size[i] < LG_size_ave:
                for j in range(len(LG_lam[i])):
                    LG_plotx[i].append(LG_lamave[i] + (j - LG_centerid[i]) * textwidth + LG_shift[i])
            else:
                for j in range(len(LG_lam[i])):
                    # print(i,j,(j/2. - LG_centerid[i]) * textwidth + LG_shift[i], LG_shift[i])
                    LG_plotx[i].append(LG_lamave[i] + (j / 2. - LG_centerid[i]) * textwidth + LG_shift[i])
        else:
            for j in range(len(LG_lam[i])):
                LG_plotx[i].append(LG_lamave[i] + (j - LG_centerid[i]) * textwidth + LG_shift[i])

    crowded_y = [yend1, yend2]
    ymid1 = (yend1 - ystart) * midfactor1 + ystart
    ymid2 = (yend1 - ystart) * midfactor2 + ystart
    ytext = (yend1 - ystart) * 0.05 + lettershift

    for i in range(len(LG_lam)):
        if crowded_flag:
            if LG_size[i] < LG_size_ave:
                for j in range(len(LG_lam[i])):
                    plt.plot([LG_lam[i][j], LG_lam[i][j], LG_plotx[i][j], LG_plotx[i][j]],
                             [ystart, ymid1, ymid2, crowded_y[0]], color=linec, lw=linew)
                    plt.text(LG_plotx[i][j], crowded_y[0] + ytext, linetagsIN[LG_id[i][j]], rotation=90, ha="center",
                             va="bottom", fontsize=fs)
            else:
                for j in range(len(LG_lam[i])):
                    plt.plot([LG_lam[i][j], LG_lam[i][j], LG_plotx[i][j], LG_plotx[i][j]],
                             [ystart, ymid1, ymid2, crowded_y[j % 2]], color=linec, lw=linew)
                    plt.text(LG_plotx[i][j], crowded_y[j % 2] + ytext, linetagsIN[LG_id[i][j]], rotation=90, ha="center",
                             va="bottom", fontsize=fs)
        else:
            for j in range(len(LG_lam[i])):
                plt.plot([LG_lam[i][j], LG_lam[i][j], LG_plotx[i][j], LG_plotx[i][j]],
                         [ystart, ymid1, ymid2, crowded_y[0]], color=linec, lw=linew)
                plt.text(LG_plotx[i][j], crowded_y[0] + ytext, linetagsIN[LG_id[i][j]], rotation=90, ha="center",
                         va="bottom", fontsize=fs)

    return crowded_flag


def ShadeAFregions(ax1, AFstart, AFend, slitpos, miny, maxy, color="0.5", alpha=0.5, lettercolor="k", lettersize="12"):
    for i in range(len(AFstart)):
        ls, le, sp = AFstart[i], AFend[i], slitpos[i]
        ax1.fill([ls, ls, le, le, ls], [miny, maxy, maxy, miny, miny], color=color, alpha=alpha)
        ax1.text((ls + le) / 2., miny, sp, va="bottom", ha="center", color=lettercolor, fontsize=lettersize)


def MultiSpecPlotter(spxlist, spylist, pp, colors, righttexts, titletext, fs=9., obsdates="INDEF", order="INDEF",
                     v_helio=0., telx="INDEF", tely="INDEF", plotsig=5., ite=20, vacorair="AIR", telflag=True,
                     DIBplot=True, STEplot=True, C2plot=True, sptype="OB", zerooffset=False, model="INDEF", DIBgrid=True):
    spnum = len(spxlist)

    if telflag:
        if telx == "INDEF" or tely == "INDEF":
            telflag = False
            print("Telluric spectra is not inputted.")

    if telflag:
        ax1 = plt.axes([0.04, 0.2, 0.85, 0.7])
        lammin, lammax, ystep = MultiSpecAxis(ax1, spxlist, spylist, colors, xaxis_label="NA",
                                              yaxis_label=r"Normalized flux + offset", obsdates=obsdates, order=order,
                                              v_helio=v_helio, plotsig=plotsig, ite=ite, vacorair=vacorair,
                                              DIBplot=DIBplot, STEplot=STEplot, C2plot=C2plot, sptype=sptype,
                                              zerooffset=zerooffset, model=model, DIBgrid=DIBgrid)
    else:
        ax1 = plt.axes([0.04, 0.1, 0.85, 0.8])
        lammin, lammax, ystep = MultiSpecAxis(ax1, spxlist, spylist, colors,
                                              xaxis_label=r"Wavelength %s ($\AA$)" % vacorair,
                                              yaxis_label=r"Normalized flux + offset", obsdates=obsdates, order=order,
                                              v_helio=v_helio, plotsig=plotsig, ite=ite, vacorair=vacorair,
                                              DIBplot=DIBplot, STEplot=STEplot, C2plot=C2plot, sptype=sptype,
                                              zerooffset=zerooffset, model=model, DIBgrid=DIBgrid)

    plt.title(titletext)

    if zerooffset:
        for i in range(spnum):
            plt.text(lammax, 1. + ystep * (spnum - i - 1) * 0.1, righttexts[i], color=colors[i], fontsize=fs)
    else:
        for i in range(spnum):
            plt.text(lammax, 1. + ystep * (spnum - i - 1), righttexts[i], color=colors[i], fontsize=fs)

    if telflag:
        ax2 = plt.axes([0.04, 0.1, 0.85, 0.1])
        spectrum_plot(ax2, [telx], [tely], [lammin, lammax], [0., 1.2], xaxis_label=r"Wavelength %s ($\AA$)" % vacorair,
                      yaxis_label="Transmittance")

    plt.savefig(pp, format="pdf")
    plt.clf()


def MultiSpecAxis(ax1, spxlist, spylist, colors, xaxis_label="NA", yaxis_label="NA", obsdates="INDEF", order="INDEF",
                  v_helio=0., plotsig=5.0, ite=20, vacorair="AIR", DIBplot=True, STEplot=True, C2plot=True, sptype="OB",
                  zerooffset=True, model="INDEF", DIBgrid=True):
    spnum = len(spxlist)

    if obsdates != "INDEF":
        if type(obsdates) == list:
            if type(v_helio) != list:
                print("'v_helio' must be list if 'obsdates' is list.")
                sys.exit()
            elif len(v_helio) != len(obsdates):
                print("The numbers of components of 'v_helio' and 'obsdates' are different.")
                sys.exit()

    lammin, lammax, ymin, ymax, yrange, ymed, ystd = [], [], [], [], [], [], []
    spylist_norm = []

    for i in range(spnum):
        ymed1 = np.median(spylist[i])
        ystd1 = np.std(spylist[i])
        yclip = np.absolute(spylist[i] - ymed1) < plotsig * ystd1
        for k in range(ite):
            ystdprev = copy.copy(ystd1)
            ymed1 = np.median(spylist[i][yclip])
            ystd1 = np.std(spylist[i][yclip])
            if ystdprev == ystd1:
                break
            yclip = np.absolute(spylist[i] - ymed1) < plotsig * ystd1

        ymed.append(ymed1)
        ystd.append(ystd1 / ymed1)
        lammin.append(np.amin(spxlist[i]))
        lammax.append(np.amax(spxlist[i]))
        ymin.append(np.amin(spylist[i] / ymed1))
        ymax.append(np.amax(spylist[i] / ymed1))
        yrange.append(ymax[i] - ymin[i])
        spylist_norm.append(spylist[i] / ymed[i])

    ystdave = np.average(ystd)
    ystep = 0.5 if zerooffset else max(ystdave * plotsig, 0.05)
    ylimlow = max(min(ymin), 1.0 - ystep * 3.) - 0.1
    ylimupp = 2. - ylimlow if zerooffset else 1. + ystep * (spnum - 1.)

    DIBSTEflag = False
    wavDIBSTE = np.array([])
    DIBSTEtags = np.array([])

    if DIBplot:
        DIBreturn = GetDIBList(min(lammin), max(lammax), vacorair=vacorair)
        if DIBreturn != None:
            DIBSTEflag = True
            [DIBID, wavDIB, refDIB, cateDIB, _, _, _] = DIBreturn
            DIBtags = np.array(["DIB%d" % (wavDIB[k]) for k in range(len(DIBID))])
            wavDIBSTE = np.concatenate((wavDIBSTE, wavDIB))
            DIBSTEtags = np.concatenate((DIBSTEtags, DIBtags))
            # linelist_marker(ax1, wavDIB, DIBtags, 0., min(lammin), max(lammax), 1.0 + ystep * (spnum - 0.8),
            #                 1.0 + ystep * (spnum - 0.4), 1.0 + ystep * (spnum + 0.5))

    if STEplot:
        STEreturn = GetLineList(min(lammin), max(lammax), sptype=sptype, vacorair=vacorair)
        if STEreturn != None:
            DIBSTEflag = True
            [StellarID, wavSte, atomSte, ionSte, sptypeSte] = STEreturn
            STEtags = np.array(["%s%s %.1f" % (atomSte[k], ionSte[k], wavSte[k]) for k in range(len(StellarID))])
            wavDIBSTE = np.concatenate((wavDIBSTE, wavSte))
            DIBSTEtags = np.concatenate((DIBSTEtags, STEtags))
            # linelist_marker(ax1, wavSte, STEtags, 0., min(lammin), max(lammax), 1.0 + ystep * (spnum + 0.2),
            #                 1.0 + ystep * (spnum + 0.6), 1.0 + ystep * (spnum + 1.5))

    if DIBSTEflag:
        if not zerooffset:
            crw = linelist_marker(ax1, wavDIBSTE, DIBSTEtags, 0., min(lammin), max(lammax), ylimupp + 0.2 * ystep * spnum ** 0.5,
                                  ylimupp + 0.4 * ystep * spnum ** 0.5, ylimupp + 0.8 * ystep * spnum ** 0.5)
            if crw:
                ylimupp += 1.2 * ystep * spnum ** 0.5
            else:
                ylimupp += 0.8 * ystep * spnum ** 0.5
        else:
            crw = linelist_marker(ax1, wavDIBSTE, DIBSTEtags, 0., min(lammin), max(lammax), 1. + (1.-ylimlow) * 0.2,
                                  1. + (1.-ylimlow) * 0.4, 1. + (1.-ylimlow) * 0.8)

    C2flag = False

    if C2plot:
        C2params = moleclinelist("C2")
        bandnames = ["AX00", "AX10", "AX20"]
        C2bands = [C2params.bandlist[i] for i in bandnames]
        C2range = [[i.minimumlambda(), i.maximumlambda()] for i in C2bands]

        C2y = ylimupp
        for k in range(3):
            if min(lammin) < C2range[k][1] and max(lammax) > C2range[k][0]:
                C2flag = True
                if zerooffset:
                    C2marks(ax1, bandnames[k], 0., 0.,
                            ycoord=[C2y + 0.85 * ystep * spnum ** 0.5, C2y + 0.55 * ystep * spnum ** 0.5, C2y + 0.25 * ystep * spnum ** 0.5, C2y])
                else:
                    C2marks(ax1, bandnames[k], 0., 0.,
                            ycoord=[C2y + 0.85 * ystep * spnum ** 0.5, C2y + 0.55 * ystep * spnum ** 0.5, C2y + 0.25 * ystep * spnum ** 0.5, C2y])

    if C2flag:
        if not zerooffset:
            ylimupp += 1.15 * ystep * spnum ** 0.5

    if not DIBSTEflag:
        if not C2flag:
            if not zerooffset:
                ylimupp += 0.3 * ystep * spnum ** 0.5

    if obsdates != "INDEF" and order != "INDEF":
        AFreturn = GetAFList(min(lammin), max(lammax), obsdates, order, v_helio, vacorair=vacorair)
        if AFreturn != None:
            ShadeAFregions(ax1, AFreturn[0], AFreturn[1], AFreturn[2], ylimlow,
                           ylimupp + 0.2)

    if model != "INDEF":
        modelx1, modely1 = model_opener(model)
        modely = [modely1]
        if vacorair == "AIR":
            modelx = [modelx1]
        elif vacorair == "VAC":
            modelx = [air2vac(modelx1)]
    else:
        modelx, modely = "NA", "NA"

    if DIBgrid:
        if DIBreturn != None:
            for i in wavDIB:
                ax1.plot([i,i], [ylimlow, ylimupp+0.2], "y")

    if zerooffset:
        spectrum_plot(ax1, spxlist, spylist_norm, [min(lammin), max(lammax)],
                      [ylimlow, ylimupp], colors=colors,
                      yshift=[0. for k in range(spnum)], xaxis_label=xaxis_label,
                      yaxis_label=yaxis_label, spx_func=modelx, spy_func=modely)

    else:
        spectrum_plot(ax1, spxlist, spylist_norm, [min(lammin), max(lammax)],
                      [ylimlow, ylimupp+0.2], colors=colors,
                      yshift=[ystep * (spnum - k - 1) for k in range(spnum)], xaxis_label=xaxis_label,
                      yaxis_label=yaxis_label, spx_func=modelx, spy_func=modely)

    return min(lammin), max(lammax), ystep


if __name__ == "__main__":
    fig = plt.figure(figsize=(10, 6))

    # [origin_x, origin_y, width, height]
    pp = PdfPages("/Users/hamano/DIB_analysis/TMP_telluric_dirs/29Vul_140913_v2.pdf")

    # open spectrum file
    plpath = "/Users/hamano/DIB_analysis/Standard_pipeline_dir/2014_09_13/2014_09_13-29Vul_pipeline_ver3.5/"
    advpath = "/Users/hamano/DIB_analysis/TMP_telluric_dirs/29Vul_140913_v2/telluric/FITS/fsr1.30/"
    plspfiles = glob.glob(plpath + "29Vul_sum/VAC_norm/fsr1.30/*norm.fits")
    plspfiles.sort()
    advfiles = glob.glob(advpath + "*fits")
    advfiles.sort()

    for i in range(len(advfiles)):
        ax1 = plt.axes([0.08, 0.1, 0.9, 0.85])

        spx1, spy1, _, _, _ = openspecfits(plspfiles[i])
        spx2, spy2, _, _, _ = openspecfits(advfiles[i])

        spectrum_plot(ax1, [spx1, spx2], [spy1, spy2], [min(spx1), max(spx2)], [0.0, 1.2], colors=["b", "k"],
                      legend_flag=True, labels=["orig.", "std. lines removed"])

        plt.savefig(pp, format="pdf")
        plt.clf()

    pp.close()

    #
    # # plot spectra
    # spectrum_plot(ax1, [spx], [spy], [12100, 12300], [0.5, 1.2], xaxis_label="X jiku", yaxis_label="Y jiku ($\mu$)")
    #
    # spectrum_plot(ax2, [spx, spx], [spy, spy], [12100, 12120], [0.5, 1.2], colors=["y", "r"], lines=["--", "-"],
    #               linew=[3., 1.], legend_flag=True, labels=["base", "(x,y) (+3,+0.1)"], xlshift=[0., 3.], yshift=[0., 0.1],
    #               grid_flag=True)
    #
    # spectrum_plot(ax3, [spx, spx, spx], [spy, spy, spy], [12050, 12090], [0.9, 1.1], xvshift=[0., -40., 50.],
    #               yshift=[-0.01, 0., 0.01], labels=["base", "-40 km/s shift", "50 km/s shift"], legend_flag=True,
    #               xticks_val=[12050, 12080], yticks_val=[0.8, 1.0], legend_loc="upper left")
    #
    #
    # plt.savefig("test.png")
