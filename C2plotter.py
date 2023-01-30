# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy.constants
from Spec1Dtools import openspecfits
from C2_parameters import *
from C2param_parser import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]

def C2marks(ax, band, vshift, textxshift, textyshift=0.003, ycoord=[1.11, 1.08, 1.05, 1.025], marklen=0.01, linew=1.,
            mcolor="k"):
    # ycoord: P, Q, Rlong, Rshort

    C2band = [C2params.bandlist[band].branchlist("P"),
              C2params.bandlist[band].branchlist("Q"),
              C2params.bandlist[band].branchlist("R")]

    vf = (1. + vshift / (scipy.constants.c * 1.e-3))
    theta_array = np.arange(0, 1, 0.01) * math.pi + math.pi / 2.

    [lowlim, upplim] = ax.get_xlim()

    for i in range(3):
        lamtmp = []
        if i != 2:
            for line in C2band[i]:
                if lowlim < line.lamair < upplim:
                    plt.plot([line.lamair * vf, line.lamair * vf], [ycoord[i], ycoord[i] + marklen],
                             lw=linew, color=mcolor)
                    plt.text(line.lamair * vf + textxshift, ycoord[i] + marklen + textyshift, line.rotJ,
                             fontname="Arial", fontsize=10, va="bottom", ha="center")
                lamtmp.append(line.lamair * vf)
            plt.plot([min(lamtmp), max(lamtmp)], [ycoord[i] + marklen, ycoord[i] + marklen], lw=linew, color=mcolor)
        else:
            for line in C2band[i]:
                if lowlim < line.lamair < upplim:
                    if line.rotJ <= 6:
                        plt.plot([line.lamair * vf, line.lamair * vf],
                                 [ycoord[i + 1], ycoord[i + 1] + marklen], lw=linew, color=mcolor)
                        plt.text(line.lamair * vf + textxshift, ycoord[i + 1] + marklen + textyshift,
                                 line.rotJ,
                                 fontname="Arial", fontsize=10, va="bottom", ha="center")
                    else:
                        plt.plot([line.lamair * vf, line.lamair * vf], [ycoord[i], ycoord[i] + marklen],
                                 lw=linew, color=mcolor)
                        plt.text(line.lamair * vf + textxshift, ycoord[i] + marklen + textyshift, line.rotJ,
                                 fontname="Arial", fontsize=10, va="bottom", ha="center")
                lamtmp.append(line.lamair * vf)

            plt.plot([min(lamtmp), max(lamtmp)], [ycoord[i] + marklen, ycoord[i] + marklen], lw=linew,
                     color=mcolor)
            plt.plot([min(lamtmp), lamtmp[0]], [ycoord[i + 1] + marklen, ycoord[i + 1] + marklen], lw=linew,
                     color=mcolor)
            plt.plot(3. * np.cos(theta_array) + min(lamtmp),
                     0.5 * (ycoord[2] - ycoord[3]) * np.sin(theta_array) + (ycoord[2] + ycoord[3]) / 2. + marklen,
                     color=mcolor, lw=linew)


if __name__ == "__main__":

    filename = sys.argv[1:]

    [spfile, telfile, outputfig, band] = filename[0:4]
    if len(filename) > 4:
        options = filename[4:]
    else:
        options = []

    textshift = 0.

    oneflag = False
    if "one" in options:
        oneflag = True

    if band not in bandlist:
        print("Registered bands keywords are: ")
        for i in bandlist:
            print(i)
        sys.exit()

    xlim = {bandlist[0]: [12060, 12180, 12300], bandlist[1]: [10120, 10260, 10400]}
    ylim = [0.75, 1.15]
    thres = 0.6

    spx, spy, _, _, _ = openspecfits(spfile)
    telx, tely, _, _, _ = openspecfits(telfile)
    spy_tm = transmittance_resampling(spx, spy, telx, tely)
    spy_mask = np.ma.masked_where(spy_tm < thres, spy)

    if "setting" in options:
        config = ReadConfig(band=band)
        config.set("Spectrum", "target", spfile)
        config.set("Spectrum", "telluric", telfile)
        config.set("C2 parameters", "band", band)
        OverwriteConfig(config)

    figw = 0.87
    figh_tel = 0.08
    figh_obj = 0.3
    orig_w = 0.08
    orig_h = 0.05
    gap_telobj = 0.025
    gap_objtel = 0.1

    if outputfig.split(".")[-1] == "pdf":
        pp = PdfPages(outputfig)

    if oneflag:
        fig = plt.figure(figsize=(10, 5))
        ax1 = plt.axes([orig_w, orig_h * 2., figw, figh_tel * 2.])
        spectrum_plot(ax1, [telx], [tely], [xlim[band][0], xlim[band][1]], [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                      yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

        ax2 = plt.axes([orig_w, (orig_h + figh_tel + gap_telobj) * 2., figw, figh_obj * 2.])
        spectrum_plot(ax2, [spx, spx], [spy, spy_mask], [xlim[band][0], xlim[band][1]], ylim, colors=["0.8", "k"],
                      xticks_val=range(xlim[band][0], xlim[band][1], 20),
                      xticks_label=["" for i in range(xlim[band][0], xlim[band][1], 20)],
                      yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1])
        C2marks(ax2, band, 0., textshift)

    else:
        fig = plt.figure(figsize=(10, 10))
        ax1 = plt.axes([orig_w, orig_h, figw, figh_tel])
        spectrum_plot(ax1, [telx], [tely], [xlim[band][1], xlim[band][2]], [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                      yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

        ax2 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj, figw, figh_obj])
        spectrum_plot(ax2, [spx, spx], [spy, spy_mask], [xlim[band][1], xlim[band][2]], ylim, colors=["0.8", "k"],
                      xticks_val=range(xlim[band][1], xlim[band][2], 20),
                      xticks_label=["" for i in range(xlim[band][1], xlim[band][2], 20)],
                      yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1])
        C2marks(ax2, band, 0., textshift)

        ax3 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj + gap_objtel + figh_obj, figw, figh_tel])
        spectrum_plot(ax3, [telx], [tely], [xlim[band][0], xlim[band][1]], [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                      yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

        ax4 = plt.axes(
            [orig_w, orig_h + figh_tel + gap_telobj + gap_objtel + figh_tel + figh_obj + gap_telobj, figw, figh_obj])
        spectrum_plot(ax4, [spx, spx], [spy, spy_mask], [xlim[band][0], xlim[band][1]], ylim, colors=["0.8", "k"],
                      xticks_val=range(xlim[band][0], xlim[band][1], 20),
                      xticks_label=["" for i in range(xlim[band][0], xlim[band][1], 20)],
                      yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1])
        C2marks(ax4, band, 0., textshift)

    if outputfig.split(".")[-1] == "pdf":
        plt.savefig(pp, format="pdf")
        pp.close()
    else:
        plt.savefig(outputfig)
