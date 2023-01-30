#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import numpy
import matplotlib.pyplot as plt

if __name__ == '__main__':
    araki = open("WINERED_ad_Araki.dat", "r")
    rl = araki.readlines()
    araki.close()

    pp_ara = [i.split()[0] for i in rl]
    fnum_ara = [i.split()[1] for i in rl]
    am_ara = numpy.array([float(i.split()[2]) for i in rl])
    a_ara = numpy.array([float(i.split()[3]) for i in rl])
    b_ara = numpy.array([float(i.split()[4]) for i in rl])
    std_ara = numpy.array([float(i.split()[5]) for i in rl])

    ntt = open("WINERED_ad_NTT.dat", "r")
    rl = ntt.readlines()
    ntt.close()
    pp_ntt = [i.split()[0] for i in rl]
    fnum_ntt = [i.split()[1] for i in rl]
    am_ntt = numpy.array([float(i.split()[2]) for i in rl])
    a_ntt = numpy.array([float(i.split()[3]) for i in rl])
    b_ntt = numpy.array([float(i.split()[4]) for i in rl])
    std_ntt = numpy.array([float(i.split()[5]) for i in rl])

    plt.figure()
    plt.scatter(am_ara, a_ara * 1000, label="Araki")
    plt.scatter(am_ntt, a_ntt * 1000, label="NTT")
    plt.legend()
    plt.xlabel("Airmass")
    plt.ylabel("Peak shift/Wavelength [pix/1000A]")
    plt.xlim(1., 2.)
    plt.grid()
    plt.title("WINERED, WIDE-mode, S/N>200")
    plt.savefig("winered_atmospheric_dispersion.png", format="png")
    plt.clf()

    am_ave = [1.025 + i * 0.05 for i in range(20)]
    aave_ara = []
    aave_ntt = []
    amed_ara = []
    amed_ntt = []
    astd_ara = []
    astd_ntt = []
    for i in range(20):
        aave_ara.append(numpy.average(a_ara[numpy.logical_and(am_ara > am_ave[i] - 0.25, am_ara < am_ave[i] + 0.25)] * 1000))
        aave_ntt.append(numpy.average(a_ntt[numpy.logical_and(am_ntt > am_ave[i] - 0.25, am_ntt < am_ave[i] + 0.25)] * 1000))
        amed_ara.append(numpy.median(a_ara[numpy.logical_and(am_ara > am_ave[i] - 0.25, am_ara < am_ave[i] + 0.25)] * 1000))
        amed_ntt.append(numpy.median(a_ntt[numpy.logical_and(am_ntt > am_ave[i] - 0.25, am_ntt < am_ave[i] + 0.25)] * 1000))
        astd_ara.append(numpy.std(a_ara[numpy.logical_and(am_ara > am_ave[i] - 0.25, am_ara < am_ave[i] + 0.25)] * 1000))
        astd_ntt.append(numpy.std(a_ntt[numpy.logical_and(am_ntt > am_ave[i] - 0.25, am_ntt < am_ave[i] + 0.25)] * 1000))

    plt.scatter(am_ave, aave_ara, label="Araki")
    plt.scatter(am_ave, aave_ntt, label="NTT")
    plt.legend()
    plt.xlabel("Airmass")
    plt.ylabel("Peak shift/Wavelength [pix/1000A]")
    plt.xlim(1., 2.)
    plt.grid()
    plt.title("WINERED, WIDE-mode, S/N>200")
    plt.savefig("winered_atmospheric_dispersion_average.png", format="png")
    plt.clf()

    plt.scatter(am_ave, amed_ara, label="Araki")
    plt.scatter(am_ave, amed_ntt, label="NTT")
    plt.legend()
    plt.xlabel("Airmass")
    plt.ylabel("Peak shift/Wavelength [pix/1000A]")
    plt.xlim(1., 2.)
    plt.grid()
    plt.title("WINERED, WIDE-mode, S/N>200")
    plt.savefig("winered_atmospheric_dispersion_median.png", format="png")
    plt.clf()

    plt.scatter(am_ave, astd_ara, label="Araki")
    plt.scatter(am_ave, astd_ntt, label="NTT")
    plt.legend()
    plt.xlabel("Airmass")
    plt.ylabel("Peak shift/Wavelength [pix/1000A]")
    plt.xlim(1., 2.)
    plt.grid()
    plt.title("WINERED, WIDE-mode, S/N>200")
    plt.savefig("winered_atmospheric_dispersion_std.png", format="png")
    plt.clf()

