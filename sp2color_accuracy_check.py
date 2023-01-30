#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
from sp2color import sp2color
import matplotlib.pyplot as plt
import numpy

if __name__ == '__main__':
    rf = open("temporaly_files/EBV_Friedman2011_formatted.txt")
    rl = rf.readlines()
    rf.close()

    ldic = {"I":1, "II":2, "III":3, "IV":4, "V":5}
    colors = {1:"b", 2:"g", 3:"cyan", 4:"k", 5:"r"}

    id = numpy.array([int(l.split("\t")[0]) for l in rl[1:]])
    name = numpy.array([l.split("\t")[1] for l in rl[1:]])
    type = numpy.array([l.split("\t")[2] for l in rl[1:]])
    subtype = numpy.array([float(l.split("\t")[3]) for l in rl[1:]])
    lclass = numpy.array([ldic[l.split("\t")[4]] for l in rl[1:]])
    bvobs = numpy.array([float(l.split("\t")[5]) for l in rl[1:]])
    ebv_F11 = numpy.array([float(l.split("\t")[6]) for l in rl[1:]])

    bvcal = numpy.array([sp2color(type[i], subtype[i], lclass[i]) for i in range(len(id))])

    plt.figure()
    for i in range(1,6):
        plt.scatter(bvobs[lclass==i] - bvcal[lclass==i], bvobs[lclass==i] - bvcal[lclass==i] - ebv_F11[lclass==i], color=colors[i], label=str(i))

    plt.legend()
    plt.xlabel("B-V (simbad) - B-V (calc.)")
    plt.ylabel("E(B-V) - E(B-V) (Friedman et al. 2011")
    plt.grid()
    plt.savefig("temporaly_files/EBV_Friedman2011_formatted2.png")
    plt.clf()
