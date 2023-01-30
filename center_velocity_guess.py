#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import scipy.constants
from add_lineDIB_mysql import GetDIBLine
import numpy

if __name__ == "__main__":
    refids = [35, 37, 38]
    vacorair = "AIR"
    inputwav = float(sys.argv[1])
    refwav = []
    for id in refids:
        DIBinfo = GetDIBLine(id, vacorair=vacorair)
        refwav.append(DIBinfo[1])
    refwav = numpy.array(refwav)
    choosedwav = refwav[numpy.argmin(numpy.absolute(refwav - inputwav))]
    cvel = (inputwav - choosedwav) / choosedwav * scipy.constants.c * 1.e-3

    print("Wavelength system: %s" % vacorair)
    print("Input wavelength: %.3f" % inputwav)
    print("Reference wavelength: %.3f" % choosedwav)
    print("------------------------------")
    print("Velocity: %.2f" % cvel)
