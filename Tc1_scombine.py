#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sys
import matplotlib.pyplot as plt
import glob
from combine_MySQL import pycombine
from combine_normalize import normalizeCombinedspec
from DIBanalysis import maskRegion_line, continuum

if __name__ == '__main__':
    ppdir = "/Users/hamano/DIB_analysis/DIB_pipeline_dir/HD1994786/2017_07_31-Tc_1_pipeline_ver3.7_manual/"
    frnum = 12
    spdirs = [ppdir + "Tc_1_NO{}/AIR_flux/fsr1.30/Tc_1_NO{}_m53_fsr1.30_AIR.fits".format(i+1, i+1) for i in range(frnum)]
    outputdir = "/Users/hamano/DIB_analysis/DIB_pipeline_dir/HD1994786/2017_07_31-Tc_1_pipeline_ver3.7_manual/Tc_1_small_data/"

    for i in range(2,7,2):
        for j in range(int(frnum/i)):
            outputfits = outputdir + "Tc_1_{}fr_{}_m53_fsr1.30_AIR.fits".format(i,j)
            normfits = outputdir + "Tc_1_{}fr_{}_m53_fsr1.30_AIRnorm.fits".format(i,j)
            pycombine(spdirs[i*j:i*(j+1)], outputfits, [1.0 for _ in range(i)])
            maskHeI = maskRegion_line(outputfits, 10830., 0.)
            continuum(outputfits, normfits, order=15, high_rej=2., low_rej=2.)


