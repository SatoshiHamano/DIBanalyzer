#!/usr/bin/env python
# -*- coding:utf-8 -*-

from open_mysql_project import openproject
from add_lineDIB_mysql import GetDIBList, UpdateDIBDB, AddDIBDB
import sys
import scipy.constants

if __name__ == "__main__":
    rf = open("analysis/linelist/FanHaoyu_apjab1b74t2_mrt.txt", "r")
    rl = rf.readlines()
    rf.close()

    vl = scipy.constants.c * 1.e-3

    optlist = GetDIBList(4000., 9000., vacorair="AIR")
    [DIBIDs, wav_air_db, reference_db, category_db, fwhm_db, wavenumber_db, comment_db] = optlist
    wav_air_db = list(wav_air_db)

    ids = [rl[i][0:3] for i in range(43, len(rl))]
    wave = [rl[i][4:13] for i in range(43, len(rl))]
    waveerr = [rl[i][14:18] for i in range(43, len(rl))]
    velerr = [rl[i][19:24] for i in range(43, len(rl))]
    wave_prev = [rl[i][25:40] for i in range(43, len(rl))]
    fwhm = [rl[i][41:48] for i in range(43, len(rl))]
    fwhm_err = [rl[i][49:55] for i in range(43, len(rl))]
    fwhm_prev = [rl[i][56:69] for i in range(43, len(rl))]
    xm2 = [rl[i][70:76] for i in range(43, len(rl))]
    xm2err = [rl[i][77:82] for i in range(43, len(rl))]
    ebvratio = [rl[i][83:90] for i in range(43, len(rl))]
    dlr = [rl[i][91:98] for i in range(43, len(rl))]
    detpct = [rl[i][99:112] for i in range(43, len(rl))]
    comment = [rl[i][113:176] for i in range(43, len(rl))]

    reference = "Fan et al. (2019), ApJ, 878, 151"

    newdib = []
    updatedib = []
    insertdib = []
    registereddib = []
    for i in range(len(ids)):
        if ids[i] != "   ":
            if wave_prev[i].find("/") != -1:
                if wave_prev[i].split("/")[1].find("*") == -1:
                    hd204827 = wave_prev[i].split("/")[1]
                    if waveerr[i] != "    ":
                        print("%d\t%s" % (DIBIDs[wav_air_db.index(float(wave[i]))], hd204827))

