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

    hobbslist = GetDIBList(1., 1.e+5, reference="Hobbs et al. (2009), ApJ, 705, 32", vacorair="AIR")
    [DIBIDs, wav_air_db, reference_db, category_db, fwhm_db, wavenumber_db, comment_db] = hobbslist
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
            if wave_prev[i].find("/") == -1:
                newdib.append(ids[i])
                if waveerr[i] != "    ":
                    AddDIBDB(float(wave[i]), float(fwhm[i]) / vl * float(wave[i]), reference, "weak",
                             comment[i])
                    print("Add: ", float(wave[i]), float(fwhm[i]) / vl * float(wave[i]), reference, "weak",
                          comment[i])
            elif wave_prev[i].split("/")[0].find("*") == -1:
                hd183143 = float(wave_prev[i].split("/")[0])
                if hd183143 in wav_air_db:
                    updatedib.append(ids[i])
                    if waveerr[i] != "    ":
                        UpdateDIBDB(DIBIDs[wav_air_db.index(hd183143)], lam=float(wave[i]),
                                    fwhm=float(fwhm[i]) / vl * float(wave[i]), reference=reference,
                                    comment=comment[i], confirm=False)
                        print("Update: DIBID%d" % (DIBIDs[wav_air_db.index(hd183143)]), float(wave[i]),
                              float(fwhm[i]) / vl * float(wave[i]), reference, "weak", comment[i])
                else:
                    insertdib.append(ids[i])
                    if waveerr[i] != "    ":
                        AddDIBDB(float(wave[i]), float(fwhm[i]) / vl * float(wave[i]), reference, "weak",
                                 comment[i])
                        print("Add: ", float(wave[i]), float(fwhm[i]) / vl * float(wave[i]), reference, "weak",
                              comment[i])
            else:
                insertdib.append(ids[i])
                if waveerr[i] != "    ":
                    AddDIBDB(float(wave[i]), float(fwhm[i]) / vl * float(wave[i]), reference, "weak",
                             comment[i])
                    print("Add: ", float(wave[i]), float(fwhm[i]) / vl * float(wave[i]), reference, "weak",
                          comment[i])

    print(len(DIBIDs))
    print(newdib, len(newdib))
    print(updatedib, len(updatedib))
    print(insertdib, len(insertdib))
    print(len(newdib) + len(updatedib) + len(insertdib))
