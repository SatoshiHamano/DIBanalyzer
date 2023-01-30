#!/usr/bin/env python
# -*- coding:utf-8 -*-

import mysql.connector
from urllib.parse import urlparse
import sys
from vac2air_spec import *
import numpy
from open_mysql_project import openproject
import scipy.constants


def ionstate_converter(ionint):
    if ionint == 1:
        return "I"
    elif ionint == 2:
        return "II"
    elif ionint == 3:
        return "III"
    elif ionint == 4:
        return "IV"
    elif ionint == 5:
        return "V"
    else:
        return str(ionint)


def ionstate_converter_reverse(ionstr):
    if ionstr == "I":
        return 1
    elif ionstr == "II":
        return 2
    elif ionstr == "III":
        return 3
    elif ionstr == "IV":
        return 4
    elif ionstr == "V":
        return 5
    else:
        return ionstr


def AddLineDB(lam, atom, ionstate, sptype, reference, vacorair="AIR"):
    conn, cur = openproject()
    cur.execute("select lineID,wavelength_vac,wavelength_air,atom,ionstate,sptype from StellarLine;")
    rows = cur.fetchall()
    if rows == []:
        curID = 1
        newflag = True
    else:
        prevIDs = [int(i[0]) for i in rows]
        wav_vac_prev = numpy.array([float(i[1]) for i in rows])
        wav_air_prev = numpy.array([float(i[2]) for i in rows])
        atom_prev = [i[3] for i in rows]
        ionstate_prev = [int(i[4]) for i in rows]
        sptype_prev = [i[5] for i in rows]
        curID = max(prevIDs) + 1
        newflag = False

    if type(lam) is list:
        for i in range(len(lam)):
            if vacorair == "VAC":
                lamvac, lamair = lam[i], vac2air(lam[i])
            else:
                lamvac, lamair = air2vac(lam[i]), lam[i]
            if not newflag:
                if numpy.amin(numpy.absolute(wav_vac_prev - lamvac)) < 0.01:
                    indx = numpy.argmin(numpy.absolute(wav_vac_prev - lamvac))
                    if atom[i] == atom_prev[indx] and ionstate[i] == ionstate_prev[indx]:
                        if sptype[i] in sptype_prev[indx]:
                            print("%s%s (%.2f) is already registered as lineID=%d." % (
                                atom[i], ionstate_converter(ionstate[i]), lam[i], prevIDs[indx]))
                            continue
                        else:
                            cur.execute(
                                "UPDATE StellarLine SET sptype='%s' where lineID=%d;" % (
                                    sptype_prev[indx] + sptype[i], prevIDs[indx])
                            )
                            print("Sptype is updated from '%s' to '%s' for lineID=%d." % (
                                sptype_prev[indx], sptype_prev[indx] + sptype[i], prevIDs[indx]))

                cur.execute(
                    "INSERT IGNORE INTO StellarLine (lineID, wavelength_vac,wavelength_air,atom,ionstate,sptype,reference) VALUES(%d,%.3f,%.3f,'%s',%d,'%s','%s');" % (
                        curID, lamvac, lamair, atom[i], ionstate[i], sptype[i], reference[i]))
                curID += 1

    else:
        if vacorair == "VAC":
            lamvac, lamair = lam, vac2air(lam)
        else:
            lamvac, lamair = air2vac(lam), lam

        if not newflag:
            if numpy.amin(numpy.absolute(wav_vac_prev - lamvac)) < 0.01:
                indx = numpy.argmin(numpy.absolute(wav_vac_prev - lamvac))
                if atom == atom_prev[indx] and ionstate == ionstate_prev[indx]:
                    if sptype in sptype_prev[indx]:
                        print("%s%s (%.2f) is already registered as lineID=%d." % (
                            atom, ionstate_converter(ionstate), lam, prevIDs[indx]))
                    else:
                        cur.execute(
                            "UPDATE StellarLine SET sptype='%s' where lineID=%d;" % (
                                sptype_prev[indx] + sptype, prevIDs[indx])
                        )
                        print("Sptype is updated from '%s' to '%s' for lineID=%d." % (
                            sptype_prev[indx], sptype_prev[indx] + sptype, prevIDs[indx]))
        cur.execute(
            "INSERT IGNORE INTO StellarLine (lineID, wavelength_vac,wavelength_air,atom,ionstate,sptype) VALUES(%d,%.3f,%.3f,'%s',%d,'%s');" % (
                curID, lamvac, lamair, atom, ionstate, sptype))

    conn.commit()
    conn.close()


def GetLineList(lamstart, lamend, atom="all", ionstate="all", sptype="all", vacorair="AIR"):
    conn, cur = openproject()
    if vacorair == "AIR":
        cur.execute(
            "select lineID,wavelength_vac,wavelength_air,atom,ionstate,sptype from StellarLine where wavelength_air between %.3f and %.3f order by wavelength_air;" % (
                lamstart, lamend))
    elif vacorair == "VAC":
        cur.execute(
            "select lineID,wavelength_vac,wavelength_air,atom,ionstate,sptype from StellarLine where wavelength_vac between %.3f and %.3f order by wavelength_air;" % (
                lamstart, lamend))
    rows = cur.fetchall()

    conn.close()

    if rows == []:
        print("No stellar lines found between %.1f and %.1f." % (lamstart, lamend))
        return None
    else:
        lineIDs = numpy.array([int(i[0]) for i in rows])
        wav_vac = numpy.array([float(i[1]) for i in rows])
        wav_air = numpy.array([float(i[2]) for i in rows])
        atom_db = numpy.array([i[3] for i in rows])
        ionstate_db = numpy.array([int(i[4]) for i in rows])
        sptype_db = numpy.array([i[5] for i in rows])

        if atom == "all":
            atom_bool = numpy.array([True for i in atom_db])
        else:
            atom_bool = atom_db == atom

        if ionstate == "all":
            ionstate_bool = numpy.array([True for i in ionstate_db])
        else:
            ionstate_bool = ionstate_db == ionstate

        if sptype == "all":
            sptype_bool = numpy.array([True for i in ionstate_db])
        else:
            sptype_bool = numpy.array([sptype in sptype_db[i] for i in range(len(sptype_db))])

        all_bool = sptype_bool * atom_bool * ionstate_bool

        if not any(all_bool):
            print("No stellar lines satisfying atom=%s, ionstate=%s, and sptype=%s." % (atom, ionstate, sptype))
            return None

        if vacorair == "AIR":
            return [lineIDs[all_bool], wav_air[all_bool], atom_db[all_bool], ionstate_db[all_bool], sptype_db[
                all_bool]]
        elif vacorair == "VAC":
            return [lineIDs[all_bool], wav_vac[all_bool], atom_db[all_bool], ionstate_db[all_bool], sptype_db[
                all_bool]]


def AddDIBDB(lam, fwhm, reference, category, comment, vacorair="AIR"):
    conn, cur = openproject()
    cur.execute("select DIBID,wavelength_vac,wavelength_air,reference,category,FWHM from DIBlist;")
    rows = cur.fetchall()
    if rows == []:
        curID = 1
        newflag = True
    else:
        prevIDs = [int(i[0]) for i in rows]
        wav_vac_prev = [float(i[1]) for i in rows]
        wav_air_prev = [float(i[2]) for i in rows]
        reference_prev = [i[3] for i in rows]
        category_prev = [i[4] for i in rows]
        FWHM_prev = [float(i[5]) for i in rows]
        curID = max(prevIDs) + 1
        newflag = False

    if type(lam) is list:
        for i in range(len(lam)):
            if vacorair == "VAC":
                if not newflag:
                    if lam[i] in wav_vac_prev:
                        indx = wav_vac_prev.index(lam[i])
                        if reference[i] == reference_prev[indx] and category[i] == category_prev[indx]:
                            print("%.2f (%s, %s) is already registered as DIBID=%d." % (
                                lam[i], reference[i], category[i], prevIDs[indx]))
                            continue
                cur.execute(
                    "INSERT IGNORE INTO DIBlist (DIBID, wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment) VALUES(%d,%.3f,%.3f,'%s','%s',%.2f,%.3f,'%s');" % (
                        curID, lam[i], vac2air(lam[i]), reference[i], category[i], fwhm[i], wavenumber(lam[i]),
                        comment[i]))
                curID += 1
            if vacorair == "AIR":
                if not newflag:
                    if lam[i] in wav_air_prev:
                        indx = wav_air_prev.index(lam[i])
                        if reference[i] == reference_prev[indx] and category[i] == category_prev[indx]:
                            print("%.2f (%s, %s) is already registered as DIBID=%d." % (
                                lam[i], reference[i], category[i], prevIDs[indx]))
                            continue
                cur.execute(
                    "INSERT IGNORE INTO DIBlist (DIBID, wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment) VALUES(%d,%.3f,%.3f,'%s','%s',%.2f,%.3f,'%s');" % (
                        curID, air2vac(lam[i]), lam[i], reference[i], category[i], fwhm[i], wavenumber(air2vac(lam[i])),
                        comment[i]))
                curID += 1

    else:
        alreadyflag = False
        if vacorair == "VAC":
            if not newflag:
                if lam in wav_vac_prev:
                    indx = wav_vac_prev.index(lam)
                    if reference == reference_prev[indx] and category == category_prev[indx]:
                        print("%.2f (%s, %s) is already registered as DIBID=%d." % (
                            lam, reference, category, prevIDs[indx]))
                        alreadyflag = True
            if not alreadyflag:
                cur.execute(
                    "INSERT IGNORE INTO DIBlist (DIBID, wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment) VALUES(%d,%.3f,%.3f,'%s','%s',%.2f,%.3f,'%s');" % (
                    curID, lam, vac2air(lam), reference, category, fwhm, wavenumber(lam),
                    comment))
        if vacorair == "AIR":
            if not newflag:
                if lam in wav_air_prev:
                    indx = wav_air_prev.index(lam)
                    if reference == reference_prev[indx] and category == category_prev[indx]:
                        print("%.2f (%s, %s) is already registered as DIBID=%d." % (
                            lam, reference, category, prevIDs[indx]))
                        alreadyflag = True
            if not alreadyflag:
                cur.execute(
                    "INSERT IGNORE INTO DIBlist (DIBID, wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment) VALUES(%d,%.3f,%.3f,'%s','%s',%.2f,%.3f,'%s');" % (
                    curID, air2vac(lam), lam, reference, category, fwhm, wavenumber(air2vac(lam)),
                    comment))

    conn.commit()
    conn.close()


def UpdateDIBDB(DIBID, lam="INDEF", fwhm="INDEF", reference="INDEF", category="INDEF", comment="INDEF", vacorair="AIR", confirm=True):
    conn, cur = openproject()

    [DIBID_db, wav, reference_db, category_db, fwhm_db, wavenumber_db, comment_db] = GetDIBLine(DIBID, vacorair=vacorair)
    updateconfirm = "Updated DIB = %d\n\n" % DIBID

    if lam != "INDEF":
        if vacorair == "VAC":
            lamvac, lamair = lam, vac2air(lam)
        elif vacorair == "AIR":
            lamvac, lamair = air2vac(lam), lam
        wavenum = wavenumber(lamvac)
        cur.execute("UPDATE DIBlist set wavelength_vac = %.4f where DIBID=%d;" % (lamvac, DIBID))
        cur.execute("UPDATE DIBlist set wavelength_air = %.4f where DIBID=%d;" % (lamair, DIBID))
        cur.execute("UPDATE DIBlist set wavenumber = %.4f where DIBID=%d;" % (wavenum, DIBID))

        if vacorair == "VAC":
            updateconfirm += "wavelength_vac: %.4f --> %.4f\n" % (wav, lamvac)
        elif vacorair == "AIR":
            updateconfirm += "wavelength_air: %.4f --> %.4f\n" % (wav, lamair)

    if fwhm != "INDEF":
        cur.execute("UPDATE DIBlist set FWHM = %.2f where DIBID=%d;" % (fwhm, DIBID))
        updateconfirm += "FWHM: %.2f --> %.2f\n" % (fwhm_db, fwhm)

    if reference != "INDEF":
        cur.execute("UPDATE DIBlist set reference = '%s' where DIBID=%d;" % (reference, DIBID))
        updateconfirm += "Reference: %s --> %s\n" % (reference_db, reference)

    if category != "INDEF":
        cur.execute("UPDATE DIBlist set category = '%s' where DIBID=%d;" % (category, DIBID))
        updateconfirm += "Category: %s --> %s\n" % (category_db, category)

    if comment != "INDEF":
        cur.execute("UPDATE DIBlist set comment='%s' where DIBID=%d;" % (comment, DIBID))
        updateconfirm += "Comment: %s --> %s\n" % (comment_db, comment)

    if confirm:
        print(updateconfirm)
        ans = input("Update? ('y'es/no): ")
        if ans == "yes" or ans == "y":
            conn.commit()
            conn.close()
            return True
        else:
            conn.close()
            return False
    else:
        conn.commit()
        conn.close()
        return True


def GetDIBList(lamstart, lamend, reference="all", category="all", vacorair="AIR"):
    conn, cur = openproject()
    if vacorair == "AIR":
        cur.execute(
            "select DIBID,wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment from DIBlist where wavelength_air between %.3f and %.3f order by wavelength_air;" % (
                lamstart, lamend))
    elif vacorair == "VAC":
        cur.execute(
            "select DIBID,wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment from DIBlist where wavelength_vac between %.3f and %.3f order by wavelength_air;" % (
                lamstart, lamend))
    rows = cur.fetchall()

    conn.close()

    if rows == []:
        print("No lines found between %.1f and %.1f." % (lamstart, lamend))
        return None
    else:
        DIBIDs = numpy.array([int(i[0]) for i in rows])
        wav_vac = numpy.array([float(i[1]) for i in rows])
        wav_air = numpy.array([float(i[2]) for i in rows])
        reference_db = numpy.array([i[3] for i in rows])
        category_db = numpy.array([i[4] for i in rows])
        fwhm_db = numpy.array([float(i[5]) for i in rows])
        wavenumber_db = numpy.array([float(i[6]) for i in rows])
        comment_db = numpy.array([i[7] for i in rows])

        if reference == "all":
            reference_bool = numpy.array([True for i in reference_db])
        else:
            reference_bool = reference_db == reference

        if category == "all":
            category_bool = numpy.array([True for i in category_db])
        else:
            category_bool = category_db == category

        all_bool = reference_bool * category_bool

        if vacorair == "AIR":
            return [DIBIDs[all_bool], wav_air[all_bool], reference_db[all_bool], category_db[all_bool], fwhm_db[
                all_bool], wavenumber_db[all_bool], comment_db[all_bool]]
        elif vacorair == "VAC":
            return [DIBIDs[all_bool], wav_vac[all_bool], reference_db[all_bool], category_db[all_bool], fwhm_db[
                all_bool], wavenumber_db[all_bool], comment_db[all_bool]]


def GetDIBLine(dibid, vacorair="AIR"):
    conn, cur = openproject()
    cur.execute(
        "select DIBID,wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment from DIBlist where DIBID=%d;" % (
            dibid))
    rows = cur.fetchall()

    conn.close()

    if rows == []:
        print("No lines found for DIBID=%d." % (dibid))
        return None
    else:
        i = rows[0]
        DIBIDs = int(i[0])
        wav_vac = float(i[1])
        wav_air = float(i[2])
        reference_db = i[3]
        category_db = i[4]
        fwhm_db = float(i[5])
        wavenumber_db = float(i[6])
        comment_db = i[7]

    if vacorair == "AIR":
        return [DIBIDs, wav_air, reference_db, category_db, fwhm_db, wavenumber_db, comment_db]
    elif vacorair == "VAC":
        return [DIBIDs, wav_vac, reference_db, category_db, fwhm_db, wavenumber_db, comment_db]


def GetDIBdict(vacorair="AIR"):
    conn, cur = openproject()
    cur.execute(
        "select DIBID,wavelength_vac,wavelength_air,reference,category,FWHM,wavenumber,comment "
        "from DIBlist where category NOT IN ('Cs I', 'fake');")
    rows = cur.fetchall()

    conn.close()

    DIBIDlist = []
    wav_vac, wav_air, reference_db,category_db,fwhm_db,wavenumber_db,comment_db = {}, {}, {}, {}, {}, {}, {}

    if rows == []:
        print("No results, somehow...")
        return None
    else:
        for i in rows:
            d = int(i[0])
            DIBIDlist.append(d)
            wav_vac[d] = float(i[1])
            wav_air[d] = float(i[2])
            reference_db[d] = i[3]
            category_db[d] = i[4]
            fwhm_db[d] = float(i[5])
            wavenumber_db[d] = float(i[6])
            comment_db[d] = i[7]

    if vacorair == "AIR":
        return [DIBIDlist, wav_air, reference_db, category_db, fwhm_db, wavenumber_db, comment_db]
    elif vacorair == "VAC":
        return [DIBIDlist, wav_vac, reference_db, category_db, fwhm_db, wavenumber_db, comment_db]


def GetAFList(lambdamin, lambdamax, obsdate, echelleorder, vshift, slitposition="all", vacorair="AIR"):
    conn, cur = openproject()
    if type(obsdate) is not list:
        obsdate = [obsdate]
    if type(vshift) is not list:
        vshift = [vshift]

    lambda_vac_start_list = []
    lambda_vac_end_list = []
    slitpos = []
    AFID_list = []
    for i in range(len(obsdate)):
        if slitposition == "all":
            cur.execute(
                "select lambda_vac_start,lambda_vac_end,slitposition,AFID from Artificialfeatures where DATESTART < '%s' and DATEEND > '%s' and echelleorder=%d;" % (
                    obsdate[i], obsdate[i], echelleorder))
        else:
            cur.execute(
                "select lambda_vac_start,lambda_vac_end,slitposition,AFID from Artificialfeatures where DATESTART < '%s' and DATEEND > '%s' and slitposition='%s' and echelleorder=%d;" % (
                    obsdate[i], obsdate[i], slitposition, echelleorder))
        rows = cur.fetchall()
        if rows == []:
            return None
        else:
            lambda_vac_start_date = [float(i[0]) for i in rows]
            lambda_vac_end_date = [float(i[1]) for i in rows]
            slitpos_date = [i[2] for i in rows]
            AFID_date = [int(i[3]) for i in rows]
            for j in range(len(AFID_date)):
                # if not AFID_date[j] in AFID_list:
                lambda_vac_start_list.append(lambda_vac_start_date[j] * (1. + vshift[i] / (scipy.constants.c / 1.e+3)))
                lambda_vac_end_list.append(lambda_vac_end_date[j] * (1. + vshift[i] / (scipy.constants.c / 1.e+3)))
                slitpos.append(slitpos_date[j])
                AFID_list.append((AFID_date[j]))

    lambda_vac_start = numpy.array(lambda_vac_start_list)
    lambda_vac_end = numpy.array(lambda_vac_end_list)

    if vacorair == "VAC":
        lambda_start, lambda_end = lambda_vac_start, lambda_vac_end
    elif vacorair == "AIR":
        lambda_start, lambda_end = vac2air(lambda_vac_start), vac2air(lambda_vac_end)
    else:
        print("vacorair: %s" % vacorair)
        return None

    lambda_start_range, lambda_end_range = [], []
    for ls, le in zip(lambda_start, lambda_end):
        if lambdamin < le and lambdamax > ls:
            lambda_start_range.append(max(lambdamin, ls))
            lambda_end_range.append(min(lambdamax, le))

    return [lambda_start, lambda_end, slitpos]


if __name__ == "__main__":
    rf = open(sys.argv[1], "r")
    rl = rf.readlines()
    rf.close()

    if False:
        lamair, ionstate, atom = [], [], []
        for i in rl:
            comps = i.split("|")
            if len(comps) > 8:
                if comps[0].find("-----") == -1 and comps[0].find("lam") == -1 and comps[0].find("nm") == -1:
                    lamair.append(float(comps[0]))
                    atom.append(comps[2].split()[0])
                    ionstate.append(ionstate_converter_reverse(comps[2].split()[1]))
        sptype = ["A" for i in range(len(lamair))]
        reference = ["Sameshima et al. (2018), ApJS, 239, 19" for i in range(len(lamair))]
        AddLineDB(lamair, atom, ionstate, sptype, reference)

    if False:
        lamair, ew, fwhm, category = [], [], [], []
        for i in rl:
            rlcomp = i.split()
            if len(rlcomp) > 2:
                lamair.append(float(rlcomp[0]))
                ew.append(float(rlcomp[-1]))
                fwhm.append(float(rlcomp[1]))
                catetmp = "strong" if ew[-1] > 300. else "medium"
                category.append(catetmp)
        comment = ["Measured for HD183143" for i in lamair]
        reference = ["Friedman et al. (2011), ApJ, 727, 33" for i in lamair]
        AddDIBDB(lamair, fwhm, reference, category, comment)

    if False:
        lamair = [float(i.split()[1]) for i in rl]
        ew = [float(i.split()[3]) for i in rl]
        fwhm = [float(i.split()[2]) for i in rl]
        category = ["strong" if i > 300. else "medium" for i in ew]
        comment = ["Measured for Cyg OB2 No.12 data obtained in 2014-10-17." for i in rl]
        reference = ["New DIBs in SH et al. in prep." for i in rl]
        # refdict = {"Foing1994":"Foing et al. (1994), Nature, 369, 6478, 296.",
        #            "Groh2007":"Groh et al. (2007), A&A, 465, 993.",
        #            "Joblin1990":"Joblin et al. (1990), Nature, 346, 6286, 729"}
        # reference = [refdict[i.split()[5]] for i in rl]

        AddDIBDB(lamair, fwhm, reference, category, comment)

    if True:  # For stellar line
        atom, ionstate, lamair = [], [], []
        for i in rl:
            if i[0] != "#":
                atom.append(i[0:2].split()[0])
                ionstate.append(int(i[2]))
                lamair.append(float(i[4:].split()[0]))
        sptype = ["all" if a == "H" else "OB" for a in atom]
        reference = ["Sameshima et al. in prep., based on theoretical SPTOOL spectrum." for a in atom]

        AddLineDB(lamair, atom, ionstate, sptype, reference)

    if False:  # For HK-band DIBs (Galazutdinov et al. 2017)
        lamair = [float(i.split()[0]) for i in rl]
        ew = [float(i.split()[1]) for i in rl]
        fwhm = [float(i.split()[2]) for i in rl]
        category = ["strong" if i > 300. else "medium" for i in ew]
        comment = ["Measured for BD+40 4220." for i in rl]
        reference = ["Galazutdinov et al. (2017), MNRAS, 467, 3099." for i in rl]

        AddDIBDB(lamair, fwhm, reference, category, comment)

    if False:  # For H-band DIBs (Cox et al. 2014)
        lamair = [float(i.split()[0]) for i in rl]
        ew = [float(i.split()[1]) for i in rl]
        fwhm = [float(i.split()[2]) for i in rl]
        obj = [i.split()[3] for i in rl]
        category = ["strong" if i > 300. else "medium" for i in ew]
        comment = ["Measured for HD183143." if i == "HD183143" else "Measured for 4U 1907" for i in obj]
        reference = ["Cox et al. (2014), A&A, 569, A117." for i in rl]

        AddDIBDB(lamair, fwhm, reference, category, comment)
