#!/usr/bin/env python
# -*- coding: utf-8 -*-

# +
# NAME:
#   make_target_list_run05.py
#
# PURPOSE:
#   Get information about the inputted object list 
#   by asking SIMBAD
#
# CALLING SEQUENCE:
#   python make_target_list_run05.py list.txt
#
# OUTPUT FILE:
#   target_list_run05.txt
#
# REQUIRED FILE:
#   list.txt (Any name is OK)
#
# REMARKS:
#   N/A
#
# AUTHOR:
#   Hiroaki Sameshima (sameshima@kyoto-su.ac.jp)
#
# REVISION HISTORY:
#   2015-04-15 Created
#   2015-04-26 Items to be output are changed.
# -

## modules
import os, subprocess, sys
import urllib.request
from optparse import OptionParser
from urllib.parse import urlparse
import mysql.connector

bandname = ["U", "B", "V", "G", "R", "I", "J", "H", "K"]
bandall = ""
for i in bandname:
    bandall += i


def alternativequestion(question, anss, defans):
    flagans = 0
    while flagans == 0:
        flag = input(question)
        if flag in anss:
            flagans += 1
        else:
            print("Answers: ", anss)

    if flag != "":
        return flag
    else:
        return defans


##############
# sub routine
##############
def simbad_query(objNameSearch, summary=0):

    ## adjust object name
    objName = objNameSearch
    objNameURL = objNameSearch.replace("+", "%2B").replace(" ", "+")

    ## get information from SIMBAD
    print('Now searching SIMBAD about "' + objName + '"...')

    url = 'http://simbad.u-strasbg.fr/simbad/sim-id?Ident=' + objNameURL
    htmldata = urllib.request.urlopen(url)
    htmllines_enc = htmldata.readlines()
    htmllines = [i.decode() for i in htmllines_enc]

    ## define flags
    flag_result = 0
    flag_specType = 0
    flag_propMot = 0
    flag_mag = 0
    flag_radvel = 0

    magnitudes = {}
    for i in bandname:
        magnitudes[i] = "Null"

    ## html analysis
    for i in range(len(htmllines)):
        ## SIMBAD name & Object Type
        if 'Basic data :' in htmllines[i]:
            flag_result = 1
            simbadName = htmllines[i + 13].rstrip()
            objType = htmllines[i + 16].rstrip()

        ## Coodinate
        if 'ICRS' in htmllines[i] and 'coord.' in htmllines[i + 2]:
            epoch = htmllines[i + 4].split()[0].replace("(ep=", "").rstrip(")")
            radec = htmllines[i + 11].split()
            ra = radec[0:3]
            dec = radec[3:6]

        ## Magnitude
        if 'Fluxes' in htmllines[i]:
            j = i + 9

            while not '</TABLE>' in htmllines[j]:
                if htmllines[j].split()[0] in bandname and len(htmllines[j].split()) > 1:
                    magnitudes[htmllines[j].split()[0]] = htmllines[j].split()[1]
                    flag_mag = 1
                j += 1

        ## Proper Motion
        if 'Proper motions' in htmllines[i]:
            flag_propMot = 1
            propMotions = htmllines[i + 5].split()  # mas/year
            propMotRA = float(propMotions[0]) / 1000.  # arcsec/year
            propMotDec = float(propMotions[1]) / 1000.  # arcsec/year

        ## Radial Velocity
        if 'Radial velocity' in htmllines[i]:
            flag_radvel = 1
            radVel = float(htmllines[i+5].split()[1])

        ## Spectral Type
        if 'Spectral type:' in htmllines[i]:
            flag_specType = 1
            tmp = htmllines[i + 5].split()
            specType = str(tmp[0])
            if len(tmp) > 1:
                for j in range(len(tmp) - 1):
                    specType += " " + str(tmp[j + 1])

    ## flag processing
    if flag_result == 0:
        print('\033[93m' + "!!! Warning !!!" + '\033[0m' + "  No SIMBAD entry could be found for : " + objName)
        return [objName, "Null", "Null", "Null", "Null", "Null", "Null", "Null", magnitudes, "Null", "Null"], flag_result
        quit()

    if flag_specType == 0:
        specType = "Null"

    if flag_propMot == 0:
        propMotRA = "0"
        propMotDec = "0"

    if flag_radvel == 0:
        radVel = "0"

    ## print output
    if summary == 1:
        print("Search Name    : " + objName)
        print("SIMBAD Name    : " + simbadName)
        print("Object Type    : " + objType)
        print("Spectral Type  : " + specType)
        print("FK5 Coord.     : " + str(ra[0]) + " " + str(ra[1]) + " " + str(ra[2]) + " " + str(dec[0]) + " " + str(
            dec[1]) + " " + str(dec[2]))
        print("Proper Motions : " + str(propMotRA) + " " + str(propMotDec))
        print("Radial velocity : " + str(radVel))
        print("Epoch          : " + str(epoch))
        print("J mag          : " + magnitudes["J"])
        print("R mag          : " + magnitudes["R"])
        print("SIMBAD URL     : " + url)

    return [objName, simbadName, objType, specType, ra, dec, [propMotRA, propMotDec], epoch, magnitudes,
            url, radVel], flag_result


if __name__ == '__main__':

    ##############
    # main
    ##############

    ## options
    usage = "usage: %prog [options] object_list"
    parser = OptionParser(description='Get information about the inputted object list by asking SIMBAD.')
    (options, objList) = parser.parse_args()

    ## error processing
    if len(objList) != 1:
        print('\033[93m' + "!!! Warning !!!" + '\033[0m' + "  Please specify only one list.")
        quit()

    ## read target list
    f = open(objList[0])
    lines = f.readlines()
    f.close()

    objstd = alternativequestion("OBJECT or STANDARD (default: OBJECT): ", ["OBJECT", "object", "STANDARD", "standard", ""],
                                 "OBJECT")

    if objstd == "OBJECT" or objstd == "object":
        obstype = "OBJECT"
    else:
        obstype = "STANDARD"

    ## prepare output files
    url = urlparse('mysql://root:kwmjbqb9py@localhost:3306/DIBproject')

    conn = mysql.connector.connect(
        host=url.hostname or 'localhost',
        port=url.port or 3306,
        user=url.username or 'root',
        password=url.password or 'kwmjbqb9py',
        database=url.path[1:],
    )

    cur = conn.cursor()

    cur.execute("select objectname from object;")
    rows = cur.fetchall()
    if rows == []:
        objlist = []
    else:
        objlist = [rows[i][0] for i in range(len(rows))]

    cur.execute("select registeredname from objectdict;")
    rows = cur.fetchall()
    if rows == []:
        objdictlist = []
    else:
        objdictlist = [rows[i][0] for i in range(len(rows))]

    prevflag = False

    ## process each line
    for line in lines:

        ## get object name,S/N,PI,theme and comment
        inputs = line.rstrip().split("#")
        if len(inputs) > 1:
            if inputs[1] == "same":
                flagsame = True
            else:
                flagsame = False
        else:
            flagsame = False

        name = inputs[0].rstrip()

        if not name in objdictlist:
            if flagsame:
                cur.execute("select objectid from objectdict where registeredname='%s';" % prevname)
                rows = cur.fetchall()
                objid = rows[0][0]

                cur.execute("INSERT IGNORE INTO objectdict (registeredname,objectid,priority) VALUES ('%s','%s',0);" % (name, objid))
                conn.commit()

            else:
                ## check if the selected line is NOT a comment
                if inputs[0][0] != "#":

                    ## query SIMBAD
                    result, flagsimbad = simbad_query(name)

                    ## prepare output
                    if flagsimbad != 0:
                        prevflag = True

                        name = result[0]
                        simbname = result[1]

                        if result[4] == "Null":
                            ra = "Null"
                        else:
                            ra = result[4][0] + ":" + result[4][1] + ":" + result[4][2][0:5]

                        if result[5] == "Null":
                            dec = "Null"
                        else:
                            dec = result[5][0] + ":" + result[5][1] + ":" + result[5][2][0:4]

                        if result[6] == "Null":
                            propMotRA = "Null"
                            propMotDec = "Null"
                        else:
                            propMotRA = str(result[6][0])  # arcsec/year
                            propMotDec = str(result[6][1])  # arcsec/year

                        if result[7] == "Null":
                            epoch = "Null"
                        else:
                            epoch = result[7].replace("J", "") + ".0"

                        mag_simbad = result[8]
                        magstr = ""
                        for i in bandname:
                            if i != bandname[-1]:
                                magstr += mag_simbad[i] + ","
                            else:
                                magstr += mag_simbad[i]

                        objType = result[2]
                        specType = result[3]

                        if specType == "Null":
                            note = objType
                        else:
                            note = objType + " (" + specType + ")"

                        URL = result[9]
                        if URL == "Null":
                            linkColumn = "Null"
                        else:
                            linkColumn = "[link|" + URL + "]"

                        if not simbname in objlist:
                            cur.execute(
                                "INSERT IGNORE INTO object (objectname,type,objtype,sptype,ra,decli,pmra,pmdec,epoch,Umag,Bmag,Vmag,Gmag,Rmag,Imag,Jmag,Hmag,Kmag) VALUES('%s','%s','%s','%s','%s','%s',%s,%s,%s,%s);" % (
                                    simbname, obstype, objType, specType, ra, dec, propMotRA, propMotDec, epoch, magstr))
                            conn.commit()
                            objlist.append(simbname)
                            priority = 1
                        else:
                            priority = 0

                        cur.execute("select objectid from object where objectname='%s';" % simbname)
                        rows = cur.fetchall()
                        objid = rows[0][0]
                        cur.execute("INSERT IGNORE INTO objectdict (registeredname,objectid,priority) VALUES ('%s','%s',%d);" % (name, objid, priority))
                        conn.commit()

                    else:
                        print("Not found in Simbad.")
                        cur.execute(
                            "INSERT IGNORE INTO object (objectname,type) VALUES('%s','%s');" % (
                                name, obstype))
                        conn.commit()
                        objlist.append(name)

                        cur.execute("select objectid from object where objectname='%s';" % name)
                        rows = cur.fetchall()
                        objid = rows[0][0]
                        cur.execute("INSERT IGNORE INTO objectdict (registeredname,objectid) VALUES ('%s','%s');" % (name, objid))
                        conn.commit()

        else:
            print("%s was registered." % name)

        prevname = name

    ## print remarks
    print('Finished! "target_list_run05.txt" is created.')
    conn.close()
