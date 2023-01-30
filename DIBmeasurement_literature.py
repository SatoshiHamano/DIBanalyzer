#!/usr/bin/env python
# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import mysql.connector
from urllib.parse import urlparse
import sys

if __name__ == "__main__":
    rf = open(sys.argv[1], "r")
    rl = rf.readlines()
    rf.close()

    boolstr = {"TRUE": True, "FALSE": False}
    conn, cur = openproject()

    reso = {'F':48000, 'E':68000, 'H':115000, 'B':45000}
    spectro = {'F':'The 2.2m ESO LaSilla telescope/FEROS', 'E': 'The 3.58m CFHT/ESPaDOnS', 'H':'The 3.6m ESO LaSilla telescope/HARPS', 'B':'The 1.8m telescope of the Bohyunsan Observatory/BOES'}


    for i in rl:
        a = i.rstrip("\n").split("\t")
        if len(a) == 9:
            [reference, dibid, objid, ew, ewerr, upperlimitflag, obsdate, resolution, telesspec] = a

            # cur.execute(
            #     "INSERT IGNORE INTO DIBEWsummary (DIBID, objectID, primaryflag, reference, EW, EWerr, obsdate, resolution, teles_spectrograph, upperlimit) VALUES (%d, %d, %s, '%s', %.3f, %.3f, '%s', %.0f, '%s', %s);" % (
            #         int(dibid), int(objid), "True", reference, float(ew), float(ewerr), obsdate, float(resolution),
            #         telesspec, boolstr[upperlimitflag]))
            cur.execute(
                "INSERT IGNORE INTO DIBEWsummary (DIBID, objectID, primaryflag, reference, EW, EWerr, obsdate, resolution, teles_spectrograph, upperlimit) VALUES (%d, %d, %s, '%s', %.3f, %.3f, '%s', %.0f, '%s', %s);" % (
                    int(dibid), int(objid), "True", reference, float(ew), float(ewerr), obsdate, reso[telesspec],
                    spectro[telesspec], boolstr[upperlimitflag]))
        elif len(a) == 10:
            [reference, dibid, objid, ew, ewerr, upperlimitflag, obsdate, resolution, telesspec, comment] = a
            cur.execute(
                "INSERT IGNORE INTO DIBEWsummary (DIBID, objectID, primaryflag, reference, EW, EWerr, obsdate, resolution, teles_spectrograph, upperlimit, comment) VALUES (%d, %d, %s, '%s', %.3f, %.3f, '%s', %.0f, '%s', %s, '%s');" % (
                    int(dibid), int(objid), "True", reference, float(ew), float(ewerr), obsdate, reso[telesspec],
                    spectro[telesspec], boolstr[upperlimitflag], comment))
            # cur.execute(
            #     "INSERT IGNORE INTO DIBEWsummary (DIBID, objectID, primaryflag, reference, EW, EWerr, obsdate, resolution, teles_spectrograph, upperlimit, comment) VALUES (%d, %d, %s, '%s', %.3f, %.3f, '%s', %.0f, '%s', %s, '%s');" % (
            #         int(dibid), int(objid), "True", reference, float(ew), float(ewerr), obsdate, float(resolution),
            #         telesspec, boolstr[upperlimitflag], comment))
        else:
            print(a)

    conn.commit()
    conn.close()
