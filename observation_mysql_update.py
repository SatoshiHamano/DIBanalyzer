#!/usr/bin/env python
# -*- coding: utf-8 -*-
## modules
import urllib.request
from urllib.parse import urlparse
import mysql.connector
from open_mysql_project import openproject

if __name__ == "__main__":

    url = 'http://merlot.kyoto-su.ac.jp/WINERED/WODB/obslog/obslog.php?NAME=&THEME=&FROM_YEAR=2015&FROM_MONTH=1&FROM_DAY=1&TO_YEAR=2030&TO_MONTH=12&TO_DAY=31&DTYPES%5B%5D=OBJECT&DTYPES%5B%5D=STANDARD&MODES%5B%5D=WIDE&MODES%5B%5D=HIRES-Y&MODES%5B%5D=HIRES-J&MODES%5B%5D=NULL&SLITS%5B%5D=100&SLITS%5B%5D=140&SLITS%5B%5D=200&SLITS%5B%5D=400&SLITS%5B%5D=0&LIMIT=100000'
    htmldata = urllib.request.urlopen(url)
    htmllines_enc = htmldata.readlines()
    htmllines = [i.decode() for i in htmllines_enc]

    conn, cur = openproject()

    cur.execute("select registeredname,objectid from objectdict;")
    rows = cur.fetchall()
    if rows == []:
        objlist = []
        objidlist = []
    else:
        objlist = [rows[i][0].upper() for i in range(len(rows))]
        objidlist = [rows[i][1] for i in range(len(rows))]
        # for i in range(len(objlist)):
        #     print(objidlist[i],objlist[i])

    cur.execute("select frame,objectid,object from observation;")
    rows = cur.fetchall()
    if rows == []:
        framelist = []
        objidobslist = []
    else:
        framelist = [rows[i][0] for i in range(len(rows))]
        objidobslist = [rows[i][1] for i in range(len(rows))]
        objname = [rows[i][2] for i in range(len(rows))]
        # for i in range(len(framelist)):
        #     print(framelist[i],objidobslist[i],objname[i])

    for i in range(len(htmllines)):
        if htmllines[i].find("<tr class=") != -1:
            rlsplit = [htmllines[i].split("<td>")[j].rstrip("</tr>\n").rstrip("td>").rstrip("</") for j in
                       range(1, len(htmllines[i].split("<td>")))]
            [frameid, target, theme, type, acqtime, mode, slit, exptime, position, airmass, quality, memo] = rlsplit
            if airmass == "": airmass = "10"

            if target.upper() in objlist:
                objid = objidlist[objlist.index(target.upper())]
                if not frameid in framelist:
                    cur.execute(
                        "INSERT IGNORE INTO observation (frame, objectID, object, theme, type, obsdate, instmode, slit, exptime, position, airmass, flag, memo) VALUES('%s',%s,'%s','%s','%s',cast('%s' as datetime),'%s',%s,%s,'%s',%s,'%s','%s');" % (
                        frameid, objid, target, theme, type, acqtime, mode, slit, exptime, position, airmass, quality,
                        memo))
                    print("Frame %s is inserted Object ID %s (Theme: %s)." % (frameid, objid, theme))
                else:
                    if objidobslist[framelist.index(frameid)] != objid:
                        cur.execute("UPDATE observation SET objectid = %s where frame = '%s';" % (objid, frameid))
                        print("ObjectID is updated (from %s to %s) for frame %s (%s)." % (
                        objidobslist[framelist.index(frameid)], objid, frameid, target))
            elif theme.find("DIB") != -1 or type == "standard" or theme.find("C2") != -1:
                print("%s (%s) is not registered in object table." % (target, theme))


    conn.commit()
    conn.close()
