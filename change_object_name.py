# -*- coding: utf-8 -*-

import argparse
import os
import sys
from open_mysql_project import openproject
from telluric_auto import alternativequestion

if __name__ == '__main__':
    parantdir = "/Users/hamano/DIB_analysis/DIB_pipeline_dir/"

    parser = argparse.ArgumentParser()
    parser.add_argument("objectid", type=int, help="objectid")
    parser.add_argument("newname", type=str, help="new short name")

    args = parser.parse_args()

    objectid = args.objectid
    newname = args.newname

    conn, cur = openproject()

    cur.execute("SELECT objectname from object where objectid=%d;" % objectid)
    rows = cur.fetchall()
    if rows == []:
        print("No object is registered with ID=%d." % objectid)
        sys.exit()
    else:
        oldname = rows[0][0]
        oldname_formatted = oldname.replace("*", "").lstrip(" ").replace("   ", " ").replace("  ", " ").replace(" ",
                                                                                                                "_").replace(
            "[", "").replace("]", "")

    ans = alternativequestion("Change name of ID=%d from '%s' to '%s'?:  " % (objectid, oldname, newname),
                              ["yes", "no"], "yes")
    if ans == "no":
        print("bye/~~~")
        sys.exit()

    print("## Change database records ##")

    if os.path.exists(parantdir + oldname_formatted):
        os.rename(parantdir + oldname_formatted, parantdir + newname)
        print("Directory was renamed to %s." % (parantdir + newname))

    cur.execute("UPDATE object SET objectname='%s' where objectid = %d;" % (newname, objectid))
    conn.commit()
    print("object.objectname updated.")
    print(newname)

    cur.execute("UPDATE objectdict SET registeredname='%s' where objectid = %d and registeredname='%s';" % (
    newname, objectid, oldname))
    conn.commit()
    print("objectdict.registeredname updated.")
    print(newname)

    cur.execute(
        "SELECT x.telluricID from telluriccorrection as x join datareduction as y on x.pipelineIDobj=y.pipelineID where y.objectID=%d;" % objectid)
    rows = cur.fetchall()
    telluricID = [i[0] for i in rows]
    for t in telluricID:
        cur.execute(
            "SELECT telluricpath from telluriccorrection where telluricID='%s';" % t)
        rows = cur.fetchall()
        telluricpath = rows[0][0]
        cur.execute("UPDATE telluriccorrection SET telluricpath='%s' where telluricID='%s';" % (
        telluricpath.replace(oldname_formatted, newname), t))
        conn.commit()
        print("telluriccorrection.telluricpath updated")
        print(telluricpath.replace(oldname_formatted, newname))

        cur.execute(
            "SELECT telluricfilepath from telluricresult where telluricID='%s';" % t)
        rows = cur.fetchall()
        telluricfilepath = [i[0] for i in rows]
        for f in telluricfilepath:
            cur.execute(
                "UPDATE telluricresult SET telluricfilepath='%s' where telluricID='%s' and telluricfilepath='%s';" % (
                    f.replace(oldname_formatted, newname), t, f))
            conn.commit()
            print("telluriccorrection.telluricfilepath updated")
            print(f.replace(oldname_formatted, newname))

    cur.execute("SELECT combineID, combinepath from combinesummary where objectid=%d;" % objectid)
    rows = cur.fetchall()
    if rows == []:
        print("No combine dataset.")
    else:
        combineID = [i[0] for i in rows]
        combinepath = [i[1] for i in rows]

        for c in range(len(combineID)):
            cur.execute(
                "UPDATE combinesummary SET combinepath='%s' where combineID='%s';" % (
                    combinepath[c].replace(oldname_formatted, newname), combineID[c]))
            conn.commit()
            print("combinesummary.combinepath updated")
            print(combinepath[c].replace(oldname_formatted, newname))
            combdirold = combinepath[c].split("/")[-2]
            combdirnew = combinepath[c].replace(oldname_formatted, newname).split("/")[-2]
            combparent = combinepath[c].split(combdirold)[0].replace(oldname_formatted, newname)
            if combdirold != combdirnew:
                os.rename(combparent + combdirold, combparent + combdirnew)
                print("Directory was renamed to %s." % (combparent + combdirnew))

            cur.execute(
                "UPDATE combinesummary SET combineID='%s' where combineID='%s';" % (
                    combineID[c].replace(oldname_formatted, newname), combineID[c]))
            conn.commit()
            print("combinesummary.combineID updated")
            print(combineID[c].replace(oldname_formatted, newname))

            cur.execute(
                "UPDATE combinedataset SET combineID='%s' where combineID='%s';" % (
                    combineID[c].replace(oldname_formatted, newname), combineID[c]))
            conn.commit()
            print("combinedataset.combineID updated")
            print(combineID[c].replace(oldname_formatted, newname))

            cur.execute(
                "UPDATE combinedspectrum SET combineID='%s' where combineID='%s';" % (
                    combineID[c].replace(oldname_formatted, newname), combineID[c]))
            conn.commit()
            print("combinedspectrum.combineID updated")
            print(combineID[c].replace(oldname_formatted, newname))

            cur.execute("SELECT combinefilepath from combinedspectrum where combineID='%s';" % combineID[c])
            rows = cur.fetchall()
            combinefilepath = [i[0] for i in rows]
            for f in combinefilepath:
                cur.execute(
                    "UPDATE combinedspectrum SET combinefilepath='%s' where combineID='%s' and combinefilepath='%s';" % (
                        f.replace(oldname_formatted, newname), combineID[c], f))
                conn.commit()
                print("combinedspectrum.combinefilepath updated")
                print(f.replace(oldname_formatted, newname))

                combfilenameold = os.path.basename(f)
                combfilenamenew = os.path.basename(f.replace(oldname_formatted, newname))
                combfiledir = os.path.dirname(f.replace(oldname_formatted, newname)) + "/"
                if combfilenameold != combfilenamenew:
                    os.rename(combfiledir + combfilenameold, combfiledir + combfilenamenew)
                    print("File was renamedto %s" % combfilenamenew)

    conn.close()
