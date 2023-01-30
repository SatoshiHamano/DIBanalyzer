#!/usr/bin/env python
# -*- coding: utf-8 -*-
## modules
import os, time, glob, sys
import mysql.connector
from urllib.parse import urlparse
from open_mysql_project import openproject


if __name__ == "__main__":
    conn, cur = openproject()

    rf = open(sys.argv[1], "r")
    rl = rf.readlines()
    rf.close()
    pipid = [i.split()[0] for i in rl]

    wf = open(sys.argv[2], "w")
    wf.write("--- Data reduction info of each dataset. ---\n\n")
    wf.write(
        "pipelineID, FrameNum, totalexp, totalSNR, mode, pipelinever, scatteredlight, manual, background, hotpix\n")

    for i in range(len(pipid)):
        cur.execute("select pipelineID, FrameNum, totalexp, totalSNR, mode, pipelinever, scatteredlight, manual, background, hotpix from datareduction where pipelineID='%s';" % pipid[i])

        rows = cur.fetchall()
        if len(rows) == 1:
            [pid, frn, exp, snr, mode, ver, sclight, manual, bg, hp] = rows[0]
        else:
            print(len(rows), "records are hit! No info is shown here.\n\n")
            break

        wf.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n" % (pid, frn, exp, snr, mode, ver, sclight, manual, bg, hp))

    wf.write("\n\n--- Data reduction info of each frame. ---\n\n")
    wf.write("pipelineID, objectframe, skyframe, exptime, aperture, background, badpix, wsave, wsstd, wsnum\n")

    for i in range(len(pipid)):
        cur.execute("select objectframe, skyframe, exptime, aperture, background, badpix, wsave, wsstd, wsnum from reducedframe where pipelineID='%s';" % pipid[i])
        rows = cur.fetchall()
        obj = [j[0] for j in rows]
        sky = [j[1] for j in rows]
        expf = [j[2] for j in rows]
        ap = [j[3] for j in rows]
        bgf = [j[4] for j in rows]
        bp = [j[5] for j in rows]
        wsave = [j[6] for j in rows]
        wsstd = [j[7] for j in rows]
        wsnum = [j[8] for j in rows]
        wf.write(pipid[i] + "\n")
        for j in range(len(obj)):
            wf.write("%s, %s, %s, %s, %s, %s, %s, %s, %s\n" % (obj[j], sky[j], expf[j], ap[j], bgf[j], bp[j], wsave[j], wsstd[j], wsnum[j]))
        wf.write("\n")

    conn.close()
    wf.close()
