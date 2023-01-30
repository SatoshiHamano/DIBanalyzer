#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import numpy

from open_mysql_project import openproject
from Spec1Dtools import FSR_angstrom

def read_centersearch_log(centersearch_log, apnum, objnum):
    rf = open(centersearch_log, "r")
    rl = rf.readlines()
    rf.close()

    readlist = []

    counter = 0
    for i in range(len(rl)):
        if rl[i].find("Files") != -1:
            readlist.append([[] for n in range(objnum)])
            for j in range(i + 2, i + apnum + 2):
                rl1comp = rl[j].split()
                for k in range(objnum):
                    readlist[counter][k].append(float(rl1comp[k + 1]))
            counter += 1

    [xs, gw] = readlist

    return xs, gw

if __name__ == '__main__':
    conn, cur = openproject()
    cur.execute("select pipelineID, path, frameNum from datareduction where totalSNR > 200 and mode = 'WIDE' and obsdate < '2017-01-01' and obsdate > '2015-01-01';")
    rows = cur.fetchall()
    ppid = [i[0] for i in rows]
    path = [i[1] for i in rows]
    n = [i[2] for i in rows]

    fsr = FSR_angstrom()
    m = [42 + i for i in range(20)]
    lam = numpy.array([numpy.average(fsr[i]) for i in m])

    wf = open("WINERED_ad_Araki.dat", "w")

    for i in range(len(ppid)):
        if os.path.exists(path[i] + "reduction_log/centersearch_log.txt"):
            try:
                xs, gw = read_centersearch_log(path[i] + "reduction_log/centersearch_log.txt", 20, n[i])
                cur.execute(
                    "select airmass from reducedframe join observation on reducedframe.objectframe = observation.frame  where pipelineID = '{}';".format(ppid[i]))
                rows = cur.fetchall()
                airmass = [a[0] for a in rows]

                for j in range(n[i]):
                    a, b = numpy.polyfit(lam, xs[j], 1)
                    wf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(ppid[i], j, airmass[j], a, b, numpy.std(xs[j]-lam*a-b)))
            except:
                pass

    wf.close()
    conn.close()