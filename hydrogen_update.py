#!/usr/bin/env python
# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import sys
import numpy as np

if __name__ == '__main__':
    conn, cur = openproject()

    rf = open(sys.argv[1], "r")
    rl = rf.readlines()
    rf.close()

    oid = np.array([int(i.split("\t")[0]) for i in rl])
    ebv = np.array([float(i.split("\t")[1]) for i in rl])
    logNHI = np.array([float(i.split("\t")[2]) for i in rl])
    logNHIerr = np.array([float(i.split("\t")[3]) for i in rl])
    logNH2 = np.array([float(i.split("\t")[4]) for i in rl])
    logNH2err = np.array([float(i.split("\t")[5]) for i in rl])
    comment = np.array([i.split("\t")[6] for i in rl])

    reference = input("Reference: ")

    for i in range(len(oid)):
        print("object ID = {}".format(oid[i]))
        cur.execute("SELECT log_N_HI, log_N_HIerr, log_N_HIupp, log_N_H2, log_N_H2err, log_N_H2upp, N_HI_reference, N_H2_reference "
                    "from object where objectid={};".format(oid[i]))
        rows = cur.fetchall()
        logNHI_old = rows[0][0]
        logNHIerr_old = rows[0][1]
        logNHIupp_old = rows[0][2]
        logNH2_old = rows[0][3]
        logNH2err_old = rows[0][4]
        logNH2upp_old = rows[0][5]
        N_HI_reference_old = rows[0][6]
        N_H2_reference_old = rows[0][7]

        if logNHI[i] != 0. and logNHIerr[i] != 0.:
            cur.execute("UPDATE object set log_N_HI={} where objectID={};".format(logNHI[i], oid[i]))
            cur.execute("UPDATE object set log_N_HIerr={} where objectID={};".format(logNHIerr[i], oid[i]))
            cur.execute("UPDATE object set log_N_HIupp=Null where objectID={};".format(oid[i]))
            cur.execute("UPDATE object set N_HI_reference={} where objectID={};".format(reference, oid[i]))
            if logNHI_old != None or logNHIerr_old != None or N_HI_reference_old != None:
                print("log_N_HI is updated from {} to {}.".format(logNHI_old, logNHI[i]))
                print("log_N_HIerr is updated from {} to {}.".format(logNHIerr_old, logNHIerr[i]))
                print("log_N_HIupp is updated from {} to {}.".format(logNHIupp_old, "Null"))
                print("N_HI_reference is updated from {} to {}.".format(N_HI_reference_old, reference))
        elif logNHI[i] != 0 and logNHIerr[i] == 0:
            cur.execute("UPDATE object set log_N_HI=Null where objectID={};".format(oid[i]))
            cur.execute("UPDATE object set log_N_HIerr=Null where objectID={};".format(oid[i]))
            cur.execute("UPDATE object set log_N_HIupp={} where objectID={};".format(logNHI[i], oid[i]))
            cur.execute("UPDATE object set N_HI_reference={} where objectID={};".format(reference, oid[i]))
            if logNHI_old != None or logNHIerr_old != None or N_HI_reference_old != None:
                print("log_N_HI is updated from {} to {}.".format(logNHI_old, "Null"))
                print("log_N_HIerr is updated from {} to {}.".format(logNHIerr_old, "Null"))
                print("log_N_HIupp is updated from {} to {}.".format(logNHIupp_old, logNHI[i]))
                print("N_HI_reference is updated from {} to {}.".format(N_HI_reference_old, reference))

        if logNH2[i] != 0. and logNH2err[i] != 0.:
            cur.execute("UPDATE object set log_N_H2={} where objectID={};".format(logNH2[i], oid[i]))
            cur.execute("UPDATE object set log_N_H2err={} where objectID={};".format(logNH2err[i], oid[i]))
            cur.execute("UPDATE object set log_N_H2upp=Null where objectID={};".format(oid[i]))
            cur.execute("UPDATE object set N_H2_reference={} where objectID={};".format(reference, oid[i]))
            if logNH2_old != None or logNH2err_old != None or N_H2_reference_old != None:
                print("log_N_H2 is updated from {} to {}.".format(logNH2_old, logNH2[i]))
                print("log_N_H2err is updated from {} to {}.".format(logNH2err_old, logNH2err[i]))
                print("log_N_H2upp is updated from {} to {}.".format(logNH2upp_old, "Null"))
                print("N_H2_reference is updated from {} to {}.".format(N_H2_reference_old, reference))
        elif logNH2[i] != 0 and logNH2err[i] == 0:
            cur.execute("UPDATE object set log_N_H2=Null where objectID={};".format(oid[i]))
            cur.execute("UPDATE object set log_N_H2err=Null where objectID={};".format(oid[i]))
            cur.execute("UPDATE object set log_N_H2upp={} where objectID={};".format(logNH2[i], oid[i]))
            cur.execute("UPDATE object set N_H2_reference={} where objectID={};".format(reference, oid[i]))
            if logNH2_old != None or logNH2err_old != None or N_H2_reference_old != None:
                print("log_N_H2 is updated from {} to {}.".format(logNH2_old, "Null"))
                print("log_N_H2err is updated from {} to {}.".format(logNH2err_old, "Null"))
                print("log_N_H2upp is updated from {} to {}.".format(logNH2upp_old, logNH2[i]))
                print("N_H2_reference is updated from {} to {}.".format(N_H2_reference_old, reference))

        conn.commit()
    conn.close()