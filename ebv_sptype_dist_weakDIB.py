# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

if __name__ == '__main__':
    conn, cur = openproject()

    cur.execute(
        "select od.registeredname,z.E_BV,z.sptype,dr.totalSNR "
        "from object as z "
        "join objectdict as od using (objectID) "
        "join datareduction as dr using (objectID) "
        "where priority=1 and objectID IN (138,50,52,34,60,23,29,22,35,24,44,65,37,28,33,43,30,123,31,25,42,26,124,90,246,32,27,63,45,128,8,59,41,66,40,39,14,10,12,13,15,135,9,132,144);")

    rows = cur.fetchall()
    objname = np.array([i[0] for i in rows])
    ebv = np.array([i[1] for i in rows])
    sptype = np.array([i[2].split()[0] for i in rows])
    totalSNR = np.array([i[3] for i in rows])

    print(len([138,50,52,34,60,23,29,22,35,24,44,65,37,28,33,43,30,123,31,25,42,26,124,90,246,32,27,63,45,128,8,59,41,66,40,39,14,10,12,13,15,135,9,132,144]), " objects")

    tempscore = {"O": 0, "B": 10, "A":20, "F":30, "G":40}

    sptypescore = []
    for i in sptype:
        try:
            sptypescore.append(tempscore[i[0]] + int(i[1]))
        except:
            print("ERROR:", i)
            sys.exit()


    plt.figure()
    # sc = plt.scatter(ebvlist_all, sptypescore_all, s=20, vmin=200, vmax=800, c=totalSNRlist_all, cmap=cm.winter_r, edgecolors=None, marker="^", label="Not analyzed")
    sc = plt.scatter(ebv, sptypescore, s=35, c=totalSNR, vmin=200, vmax=800, cmap=cm.winter_r, edgecolors="k", linewidths=1)
    cbar = plt.colorbar(sc)
    cbar.set_label("Signal-to-noise ratio")
    # plt.legend()
    plt.xlabel("E(B-V)")
    plt.ylabel("Sptype")
    plt.yticks([5, 10, 15, 20], ["O5", "B0", "B5", "A0"])
    plt.ylim(1,25)
    # plt.xlim(-0.1, 2.0)
    plt.savefig("temporaly_files/ebv_sptype_snr_dist_weakDIB.png", format="png", dpi=200)
    plt.clf()

    conn.close()