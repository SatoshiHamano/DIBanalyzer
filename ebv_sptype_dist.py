# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if __name__ == '__main__':
    conn, cur = openproject()

    cur.execute(
        "select od.registeredname,z.E_BV,z.sptype,x.combineID,x.EW,cd.datasetID,cd.weight,tc.advanced,dr.totalSNR "
        "from DIBmeasurement as x "
        "join combinesummary as y using(combineID) "
        "join object as z using (objectID) "
        "join objectdict as od using (objectID) "
        "join combinedataset as cd using (combineID) "
        "join telluriccorrection as tc on cd.datasetID=tc.telluricID "
        "join datareduction as dr on tc.pipelineIDobj=dr.pipelineID "
        "where DIBID=38 and primaryflag =1 and od.priority=1 order by datasetID;")

    rows = cur.fetchall()
    objname = np.array([i[0] for i in rows])
    ebv = np.array([i[1] for i in rows])
    sptype = np.array([i[2].split()[0] for i in rows])
    combID = np.array([i[3] for i in rows])
    ew = np.array([i[4] for i in rows])
    datasetID = np.array([i[5] for i in rows])
    weight = np.array([i[6] for i in rows])
    advanced = np.array([i[7] for i in rows])
    totalSNR = np.array([i[8] for i in rows])

    tempscore = {"O": 0, "B": 10, "A":20, "F":30, "G":40}

    combIDlist = []
    for c in combID:
        if c not in combIDlist:
            combIDlist.append(c)

    sptypescore = []
    ebvlist = []
    totalSNRlist = []
    for c in combIDlist:
        objname_c = objname[combID == c][0]
        ebv_c = ebv[combID == c][0]
        sptype_c = sptype[combID == c][0]
        ew_c = ew[combID == c][0]
        totalSNR_c = totalSNR[combID == c][0]
        try:
            sptypescore.append(tempscore[sptype_c[0]] + int(sptype_c[1]))
            ebvlist.append(ebv_c)
            totalSNRlist.append(totalSNR_c)
        except:
            print("ERROR:", objname_c,sptype_c)


    cur.execute(
        "select od.registeredname,z.E_BV,z.sptype,dr.pipelineID,dr.totalSNR "
        "from object as z "
        "join objectdict as od using (objectID) "
        "join datareduction as dr using (objectID) "
        "where z.type='OBJECT' and dr.totalSNR>100 and dr.mode='WIDE' and od.priority=1 and z.E_BV is not Null;")

    rows = cur.fetchall()
    objname = np.array([i[0] for i in rows])
    ebv = np.array([i[1] for i in rows])
    sptype = np.array([i[2].split()[0] for i in rows])
    pipelineID = np.array([i[3] for i in rows])
    totalSNR = np.array([i[4] for i in rows])


    objnamelist = []
    for o in objname:
        if o not in objnamelist:
            objnamelist.append(o)

    sptypescore_all = []
    ebvlist_all = []
    totalSNRlist_all = []
    for c in objnamelist:
        objname_c = objname[objname == c][0]
        ebv_c = ebv[objname == c][0]
        sptype_c = sptype[objname == c][0]
        totalSNR_c = totalSNR[objname == c][0]
        try:
            sptypescore_all.append(tempscore[sptype_c[0]] + int(sptype_c[1]))
            ebvlist_all.append(ebv_c)
            totalSNRlist_all.append(totalSNR_c)
        except:
            print("ERROR:", objname_c, sptype_c)

    plt.figure()
    sc = plt.scatter(ebvlist_all, sptypescore_all, s=20, vmin=200, vmax=800, c=totalSNRlist_all, cmap=cm.winter_r, edgecolors=None, marker="^", label="Not analyzed")
    plt.scatter(ebvlist, sptypescore, s=35, c=totalSNRlist, vmin=200, vmax=800, cmap=cm.winter_r, edgecolors="k", linewidths=1, label="DIB analyzed")
    print("{0}/{1}".format(len(ebvlist), len(ebvlist_all)))
    cbar = plt.colorbar(sc)
    cbar.set_label("Signal-to-noise ratio")
    plt.legend()
    plt.xlabel("E(B-V)")
    plt.ylabel("Sptype")
    plt.yticks([5, 10, 15, 20], ["O5", "B0", "B5", "A0"])
    plt.ylim(1,25)
    # plt.xlim(-0.1, 2.0)
    plt.savefig("temporaly_files/ebv_sptype_snr_dist.png", format="png", dpi=200)
    plt.clf()

    plt.hist(totalSNRlist_all, bins=20)
    plt.xlabel("S/N")
    plt.savefig("temporaly_files/SNRhist.png")
    plt.clf()

    plt.hist(sptypescore_all, bins=50, range=(0,25))
    plt.xlabel("Sp type")
    plt.xticks([5, 10, 15, 20], ["O5", "B0", "B5", "A0"])
    plt.xlim(0,25)
    plt.savefig("temporaly_files/SpTypehist.png")

    conn.close()
