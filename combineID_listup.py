# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import numpy as np

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

    combIDlist = []
    for c in combID:
        if c not in combIDlist:
            combIDlist.append(c)

    for c in combIDlist:
        objname_c = objname[combID == c][0]
        ebv_c = ebv[combID == c][0]
        sptype_c = sptype[combID == c][0]
        ew_c = ew[combID == c][0]
        datasetID_c = datasetID[combID == c]
        weight_c = weight[combID == c]
        advanced_c = advanced[combID == c]
        totalSNR_c = totalSNR[combID == c]
        datasetIDlist = list(set(datasetID_c))
        cnum = len(datasetIDlist)
        righttext = "%s(%.2f, %s)\\n13175 EW = %.1f mA\\n" % (objname_c, ebv_c, sptype_c, ew_c)
        for i in range(cnum):
            righttext += "%s(w=%.2f,snr=%.0f,a=%d)\\n" % (datasetID_c[i][:10],weight_c[i],totalSNR_c[i],advanced_c[i])
        print(objname_c, "|", ebv_c, "|", sptype_c)
        # print(c,"|", righttext)

    conn.close()
