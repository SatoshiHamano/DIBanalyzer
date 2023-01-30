# -*- coding:utf-8 -*-

import sys
from open_mysql_project import openproject

if __name__ == '__main__':
    conn, cur = openproject()

    cur.execute("SELECT x.E_BV,y.totalSNR,x.Bmag-x.Vmag,x.objectID,z.registeredname,x.sptype from object as x "
                "join datareduction as y using(objectID) "
                "join objectdict as z using (objectID) where type=\"OBJECT\" and z.priority=1;")
    rows = cur.fetchall()
    ebv = [i[0] for i in rows]
    snr = [i[1] for i in rows]
    bv = [i[2] for i in rows]
    objid = [i[3] for i in rows]
    objname = [i[4] for i in rows]
    sptype = [i[5] for i in rows]

    lookedids = []
    for i in range(len(ebv)):
        if ebv[i] == None and objid[i] not in lookedids and bv[i] != None:
            cur.execute("SELECT reference from DIBEWsummary where objectID=%d;" % objid[i])
            rows = cur.fetchall()
            ref = list(set([j[0] for j in rows]))
            lookedids.append(objid[i])
            print(objid[i],"\t", objname[i], "\t", sptype[i], bv[i], snr[i], ref)

    conn.close()