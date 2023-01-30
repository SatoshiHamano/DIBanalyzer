# -*- coding:utf-8 -*-

from open_mysql_project import openproject

if __name__ == '__main__':
    conn, cur = openproject()

    rf = open("temporaly_files/Friedman2011_EBV.txt")
    rl = rf.readlines()
    rf.close()

    objid = [int(i.split()[0]) for i in rl]
    sptype = [i.split()[1] for i in rl]
    ebv = [float(i.split()[2]) for i in rl]
    HI = [i.split()[3] for i in rl]
    HIerr = [i.split()[4] for i in rl]
    H2 = [i.split()[5] for i in rl]
    H2err = [i.split()[6] for i in rl]

    for i in range(len(objid)):
        if ebv[i] != "blank":
            cur.execute("UPDATE object set E_BV=%.3f where objectid=%d;" % (ebv[i], objid[i]))
        if HI[i] != "blank":
            cur.execute("UPDATE object set log_N_HI=%s where objectid=%d;" % (HI[i], objid[i]))
        if HIerr[i] != "blank":
            cur.execute("UPDATE object set log_N_HIerr=%s where objectid=%d;" % (HIerr[i], objid[i]))
        if H2[i] != "blank":
            cur.execute("UPDATE object set log_N_HI=%s where objectid=%d;" % (H2[i], objid[i]))
        if H2err[i] != "blank":
            cur.execute("UPDATE object set log_N_HIerr=%s where objectid=%d;" % (H2err[i], objid[i]))

    conn.commit()
    conn.close()