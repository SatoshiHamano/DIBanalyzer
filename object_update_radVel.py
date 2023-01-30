#!/usr/bin/env python
# -*- coding: utf-8 -*-

from object_add_mysql import simbad_query
from open_mysql_project import openproject

if __name__ == '__main__':
    conn, cur = openproject()

    cur.execute("SELECT objectid, objectname from object;")
    rows = cur.fetchall()
    oid = [i[0] for i in rows]
    oname = [i[1] for i in rows]

    counter = 0
    for i in range(len(oid)):
        counter += 1
        result, flagsimbad = simbad_query(oname[i])
        radvel = result[10]
        if radvel != "0" and radvel != "Null":
            cur.execute("UPDATE object set radial_velocity=%.3f where objectid=%d;" % (radvel, oid[i]))
            print("{}/{} ID={} name={} V={}km/s".format(counter, len(oid), oid[i], oname[i], radvel))
        else:
            print("{}/{} ID={} name={} V=INDEF".format(counter, len(oid), oid[i], oname[i]))

        conn.commit()
    conn.close()