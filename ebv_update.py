#!/usr/bin/env python
# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import sys

if __name__ == '__main__':
    conn, cur = openproject()

    objectid = int(sys.argv[1])
    ebv = float(sys.argv[2])
    reference = sys.argv[3]

    cur.execute("SELECT E_BV, E_BV_reference from object where objectid=%d;" % objectid)
    rows = cur.fetchall()
    ebv_old = rows[0][0]
    reference_old = rows[0][1]

    cur.execute("UPDATE object set E_BV=%.2f where objectID=%d;" % (ebv, objectid))
    print("E_BV is updated from %s to %.2f." % (ebv_old, ebv))
    cur.execute("UPDATE object set E_BV_reference='%s' where objectID=%d;" % (reference, objectid))
    print("E_BV_reference is updated from %s to %s." % (reference_old, reference))

    conn.commit()
    conn.close()