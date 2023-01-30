# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import sys

if __name__ == '__main__':
    conn, cur = openproject()

    objectID = int(sys.argv[1])
    DIBID = int(sys.argv[2])

    cur.execute(
        "select x.measurementID, x.primaryflag from DIBmeasurement as x join combinesummary as y using(combineID) where y.objectID=%d and x.DIBID=%d;" % (
        objectID, DIBID))
    rows = cur.fetchall()
    if rows == []:
        print("No record.")
        sys.exit()
    else:
        IDlist = [i[0] for i in rows]
        primaryflag = [i[1] for i in rows]

    if not 1 in primaryflag:
        print("Primaryflag was not set.")
        sys.exit()

    for i in IDlist:
        cur.execute("update DIBmeasurement set primaryflag=0 where measurementID='%s';" % i)

    conn.commit()

    conn.close()
