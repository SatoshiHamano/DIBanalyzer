# -*- coding:utf-8 -*-

from open_mysql_project import openproject
from DIBanalysis import obtainMultiIDsFromMeasurementID
import sys

if __name__ == '__main__':
    conn, cur = openproject()

    measurementID = sys.argv[1]
    # # objectID = int(sys.argv[2])
    # # combineID = measurementID.split("_DIB")[0]
    # # DIBID = int(measurementID.split("_DIB")[1].split("_m")[0])
    #
    # cur.execute("select x.DIBID, y.objectID, x.combineID from DIBmeasurement as x join combinesummary as y "
    #             "using(combineID) where measurementID='%s';" % measurementID)
    # rows = cur.fetchall()
    # DIBID = int(rows[0][0])
    # objectID = int(rows[0][1])
    # combineID = rows[0][2]
    # print("DIB ID: %d" % DIBID)
    # print("object ID: %d" % objectID)
    # print("combine ID: %s" % combineID)

    DIBID, objectID, combineID, _ = obtainMultiIDsFromMeasurementID(measurementID)

    cur.execute(
        "select x.measurementID, x.primaryflag from DIBmeasurement as x join combinesummary as y using(combineID) where y.objectID=%d and x.DIBID=%d;" % (
        objectID, DIBID))
    rows = cur.fetchall()
    IDlist = [i[0] for i in rows]
    primaryflag = [i[1] for i in rows]

    if not measurementID in IDlist:
        print("%s is not found in the database." % measurementID)
        sys.exit()

    if 1 in primaryflag:
        print("Caution: primaryflag was already set.")
        for i in range(len(IDlist)):
            print(IDlist[i], primaryflag[i])

    for i in IDlist:
        if i == measurementID:
            cur.execute("update DIBmeasurement set primaryflag=1 where measurementID='%s';" % i)
        else:
            cur.execute("update DIBmeasurement set primaryflag=0 where measurementID='%s';" % i)

    conn.commit()

    print("'%s' was set as the primary measurement result." % measurementID)

    conn.close()
