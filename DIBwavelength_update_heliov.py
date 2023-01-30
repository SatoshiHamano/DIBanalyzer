# -*- coding:utf-8 -*-

from open_mysql_project import openproject
from add_lineDIB_mysql import GetDIBLine
import argparse
import scipy.constants
import numpy


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("wavelength_air", type=float, help="Updated DIB wavelength in air")
    parser.add_argument("-d", "--DIBID", type=int, help="DIB ID")
    parser.add_argument("-o", "--objectID", type=int, help="object ID", nargs='*')

    args = parser.parse_args()
    newwav = args.wavelength_air
    objectID = args.objectID
    dibid = args.DIBID

    conn, cur = openproject()

    lightv = scipy.constants.c * 1.e-3

    for o in objectID:
        cur.execute("SELECT z.registeredname, x.EW, x.EWerr, x.centerlam_air, x.helio_velocity "
                    "from DIBmeasurement as x join combinesummary as y using(combineID) join objectdict as z using(objectID) "
                    "where y.objectID=%d and x.primaryflag=1 and z.priority = 1 and x.DIBID=%d;" % (o, dibid))
        rows = cur.fetchall()
        if rows == []:
            print("No measurement for objectID={0}.".format(o))
        else:
            name = rows[0][0]
            ew = rows[0][1]
            ewerr = rows[0][2]
            lamair = rows[0][3]
            hvel = rows[0][4]
            if ew == 0.:
                print("ID={0}: {1}, < {2} mA, {3} A, {4} km/s".format(o, name, ewerr, lamair, hvel))
            else:
                newvel = (lamair - newwav) / newwav * lightv
                print("ID={0}: {1}, {2} pm {3} mA, {4} A, {5} km/s -> {6:.2f} km/s".format(o, name, ew, ewerr, lamair, hvel, newvel))

    conn.close()