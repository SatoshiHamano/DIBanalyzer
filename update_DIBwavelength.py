#!/usr/bin/env python
# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import argparse
from add_lineDIB_mysql import *
import numpy as np
import scipy.constants

if __name__ == '__main__':
    c = scipy.constants.c * 1.e-3
    conn, cur = openproject()

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dibid", type=int, help="dibid")
    parser.add_argument("-w", "--wavelength", type=float, help="New wavelength")
    parser.add_argument("-t", "--type", choices=["AIR", "VAC"])


    args = parser.parse_args()

    DIBID = args.dibid
    newwav = args.wavelength
    vacorair = args.type

    [DIBID_db, wav, reference_db, category_db, fwhm_db, wavenumber_db, comment_db] = GetDIBLine(DIBID)
    update = UpdateDIBDB(DIBID, lam=newwav, vacorair=vacorair, confirm=False)
    if not update:
        conn.close()
        sys.exit()

    cur.execute("SELECT measurementID, centerlam_air, centerlam_vac, helio_velocity from DIBmeasurement where DIBID = '{}';".format(DIBID))
    rows = cur.fetchall()
    if rows == []:
        print("No measurements.")
    else:
        mID = np.array([i[0] for i in rows])
        clam_air = np.array([float(i[1]) for i in rows])
        clam_vac = np.array([float(i[2]) for i in rows])
        hv = np.array([float(i[3]) for i in rows])
        print("{} records were found.".format(len(mID)))
        if vacorair == "VAC":
            newhv = (clam_vac - newwav) / newwav * c
        elif vacorair == "AIR":
            newhv = (clam_air - newwav) / newwav * c

        print("DIFF: {:.3f} km/s".format(np.average(newhv - hv)))
        for i in range(len(rows)):
            cur.execute("UPDATE DIBmeasurement SET helio_velocity={:.3f} where measurementID='{}';".format(newhv[i], mID[i]))

        conn.commit()

    conn.close()