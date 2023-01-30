# -*- coding:utf-8 -*-

from open_mysql_project import openproject
import sys, os
import numpy
from add_lineDIB_mysql import GetDIBList

class Color:
    BLACK     = '\033[30m'
    RED       = '\033[31m'
    GREEN     = '\033[32m'
    YELLOW    = '\033[33m'
    BLUE      = '\033[34m'
    PURPLE    = '\033[35m'
    CYAN      = '\033[36m'
    WHITE     = '\033[37m'
    END       = '\033[0m'
    BOLD      = '\038[1m'
    UNDERLINE = '\033[4m'
    INVISIBLE = '\033[08m'
    REVERCE   = '\033[07m'

if __name__ == '__main__':
    conn, cur = openproject()

    objectID = int(sys.argv[1])

    cur.execute("SELECT x.combineID, x.measurementID, x.DIBID, x.echelleorder, x.DIBspecpath, "
                "x.primaryflag, x.autonormalizeflag, x.automeasurementflag, x.EW, x.EWerr, x.centerlam_air, "
                "x.helio_velocity, x.FWHM, x.SNR, x.integration_start, x.integration_end, x.comment, x.depth "
                "from DIBmeasurement as x join combinesummary as y using(combineID) where y.objectID=%d;" % objectID)
    rows = cur.fetchall()

    if rows == []:
        print("DIBs are not measured for Object ID = %d" % objectID)
        sys.exit()
    else:
        combineID = numpy.array([i[0] for i in rows])
        measurementID = numpy.array([i[1] for i in rows])
        DIBID = numpy.array([i[2] for i in rows])
        echelleorder = numpy.array([i[3] for i in rows])
        DIBspecpath = numpy.array([i[4] for i in rows])
        primaryflag = numpy.array([i[5] for i in rows])
        autonormalizeflag = numpy.array([i[6] for i in rows])
        automeasurementflag = numpy.array([i[7] for i in rows])
        EW = numpy.array([i[8] for i in rows])
        EWerr = numpy.array([i[9] for i in rows])
        centerlam_air = numpy.array([i[10] for i in rows])
        helio_velocity = numpy.array([i[11] for i in rows])
        FWHM = numpy.array([i[12] for i in rows])
        SNR = numpy.array([i[13] for i in rows])
        integration_start = numpy.array([i[14] for i in rows])
        integration_end = numpy.array([i[15] for i in rows])
        comment = numpy.array([i[16] for i in rows])
        depth = numpy.array([i[17] for i in rows])

    DIBinfo = GetDIBList(9000, 13500.)
    [DIBset, wav_set, reference, category, fwhm, wavenumber, DIBcomment] = DIBinfo

    # DIBset = list(set(DIBIDlist))
    # wav_set = [wav_air[i] for i in DIBset]
    # DIBzip = zip(wav_set, DIBset)
    # DIBzip_sorted = sorted(DIBzip)
    # wav_set, DIBset = zip(*DIBzip_sorted)

    print("\nDIB status for objectID=%d:\n" % objectID)

    for n in range(len(DIBset)):
        cur.execute(
            "SELECT x.measurementID, x.primaryflag from DIBmeasurement as x join combinesummary as y using(combineID) " +
            "where y.objectID=%d and x.DIBID=%d;" % (objectID, DIBset[n]))
        rows = cur.fetchall()

        if rows == []:
            print(Color.GREEN + "DIB%.1f (ID:%d): not measured" % (wav_set[n], DIBset[n]) + Color.END)
            continue
        else:
            measurementID = numpy.array([i[0] for i in rows])
            primaryflag = numpy.array([i[1] for i in rows])

        primaryNum = numpy.sum(primaryflag)
        measureNum = len(measurementID)
        if primaryNum == 1:
            print(Color.BLACK + "DIB%.1f (ID:%d): %d times measured, primary='%s'" % (wav_set[n], DIBset[n], measureNum, measurementID[primaryflag==1][0]) + Color.END)
        elif primaryNum > 1:
            print(Color.RED + "DIB%.1f (ID:%d): %d times measured, multiple primary flags." % (wav_set[n], DIBset[n], measureNum) + Color.END)
        else:
            print(Color.BLUE + "DIB%.1f (ID:%d): %d times measured, no primary flag was set." % (wav_set[n], DIBset[n], measureNum) + Color.END)


    conn.close()
