# -*- coding:utf-8 -*-

# python DIBprimary.py temporaly_files/CygOB2cluster_DIBprimary.pdf -o 12 14 13 10 15 11 9

from open_mysql_project import openproject
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from spectra_plotter import MultiSpecPlotter
import numpy
import sys, os
from add_lineDIB_mysql import GetDIBdict, GetDIBLine
from spectra_plotter import spectrum_plot
from DIBanalysis import openDIBresult
from Spec1Dtools import openspecfits
from waveshift_measure import atran_resampling
from vac2air_spec import vac2air
from correlation_ebv import ebv_correlation
from correlation_dibpair import DIBcorrelation, Vcorrelation


def telStatus(cur, combineID):
    cur.execute(
        "select telluricflag from combinesummary where combineID='%s';" % combineID)
    rows = cur.fetchall()
    telflag = rows[0][0]
    if telflag == 1:
        cur.execute("select t.advanced from combinedataset as c "
                    "join telluriccorrection as t on c.datasetID = t.telluricID "
                    "where c.combineID = '%s';" % combineID)
        rows = cur.fetchall()
        advanced = [i[0] for i in rows]
        if all(advanced):
            tstatus = "AT"
        elif any(advanced):
            tstatus = "mixT"
        else:
            tstatus = "T"
    else:
        tstatus = "R"

    return tstatus


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outputpdf", type=str, help="Output pdf")
    # parser.add_argument("-d", "--dibid", type=int, default=0, help="DIB ID")

    args = parser.parse_args()
    # objectIDinput = [522, 13,9,75,94,523,251,100,524,118,263,124,59,525,526,290,40,293,66,68,71,78,80,90,105]
    # objectIDinput = [522,75,94,523,251,100,524,118,263,124,59,525,526,290,40,293,66,68,71,78,80,90,105]
    objectIDinput = [522, 13,9,75,94,523,251,100,524,118,263,124,59,525,526,290,293,66,68,71,78,80,90,105]
    # dibid = args.dibid
    outputpdf = args.outputpdf
    # objectID.reverse()

    stylefile = "temporaly_files/stylefile.txt"
    EWref1 = 'Galazutdinov et al. (2020), ApJ, 159, 113'
    EWref2 = 'Fan et al. (2019), ApJ, 878, 151'

    objstr = ''
    for i in objectIDinput:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    conn, cur = openproject()


    [combineID, measurementID, DIBID, echelleorder, DIBspecpath, primaryflag, autonormalizeflag, automeasurementflag,
     EW, EWerr, centerlam_air, helio_velocity, FWHM, SNR, integration_start, integration_end, comment, depth,
     objIDref] = [numpy.array([]) for i in range(19)]

    cur.execute(
        "select objectID,E_BV,sptype,registeredname from object join objectdict using(objectID) "
        "where objectID IN (%s) and priority=1 order by E_BV;" % objstr)
    rows = cur.fetchall()
    objectID = [int(i[0]) for i in rows]
    objectIDdb = numpy.array([int(i[0]) for i in rows])
    ebv = numpy.array([float(i[1]) for i in rows])
    sptype = numpy.array([i[2] for i in rows])
    registeredname = numpy.array([i[3] for i in rows])

    DIBlist = [93, 101,123, 583,619] #broad DIBs
    # DIBlist = [160, 583]
    # DIBlist = [117,146,125,141,242,261,315,407,440,583,623,628,160]
    # DIBlist = [117,146,583,160,623,261]
    # DIBlist = [125,141,160,623,261,440]
    # DIBlist = [583,160,242,628,315,407]


    DIBinfo = GetDIBdict()
    [DIBIDlist, wav_air, reference, category, DIBfwhm, wavenumber, DIBcomment] = DIBinfo

    DIBset = list(set(DIBID))
    wav_set = [wav_air[i] for i in DIBlist]

    pp = PdfPages(outputpdf)
    plt.figure(figsize=(8,8))

    for i in range(len(DIBlist)):
        # cur.execute(
        #     "SELECT x.DIBID, x.EW, x.EWerr, x.upperlimit, x.objectID from DIBEWsummary as x join DIBlist as z using(DIBID) "
        #     "where x.DIBID = {} and x.objectID IN ({}) and x.reference='{}' and z.category NOT IN ('fake', 'Cs I');".format(
        #     DIBlist[i], objstr, reference))
        cur.execute(
            "SELECT x.DIBID, x.EW, x.EWerr, x.upperlimit, x.objectID from DIBEWsummary as x join DIBlist as z using(DIBID) "
            "where x.DIBID = {} and x.reference='{}' and z.category NOT IN ('fake', 'Cs I');".format(
            DIBlist[i], EWref1))
        rows = cur.fetchall()
        DIBID1 = numpy.array([k[0] for k in rows])
        EW1 = numpy.array([k[1] for k in rows])
        EWerr1 = numpy.array([k[2] for k in rows])
        objIDref1 = numpy.array([k[4] for k in rows])
        EWdict1 = {}
        EWerrdict1 = {}
        for n in range(len(DIBID1)):
            EWdict1[objIDref1[n]] = EW1[n]
            EWerrdict1[objIDref1[n]] = EWerr1[n]
        # print(EWdict1)
        for j in range(len(DIBlist)):
            # cur.execute(
            #     "SELECT x.DIBID, x.EW, x.EWerr, x.upperlimit, x.objectID from DIBEWsummary as x join DIBlist as z using(DIBID) "
            #     "where x.DIBID = {} and x.objectID IN ({}) and x.reference='{}' and z.category NOT IN ('fake', 'Cs I');".format(DIBlist[j], objstr, reference))
            cur.execute(
                "SELECT x.DIBID, x.EW, x.EWerr, x.upperlimit, x.objectID from DIBEWsummary as x join DIBlist as z using(DIBID) "
                "where x.DIBID = {} and x.reference='{}' and z.category NOT IN ('fake', 'Cs I');".format(DIBlist[j], EWref2))
            rows = cur.fetchall()
            DIBID2 = numpy.array([k[0] for k in rows])
            EW2 = numpy.array([k[1] for k in rows])
            EWerr2 = numpy.array([k[2] for k in rows])
            objIDref2 = numpy.array([k[4] for k in rows])
            EWdict2 = {}
            EWerrdict2 = {}
            for n in range(len(DIBID2)):
                EWdict2[objIDref2[n]] = EW2[n]
                EWerrdict2[objIDref2[n]] = EWerr2[n]
            # print(EWdict2)

            dibax = plt.axes()
            r_dib, p_dib, n_dib, _, _ = DIBcorrelation(dibax, EWdict1, EWerrdict1, objIDref1, wav_air[DIBlist[i]], EWdict2,
                                                  EWerrdict2, objIDref2, wav_air[DIBlist[j]], stylefile=stylefile)

            plt.title("{} - {} (R={:.3f} (N={}))".format(wav_air[DIBlist[i]], wav_air[DIBlist[j]], r_dib, n_dib))
            plt.savefig(pp, format="pdf")
            plt.clf()

    pp.close()

    conn.close()
