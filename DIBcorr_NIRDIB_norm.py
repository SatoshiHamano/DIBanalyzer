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
    parser.add_argument("-o", "--objectID", type=int, help="object ID", nargs='*')
    # parser.add_argument("-d", "--dibid", type=int, default=0, help="DIB ID")

    args = parser.parse_args()
    # objectIDinput = [37, 43, 50, 34, 90, 124, 12, 13, 14, 10, 11, 8, 23, 123, 55, 22, 59, 27, 246, 29, 30, 63, 31, 32,
    #                  33, 35, 39, 40, 41, 42, 66,
    #                  45, 75, 86, 52, 102, 128, 132, 133, 60, 44, 135, 28, 65, 24, 25, 121, 9, 15, 6, 144, 120, 149]
    objectIDinput = [50, 75, 33, 43, 30, 123, 31, 42, 124, 90, 246, 32, 27, 63, 45, 59, 41, 66, 40, 39,
                     14, 10, 12, 13, 15, 9, 132, 144]
    # dibid = args.dibid
    outputpdf = args.outputpdf
    # objectID.reverse()

    objstr = ''
    for i in objectIDinput:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    conn, cur = openproject()

    stylefile = "temporaly_files/stylefile.txt"

    [combineID, measurementID, DIBID, echelleorder, DIBspecpath, primaryflag, autonormalizeflag, automeasurementflag,
     EW, EWerr, centerlam_air, helio_velocity, FWHM, SNR, integration_start, integration_end, comment, depth,
     objIDref] = [numpy.array([]) for i in range(19)]

    cur.execute(
        "select objectID,E_BV,sptype,registeredname from object join objectdict using(objectID) "
        "where objectID IN (%s) and priority=1 order by E_BV;" % objstr)
    rows = cur.fetchall()
    objectID = [int(i[0]) for i in rows]
    objectIDdb = numpy.array([int(i[0]) for i in rows])
    ebv = numpy.array([i[1] for i in rows])
    sptype = numpy.array([i[2] for i in rows])
    registeredname = numpy.array([i[3] for i in rows])
    print(objectID, registeredname, ebv)

    for o in range(len(objectID)):
        cur.execute("SELECT x.combineID, x.measurementID, x.DIBID, x.echelleorder, x.DIBspecpath, "
                    "x.primaryflag, x.autonormalizeflag, x.automeasurementflag, x.EW, x.EWerr, x.centerlam_air, "
                    "x.helio_velocity, x.FWHM, x.SNR, x.integration_start, x.integration_end, x.comment, x.depth "
                    "from DIBmeasurement as x join combinesummary as y using(combineID) join DIBlist as z using(DIBID) "
                    "where y.objectID=%d and z.category NOT IN ('fake', 'Cs I');" % objectID[o])
        rows = cur.fetchall()

        if rows == []:
            print("DIBs are not measured for Object ID = %d" % objectID[o])
            sys.exit()
        else:
            print(objectID[o], registeredname[o], ebv[o])
            combineID = numpy.append(combineID, numpy.array([i[0] for i in rows]))
            measurementID = numpy.append(measurementID, numpy.array([i[1] for i in rows]))
            DIBID = numpy.append(DIBID, numpy.array([i[2] for i in rows]))
            echelleorder = numpy.append(echelleorder, numpy.array([i[3] for i in rows]))
            DIBspecpath = numpy.append(DIBspecpath, numpy.array([i[4] for i in rows]))
            primaryflag = numpy.append(primaryflag, numpy.array([i[5] for i in rows]))
            autonormalizeflag = numpy.append(autonormalizeflag, numpy.array([i[6] for i in rows]))
            automeasurementflag = numpy.append(automeasurementflag, numpy.array([i[7] for i in rows]))
            EW = numpy.append(EW, numpy.array([i[8] / ebv[o] for i in rows]))
            EWerr = numpy.append(EWerr, numpy.array([i[9] / ebv[o] for i in rows]))
            centerlam_air = numpy.append(centerlam_air, numpy.array([i[10] for i in rows]))
            helio_velocity = numpy.append(helio_velocity, numpy.array([i[11] for i in rows]))
            FWHM = numpy.append(FWHM, numpy.array([i[12] for i in rows]))
            SNR = numpy.append(SNR, numpy.array([i[13] for i in rows]))
            integration_start = numpy.append(integration_start, numpy.array([i[14] for i in rows]))
            integration_end = numpy.append(integration_end, numpy.array([i[15] for i in rows]))
            comment = numpy.append(comment, numpy.array([i[16] for i in rows]))
            depth = numpy.append(depth, numpy.array([i[17] for i in rows]))
            objIDref = numpy.append(objIDref, numpy.array([objectID[o] for i in rows]))


    DIBinfo = GetDIBdict()
    [DIBIDlist, wav_air, reference, category, DIBfwhm, wavenumber, DIBcomment] = DIBinfo

    DIBset = list(set(DIBID))
    wav_set = [wav_air[i] for i in DIBset]
    DIBzip = zip(wav_set, DIBset)
    DIBzip_sorted = sorted(DIBzip)
    wav_set, DIBset = zip(*DIBzip_sorted)

    pp = PdfPages(outputpdf)
    plt.figure(figsize=(8, 8))

    for i in DIBset:#[33]:
        print("Now plotting DIB %.1f" % wav_air[i])
        req = DIBID == i
        num = numpy.sum(req)

        cID = combineID[req]
        tstatus = numpy.array([telStatus(cur, c) for c in cID])
        pf = primaryflag[req]
        EWset = EW[req]
        EWerrset = EWerr[req]
        SNRset = SNR[req]
        oIDref = objIDref[req]

        objIDpf1, ewdict1, ewerrdict1 = [], {}, {}
        for n in range(len(objectID)):
            if 1 in pf[oIDref == objectID[n]]:
                repID = (oIDref == objectID[n]) & (pf == 1)
                objIDpf1.append(objectID[n])
                ewdict1[objectID[n]] = EWset[repID][0]
                ewerrdict1[objectID[n]] = EWerrset[repID][0]

        objIDpf1 = numpy.array(objIDpf1)

        for j in DIBset:#[33,34,35,37,38]:#DIBset:
            if i != j:
                req = DIBID == j
                num = numpy.sum(req)

                cID = combineID[req]
                tstatus = numpy.array([telStatus(cur, c) for c in cID])
                pf = primaryflag[req]
                EWset = EW[req]
                EWerrset = EWerr[req]
                SNRset = SNR[req]
                oIDref = objIDref[req]

                objIDpf2, ewdict2, ewerrdict2 = [], {}, {}
                for n in range(len(objectID)):
                    if 1 in pf[oIDref == objectID[n]]:
                        repID = (oIDref == objectID[n]) & (pf == 1)
                        objIDpf2.append(objectID[n])
                        ewdict2[objectID[n]] = EWset[repID][0]
                        ewerrdict2[objectID[n]] = EWerrset[repID][0]

                objIDpf2 = numpy.array(objIDpf2)

                objidcom = numpy.intersect1d(objIDpf1, objIDpf2)
                ew1com = numpy.array([ewdict1[i] for i in objidcom])
                ew2com = numpy.array([ewdict2[i] for i in objidcom])

                detection12 = (ew1com != 0.) & (ew2com != 0.)
                if numpy.sum(detection12) >= 10:
                    dibax = plt.axes()
                    r_dib, p_dib, n_dib, _, _ = DIBcorrelation(dibax, ewdict1, ewerrdict1, objIDpf1, wav_air[i],
                                                               ewdict2,
                                                               ewerrdict2, objIDpf2, wav_air[j], stylefile=stylefile)

                    if r_dib > 0.8:
                        plt.title("{} - {} (R={:.3f} (N={}))".format(wav_air[i], wav_air[j], r_dib, n_dib))
                        plt.savefig(pp, format="pdf")
                    plt.clf()

    pp.close()

    conn.close()
