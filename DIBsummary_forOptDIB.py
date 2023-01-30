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
    objectIDinput = [522, 13,9,75,94,523,251,100,524,118,263,124,59,525,526,290,40,293,66,68,71,78,80,90,105]
    # dibid = args.dibid
    outputpdf = args.outputpdf
    # objectID.reverse()

    stylefile = "temporaly_files/stylefile.txt"

    objstr = ''
    for i in objectIDinput:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    conn, cur = openproject()

    dibid5780 = 160
    cur.execute(
        "select z.objectID,x.EW,x.EWerr,z.E_BV from DIBEWsummary as x "
        "join object as z using(objectID) where x.DIBID=%d and x.reference='Fan et al. (2019), ApJ, 878, 151';" % dibid5780)
    rows = cur.fetchall()
    objectid5780 = numpy.array([int(i[0]) for i in rows])
    ew5780, ewerr5780 = {}, {}
    for i in rows:
        ew5780[i[0]] = float(i[1])
        ewerr5780[i[0]] = float(i[2])
    DIBinfo = GetDIBLine(dibid5780)
    DIBwav5780 = DIBinfo[1]

    dibid5797 = 165
    cur.execute(
        "select z.objectID,x.EW,x.EWerr,z.E_BV from DIBEWsummary as x "
        "join object as z using(objectID) where x.DIBID=%d and x.reference='Fan et al. (2019), ApJ, 878, 151';" % dibid5797)
    rows = cur.fetchall()
    objectid5797 = numpy.array([int(i[0]) for i in rows])
    ew5797, ewerr5797 = {}, {}
    for i in rows:
        ew5797[i[0]] = float(i[1])
        ewerr5797[i[0]] = float(i[2])
    DIBinfo = GetDIBLine(dibid5797)
    DIBwav5797 = DIBinfo[1]

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

    for o in objectID:
        cur.execute("SELECT x.DIBID, x.EW, x.EWerr, x.upperlimit from DIBEWsummary as x join DIBlist as z using(DIBID) "
                    "where x.objectID=%d and x.reference='Fan et al. (2019), ApJ, 878, 151' and z.category NOT IN ('fake', 'Cs I');" % o)
        rows = cur.fetchall()

        if rows == []:
            print("DIBs are not measured for Object ID = %d" % objectID)
            sys.exit()
        else:
            DIBID = numpy.append(DIBID, numpy.array([i[0] for i in rows]))
            EW = numpy.append(EW, numpy.array([i[1] for i in rows]))
            EWerr = numpy.append(EWerr, numpy.array([i[2] for i in rows]))
            objIDref = numpy.append(objIDref, numpy.array([o for i in rows]))

    DIBinfo = GetDIBdict()
    [DIBIDlist, wav_air, reference, category, DIBfwhm, wavenumber, DIBcomment] = DIBinfo

    DIBset = list(set(DIBID))
    wav_set = [wav_air[i] for i in DIBset]
    DIBzip = zip(wav_set, DIBset)
    DIBzip_sorted = sorted(DIBzip)
    wav_set, DIBset = zip(*DIBzip_sorted)

    pp = PdfPages(outputpdf)

    for i in DIBset:
        print("Now plotting DIB %.1f" % wav_air[i])
        req = DIBID == i
        num = numpy.sum(req)

        EWset = EW[req]
        EWerrset = EWerr[req]
        oIDref = objIDref[req]

        fig = plt.figure(figsize=(12, 9))
        descax = plt.axes([0.06, 0.75, 0.51, 0.2])
        tableax = plt.axes([0.6, 0.05, 0.37, 0.9])
        ebvax = plt.axes([0.06, 0.45, 0.24, 0.28])
        dibax = plt.axes([0.36, 0.45, 0.21, 0.28])
        fwhmax = plt.axes([0.06, 0.1, 0.24, 0.28])
        vax = plt.axes([0.36, 0.1, 0.21, 0.28])
        # corax = plt.axes([])

        objIDpf, objIDpfdet, ebvpf, ewlist, ewerrlist, fwhmlist, heliovdict, ewdict, ewerrdict = [], [], [], [], [], [], {}, {}, {}
        for o in objectID:
            objIDpf.append(o)
            ebvpf.append(ebv[objectIDdb == o][0])
            ewlist.append(EWset[oIDref == o][0])
            ewerrlist.append(EWerrset[oIDref == o][0])
            ewdict[o] = EWset[oIDref == o][0]
            ewerrdict[o] = EWerrset[oIDref == o][0]

        objIDpf = numpy.array(objIDpf)
        ebvpf = numpy.array(ebvpf)
        ewlist = numpy.array(ewlist)
        ewerrlist = numpy.array(ewerrlist)

        column = ["ID", "Object", "E(B-V)", r"EW (m$\AA$)"]
        rowtexts = []
        for o in objectID:
            if EWset[oIDref == o][0] > 0.:
                rowtexts.append([str(o), registeredname[objectIDdb == o][0], str(ebv[objectIDdb == o][0]),
                                 "%.1f" % EWset[oIDref == o][0] + r"$\pm$" + "%.1f" % EWerrset[oIDref == o][0]])
            else:
                rowtexts.append([str(o), registeredname[objectIDdb == o][0], str(ebv[objectIDdb == o][0]),
                                 r"$<$" + "%.1f" % EWerrset[oIDref == o][0]])

        dibtb = tableax.table(cellText=rowtexts, colLabels=column, bbox=[0, 0, 1, 1])
        dibtb.auto_set_font_size(False)
        dibtb.set_fontsize(9)
        dibtb.auto_set_column_width(col=list(range(len(column))))
        tableax.axis("off")

        for m in range(len(column)):
            dibtb[0, m].set_facecolor('#363636')
            dibtb[0, m].set_text_props(color='w')

        r_ebv, p_ebv = ebv_correlation(ebvax, ebvpf, ewlist, ewerrlist, objIDpf, wav_air[i], stylefile=stylefile)
        r_dib5780, p_dib5780, _ = DIBcorrelation(dibax, ew5780, ewerr5780, objectid5780, DIBwav5780, ewdict, ewerrdict,
                                      objIDpf, wav_air[i],
                                      stylefile=stylefile)
        fwhmax.hist(fwhmlist)
        fwhmax.set_xlabel(r'DIB 5797 ($\AA$)')
        r_dib5797, p_dib5797, _ = DIBcorrelation(fwhmax, ew5797, ewerr5797, objectid5797, DIBwav5797, ewdict, ewerrdict,
                                      objIDpf, wav_air[i],
                                      stylefile=stylefile)

        descax.set_xlim(-1, 1)
        descax.set_ylim(-1, 1)
        descax.text(0, 0.4, r"$\lambda _{air}$ = " + "%.2f" % wav_air[i] + r"$\AA$, " + "FWHM = " + "%.2f" % DIBfwhm[
            i] + r"$\AA$" + "\nreference: %s, " % reference[i] + "category: %s\n" % category[i] + "comment: %s" %
                    DIBcomment[i], ha="center", va="center", fontsize=11)
        descax.text(0, -0.4,
                    "Detection rate: %d/%d objects (%.1f percent)\nFWHM (average): %.2fA\nCorrelation coeff. w/ E(B-V): %.2f (p=%.1e)\nCorrelation coeff. w/ DIB %d: %.2f (p=%.1e)" % (
                        numpy.sum(ewlist > 0.), len(objectID), (numpy.sum(ewlist > 0.) / len(objectID) * 100.),
                        numpy.average(fwhmlist), r_ebv, p_ebv, DIBwav5780, r_dib5780, p_dib5780), ha="center",
                    va="center",
                    fontsize=11)
        descax.xaxis.set_visible(False)
        descax.yaxis.set_visible(False)

        plt.suptitle("DIB %.1f (ID=%d)" % (wav_air[i], i))
        plt.savefig(pp, format="pdf")
        plt.clf()
        plt.close(fig)

    pp.close()

    conn.close()
