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
    objectIDinput = args.objectID
    # dibid = args.dibid
    outputpdf = args.outputpdf
    # objectID.reverse()

    maxobjn = 11
    atranmodelfile = "atran.smo.11513_R28000_0ft.npz"
    samplespec = "WINERED_sample_spectrum_WIDE.fits"
    stylefile = "temporaly_files/stylefile.txt"
    telwav_vac, telflux = atran_resampling(atranmodelfile, samplespec)
    telwav_air = vac2air(telwav_vac)

    objstr = ''
    for i in objectIDinput:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    conn, cur = openproject()

    dibid13175 = 38
    cur.execute(
        "select z.objectID,x.EW,x.EWerr,z.E_BV from DIBmeasurement as x join combinesummary as y using(combineID) "
        "join object as z using(objectID) where x.DIBID=%d and x.primaryflag = 1;" % dibid13175)
    rows = cur.fetchall()
    objectid13175 = numpy.array([int(i[0]) for i in rows])
    ew13175, ewerr13175 = {}, {}
    for i in rows:
        ew13175[i[0]] = float(i[1])
        ewerr13175[i[0]] = float(i[2])
    DIBinfo = GetDIBLine(dibid13175)
    DIBwav13175 = DIBinfo[1]

    dibid10780 = 35
    cur.execute(
        "select z.objectID,x.helio_velocity from DIBmeasurement as x join combinesummary as y using(combineID) "
        "join object as z using(objectID) where x.DIBID=%d and x.primaryflag = 1;" % dibid10780)
    rows = cur.fetchall()
    objectid10780 = numpy.array([int(i[0]) for i in rows])
    heliov10780 = {}
    for i in rows:
        heliov10780[i[0]] = float(i[1])
    DIBinfo = GetDIBLine(dibid10780)
    DIBwav10780 = DIBinfo[1]

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
        cur.execute("SELECT x.combineID, x.measurementID, x.DIBID, x.echelleorder, x.DIBspecpath, "
                    "x.primaryflag, x.autonormalizeflag, x.automeasurementflag, x.EW, x.EWerr, x.centerlam_air, "
                    "x.helio_velocity, x.FWHM, x.SNR, x.integration_start, x.integration_end, x.comment, x.depth "
                    "from DIBmeasurement as x join combinesummary as y using(combineID) join DIBlist as z using(DIBID) "
                    "where y.objectID=%d and z.category NOT IN ('fake', 'Cs I');" % o)
        rows = cur.fetchall()

        if rows == []:
            print("DIBs are not measured for Object ID = %d" % objectID)
            sys.exit()
        else:
            combineID = numpy.append(combineID, numpy.array([i[0] for i in rows]))
            measurementID = numpy.append(measurementID, numpy.array([i[1] for i in rows]))
            DIBID = numpy.append(DIBID, numpy.array([i[2] for i in rows]))
            echelleorder = numpy.append(echelleorder, numpy.array([i[3] for i in rows]))
            DIBspecpath = numpy.append(DIBspecpath, numpy.array([i[4] for i in rows]))
            primaryflag = numpy.append(primaryflag, numpy.array([i[5] for i in rows]))
            autonormalizeflag = numpy.append(autonormalizeflag, numpy.array([i[6] for i in rows]))
            automeasurementflag = numpy.append(automeasurementflag, numpy.array([i[7] for i in rows]))
            EW = numpy.append(EW, numpy.array([i[8] for i in rows]))
            EWerr = numpy.append(EWerr, numpy.array([i[9] for i in rows]))
            centerlam_air = numpy.append(centerlam_air, numpy.array([i[10] for i in rows]))
            helio_velocity = numpy.append(helio_velocity, numpy.array([i[11] for i in rows]))
            FWHM = numpy.append(FWHM, numpy.array([i[12] for i in rows]))
            SNR = numpy.append(SNR, numpy.array([i[13] for i in rows]))
            integration_start = numpy.append(integration_start, numpy.array([i[14] for i in rows]))
            integration_end = numpy.append(integration_end, numpy.array([i[15] for i in rows]))
            comment = numpy.append(comment, numpy.array([i[16] for i in rows]))
            depth = numpy.append(depth, numpy.array([i[17] for i in rows]))
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

        fig = plt.figure(figsize=(12, 9))
        ax = [plt.axes([0.06, 0.25, 0.2, 0.7]), plt.axes([0.39, 0.25, 0.2, 0.7]), plt.axes([0.72, 0.25, 0.2, 0.7])]
        telax = [plt.axes([0.06, 0.1, 0.2, 0.1]), plt.axes([0.39, 0.1, 0.2, 0.1]), plt.axes([0.72, 0.1, 0.2, 0.1])]

        cID = combineID[req]
        tstatus = numpy.array([telStatus(cur, c) for c in cID])
        mID = measurementID[req]
        eo = echelleorder[req]
        sp = DIBspecpath[req]
        pf = primaryflag[req]
        anf = autonormalizeflag[req]
        amf = automeasurementflag[req]
        EWset = EW[req]
        EWerrset = EWerr[req]
        cl = centerlam_air[req]
        hv = helio_velocity[req]
        FWHMset = FWHM[req]
        SNRset = SNR[req]
        ints = integration_start[req]
        inte = integration_end[req]
        com = comment[req]
        depthset = depth[req]
        oIDref = objIDref[req]

        repIDdict = {}
        for j in range(3):
            objIDax = objectID[j * 11: min(len(objectID), (j + 1) * 11)]
            objIDax.reverse()

            spxlist = []
            spylist = []
            lammin = []
            lammax = []
            yshift = []
            yshift_factor = []
            intlist = []
            pflist = []
            colorlist = []

            for n in range(len(objIDax)):
                if 1 in pf[oIDref == objIDax[n]]:
                    repID = (oIDref == objIDax[n]) & (pf == 1)
                    pflist.append(True)
                    # colorlist += ["0.5", "0.5", "b", "k", "k"]
                    colorlist += ["0.5", "0.5", "k"]
                else:
                    SNRmax = numpy.amax(SNRset[(oIDref == objIDax[n]) & (tstatus == "AT")])
                    repID = (oIDref == objIDax[n]) & (tstatus == "AT") & (SNRset == SNRmax)
                    pflist.append(False)
                    # colorlist += ["0.25", "0.25", "0.5", "0.5", "0.5"]
                    colorlist += ["0.25", "0.25", "0.5"]

                repIDdict[objIDax[n]] = repID

                DIBdir = os.path.dirname(sp[repID][0]) + "/"
                resultnpz = sp[repID][0].rstrip("fits") + "npz"
                results = openDIBresult(resultnpz)
                [mask, spx, spy_interp, sn, error, error_interp, error_cont] = results
                intrange = (spx > ints[repID][0] - 1.e-6) & (spx < inte[repID][0] + 1.e-6)
                onespec = numpy.ones(spx.shape)
                _, spy, _, _, _ = openspecfits(sp[repID][0])
                # if n == 0:
                #     if depthset[n] != 0.:
                #         ylow = 1. - max(depthset[n] * 1.3, 1. / SNRset[n] * 10.)
                #     else:
                #         ylow = 1. - 10. / SNRset[n]
                # else:
                if depthset[repID][0] != 0.:
                    yshift.append(max(depthset[repID][0] * 1.3, 1. / SNRset[repID][0] * 10.))
                else:
                    yshift.append(1. / SNRset[repID][0] * 10.)

                # if n == num - 1:
                #     yupp = yshift + 5. / SNRset[n] + 1.

                # spylist += [onespec, spy, spy_interp, spy_interp + error, spy_interp - error]
                # spxlist += [spx, spx, spx, spx, spx]
                # yshift_factor += [n, n, n, n, n]
                spylist += [onespec, spy, spy_interp]
                spxlist += [spx, spx, spx]
                yshift_factor += [n, n, n]
                lammin.append(numpy.amin(spx))
                lammax.append(numpy.amax(spx))
                intlist.append(intrange)

            yshift_max = max(yshift)
            ylow = 1. - yshift_max
            yupp = 1. + yshift_max * len(objIDax)
            spyshifted = [spylist[k] + yshift_max * yshift_factor[k] for k in range(len(yshift_factor))]
            ylevel = [1. + yshift_max * n for n in range(len(objIDax))]

            for n in range(len(objIDax)):
                repID = repIDdict[objIDax[n]]  # repIDlist[n]
                # spx = spxlist[n * 5]
                spx = spxlist[n * 3]
                onespec = numpy.ones(spx.shape)
                # spy_interp = spylist[2 + 5 * n]
                spy_interp = spylist[2 + 3 * n]
                if EWset[repID][0] != 0. and pflist[n]:
                    ax[j].fill_between(spx[intlist[n]], onespec[intlist[n]] + yshift_max * n,
                                       spy_interp[intlist[n]] + yshift_max * n, facecolor="y", alpha=0.5)

            # spectrum_plot(ax[j], spxlist, spyshifted, [min(lammin), max(lammax)], [ylow, yupp],
            #               colors=colorlist, yaxis_label="Normalized flux + offset",
            #               lines=["--", "-", "-", "-", "-"], linew=[1., 1., 2., 1., 1.])
            spectrum_plot(ax[j], spxlist, spyshifted, [min(lammin), max(lammax)], [ylow, yupp],
                          colors=colorlist, yaxis_label="Normalized flux + offset",
                          lines=["--", "-", "-"], linew=[1., 1., 2.])
            spectrum_plot(telax[j], [telwav_air], [telflux], [min(lammin), max(lammax)], [0., 1.1],
                          yaxis_label="Transmittance", xaxis_label="Wavelength ($\AA$)")
            ax[j].plot([wav_air[i], wav_air[i]], [ylow, yupp], "lightskyblue", linestyle="--")

            for n in range(len(objIDax)):
                repID = repIDdict[objIDax[n]]  # repIDlist[n]
                dbID = objectIDdb == objIDax[n]
                if pflist[n]:
                    if EWset[repID][0] != 0.:
                        righttext = "  %s\n" % registeredname[dbID][0] + "   E(B-V)=%.2f\n" % ebv[dbID][
                            0] + "   Sp: %s\n" % \
                                    sptype[dbID][0] + "   SNR=%.1f\n" % SNRset[repID][0] + "   EW=%.2f" % EWset[repID][
                                        0] + r" $\pm$" + "%.2f mA\n" % EWerrset[repID][0] + "   comment: %s" % \
                                    com[repID][0]
                    #     righttext = " %s\n" % mID[repID][0] + " EW=%.2f" % EWset[repID][0] + r" $\pm$" + "%.2f mA\n" % \
                    #                 EWerrset[repID] + " FWHM=%.2f A\n" % FWHMset[repID][0] + " depth=%.2f\n" % \
                    #                 depthset[repID][0] + " center=%.3f A\n" % cl[repID] + " v_hel = %.2f km/s\n" % \
                    #                 hv[repID][0] + " SNR=%.1f\n" % SNRset[repID] + " range=%.3f-%.3f A\n" % (
                    #                     ints[repID][0], inte[repID][0]) + " comment: %s" % com[repID][0]
                    else:
                        righttext = "  %s\n" % registeredname[dbID][0] + "   E(B-V)=%.2f\n" % ebv[dbID][
                            0] + "   Sp: %s\n" % \
                                    sptype[dbID][0] + "   SNR=%.1f\n" % SNRset[repID][
                                        0] + "   EW" + r"$<$" + "%.2f mA\n" % \
                                    EWerrset[repID][0] + "   comment: %s" % com[repID][0]
                    #     righttext = " %s\n" % mID[repID][0] + " EW" + r"$<$" + "%.2f mA\n" % EWerrset[
                    #         repID] + " center=%.3f A\n" % cl[repID] + " v_hel = %.2f km/s\n" % hv[repID][
                    #                     0] + " SNR=%.1f\n" % SNRset[repID] + " range=%.3f-%.3f A\n" % (
                    #                     ints[repID][0], inte[repID][0]) + " comment: %s" % com[repID][0]
                else:
                    righttext = "  %s\n" % registeredname[dbID][0] + "   E(B-V)=%.2f\n" % ebv[dbID][0] + "   Sp: %s\n" % \
                                sptype[dbID][0] + "   SNR=%.1f\n" % SNRset[repID][0] + "   comment: %s" % com[repID][0]
                ax[j].text(max(lammax), ylevel[n], righttext, va="top", fontsize=5)

        plt.suptitle("DIB %.1f (ID=%d), Page 1/2" % (wav_air[i], i))
        plt.savefig(pp, format="pdf")
        plt.clf()
        plt.close(fig)

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
            if pf[repIDdict[o]][0] == 1:
                objIDpf.append(o)
                ebvpf.append(ebv[objectIDdb == o][0])
                ewlist.append(EWset[repIDdict[o]][0])
                ewerrlist.append(EWerrset[repIDdict[o]][0])
                ewdict[o] = EWset[repIDdict[o]][0]
                ewerrdict[o] = EWerrset[repIDdict[o]][0]
                if EWset[repIDdict[o]][0] > 0.:
                    objIDpfdet.append(o)
                    fwhmlist.append(FWHMset[repIDdict[o]][0])
                    heliovdict[o] = hv[repIDdict[o]][0]
        objIDpf = numpy.array(objIDpf)
        ebvpf = numpy.array(ebvpf)
        ewlist = numpy.array(ewlist)
        ewerrlist = numpy.array(ewerrlist)

        column = ["ID", "Object", "E(B-V)", r"EW (m$\AA$)", r"FWHM ($\AA$)", r"Hel. v (km/s)", "depth"]
        rowtexts = []
        for o in objectID:
            if pf[repIDdict[o]][0] == 1:
                if EWset[repIDdict[o]][0] > 0.:
                    rowtexts.append([str(o), registeredname[objectIDdb == o][0], str(ebv[objectIDdb == o][0]),
                                     "%.1f" % EWset[repIDdict[o]][0] + r"$\pm$" + "%.1f" % EWerrset[repIDdict[o]][0],
                                     "%.2f" % FWHMset[repIDdict[o]][0], "%.2f" % hv[repIDdict[o]][0]])
                else:
                    rowtexts.append([str(o), registeredname[objectIDdb == o][0], str(ebv[objectIDdb == o][0]),
                                     r"$<$" + "%.1f" % EWerrset[repIDdict[o]][0],
                                     "-", "-"])

            else:
                rowtexts.append([str(o), registeredname[objectIDdb == o][0], str(ebv[objectIDdb == o][0]),
                                 "-", "-", "-"])

        dibtb = tableax.table(cellText=rowtexts, colLabels=column, bbox=[0, 0, 1, 1])
        dibtb.auto_set_font_size(False)
        dibtb.set_fontsize(9)
        dibtb.auto_set_column_width(col=list(range(len(column))))
        tableax.axis("off")

        for m in range(6):
            dibtb[0, m].set_facecolor('#363636')
            dibtb[0, m].set_text_props(color='w')

        r_ebv, p_ebv = ebv_correlation(ebvax, ebvpf, ewlist, ewerrlist, objIDpf, wav_air[i], stylefile=stylefile)
        r_dib, p_dib = DIBcorrelation(dibax, ew13175, ewerr13175, objectid13175, DIBwav13175, ewdict, ewerrdict,
                                      objIDpf, wav_air[i],
                                      stylefile=stylefile)
        fwhmax.hist(fwhmlist)
        fwhmax.set_xlabel(r'FWHM ($\AA$)')
        Vcorrelation(vax, heliov10780, objectid10780, DIBwav10780, heliovdict, objIDpfdet, wav_air[i],
                     stylefile=stylefile)

        descax.set_xlim(-1, 1)
        descax.set_ylim(-1, 1)
        descax.text(0, 0.4, r"$\lambda _{air}$ = " + "%.2f" % wav_air[i] + r"$\AA$, " + "FWHM = " + "%.2f" % DIBfwhm[
            i] + r"$\AA$" + "\nreference: %s, " % reference[i] + "category: %s\n" % category[i] + "comment: %s" %
                    DIBcomment[i], ha="center", va="center", fontsize=11)
        descax.text(0, -0.4,
                    "Detection rate: %d/%d objects (%.1f percent)\nFWHM (average): %.2fA\nCorrelation coeff. w/ E(B-V): %.2f (p=%.1e)\nCorrelation coeff. w/ DIB %d: %.2f (p=%.1e)" % (
                        numpy.sum(ewlist > 0.), len(objectID), (numpy.sum(ewlist > 0.) / len(objectID) * 100.),
                        numpy.average(fwhmlist), r_ebv, p_ebv, DIBwav13175, r_dib, p_dib), ha="center",
                    va="center",
                    fontsize=11)
        descax.xaxis.set_visible(False)
        descax.yaxis.set_visible(False)

        plt.suptitle("DIB %.1f (ID=%d), Page 2/2" % (wav_air[i], i))
        plt.savefig(pp, format="pdf")
        plt.clf()
        plt.close(fig)

    pp.close()

    conn.close()
