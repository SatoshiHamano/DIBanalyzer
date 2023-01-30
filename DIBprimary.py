# -*- coding:utf-8 -*-

# python DIBprimary.py temporaly_files/CygOB2cluster_DIBprimary.pdf -o 12 14 13 10 15 11 9

from open_mysql_project import openproject
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from spectra_plotter import MultiSpecPlotter
import numpy
import sys, os
from add_lineDIB_mysql import GetDIBdict
from spectra_plotter import spectrum_plot
from DIBanalysis import openDIBresult
from Spec1Dtools import openspecfits

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outputpdf", type=str, help="Output pdf")
    parser.add_argument("-o", "--objectID", type=int, help="object ID", nargs='*')
    # parser.add_argument("-d", "--dibid", type=int, default=0, help="DIB ID")

    args = parser.parse_args()
    objectID = args.objectID
    # dibid = args.dibid
    outputpdf = args.outputpdf
    objectID.reverse()

    objstr = ''
    for i in objectID:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")

    conn, cur = openproject()

    [combineID, measurementID, DIBID, echelleorder, DIBspecpath, primaryflag, autonormalizeflag, automeasurementflag,
     EW, EWerr, centerlam_air, helio_velocity, FWHM, SNR, integration_start, integration_end, comment, depth] = [numpy.array([]) for
                                                                                                                 i in
                                                                                                                 range(
                                                                                                                     18)]


    for o in objectID:
        cur.execute("SELECT x.combineID, x.measurementID, x.DIBID, x.echelleorder, x.DIBspecpath, "
                    "x.primaryflag, x.autonormalizeflag, x.automeasurementflag, x.EW, x.EWerr, x.centerlam_air, "
                    "x.helio_velocity, x.FWHM, x.SNR, x.integration_start, x.integration_end, x.comment, x.depth "
                    "from DIBmeasurement as x join combinesummary as y using(combineID) join DIBlist as z using(DIBID) "
                    "where y.objectID=%d and x.primaryflag=1 and z.category != 'fake';" % o)
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

    DIBinfo = GetDIBdict()
    [DIBIDlist, wav_air, reference, category, fwhm, wavenumber, DIBcomment] = DIBinfo

    DIBset = list(set(DIBID))
    wav_set = [wav_air[i] for i in DIBset]
    DIBzip = zip(wav_set, DIBset)
    DIBzip_sorted = sorted(DIBzip)
    wav_set, DIBset = zip(*DIBzip_sorted)

    pp = PdfPages(outputpdf)

    figw = 0.63
    orig_w = 0.08

    for i in DIBset:
        req = DIBID == i
        num = numpy.sum(req)
        orig_h = 0.5 / (3 * num)
        figh = (3 * num - 1.0) / (3 * num)

        fig = plt.figure(figsize=(12, 3 * num))
        ax1 = plt.axes([orig_w, orig_h, figw, figh])

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

        spxlist = []
        spylist = []
        lammin = []
        lammax = []
        yshift = []
        yshift_factor = []
        intlist = []

        for n in range(num):
            DIBdir = os.path.dirname(sp[n]) + "/"
            resultnpz = sp[n].rstrip("fits") + "npz"
            results = openDIBresult(resultnpz)
            [mask, spx, spy_interp, sn, error, error_interp, error_cont] = results
            intrange = (spx > ints[n] - 1.e-6) & (spx < inte[n] + 1.e-6)
            onespec = numpy.ones(spx.shape)
            _, spy, _, _, _ = openspecfits(sp[n])
            # if n == 0:
            #     if depthset[n] != 0.:
            #         ylow = 1. - max(depthset[n] * 1.3, 1. / SNRset[n] * 10.)
            #     else:
            #         ylow = 1. - 10. / SNRset[n]
            # else:
            if depthset[n] != 0.:
                yshift.append(max(depthset[n] * 1.3, 1. / SNRset[n] * 10.))
            else:
                yshift.append(1. / SNRset[n] * 10.)

            # if n == num - 1:
            #     yupp = yshift + 5. / SNRset[n] + 1.

            spylist += [onespec, spy, spy_interp, spy_interp + error, spy_interp - error]
            spxlist += [spx, spx, spx, spx, spx]
            yshift_factor += [n, n, n, n, n]
            lammin.append(numpy.amin(spx))
            lammax.append(numpy.amax(spx))
            intlist.append(intrange)

        yshift_max = max(yshift)
        ylow = 1. - yshift_max
        yupp = 1. + yshift_max * num
        spyshifted = [spylist[k] + yshift_max * yshift_factor[k] for k in range(len(yshift_factor))]
        ylevel = [1. + yshift_max * n for n in range(num)]

        for n in range(num):
            spx = spxlist[n * 5]
            onespec = numpy.ones(spx.shape)
            spy_interp = spylist[2 + 5 * n]
            if EWset[n] != 0.:
                ax1.fill_between(spx[intlist[n]], onespec[intlist[n]] + yshift_max * n,
                                 spy_interp[intlist[n]] + yshift_max * n,
                                 facecolor="y",
                                 alpha=0.5)

        spectrum_plot(ax1, spxlist, spyshifted, [min(lammin), max(lammax)], [ylow, yupp],
                      colors=["0.5", "0.5", "b", "k", "k"], yaxis_label="Normalized flux + offset",
                      lines=["--", "-", "-", "-", "-"], linew=[1., 1., 2., 1., 1.], xaxis_label="Wavelength ($\AA$)")

        for n in range(num):
            if EWset[n] != 0.:
                righttext = " %s\n" % mID[n] + " EW=%.2f" % EWset[n] + r" $\pm$" + "%.2f mA\n" % EWerrset[
                    n] + " FWHM=%.2f A\n" % FWHMset[n] + " depth=%.2f\n" % depthset[n] + " center=%.3f A\n" % cl[
                                n] + " v_hel = %.2f km/s\n" % hv[n] + " SNR=%.1f\n" % SNRset[
                                n] + " range=%.3f-%.3f A\n" % (ints[n], inte[n]) + " comment: %s" % comment[n]
            else:
                righttext = " %s\n" % mID[n] + " EW" + r"$<$" + "%.2f mA\n" % EWerrset[
                    n] + " center=%.3f A\n" % cl[
                                n] + " v_hel = %.2f km/s\n" % hv[n] + " SNR=%.1f\n" % SNRset[
                                n] + " range=%.3f-%.3f A\n" % (ints[n], inte[n]) + " comment: %s" % comment[n]

            if pf[n] == 1:
                plt.text(max(lammax), ylevel[n], righttext, va="top", fontweight="bold")
            else:
                plt.text(max(lammax), ylevel[n], righttext, va="top")

        plt.title("DIB %.1f" % wav_air[i])
        plt.savefig(pp, format="pdf")
        plt.clf()
        plt.close(fig)

    pp.close()

    conn.close()
