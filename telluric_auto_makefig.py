# -*- coding:utf-8 -*-

import copy
from add_lineDIB_mysql import *
from obtain_filelist import *
from spectra_plotter import *
from open_mysql_project import openproject
from dopcor_vcorr import *

if __name__ == "__main__":
    conn, cur = openproject()

    boolstr = {True: "True", False: "False"}

    telluricID = sys.argv[1]
    if len(sys.argv) == 2:
        vacorair = "AIR"
    elif sys.argv[2] == "VAC":
        vacorair = "VAC"
    elif sys.argv[2] == "AIR":
        vacorair = "AIR"
    else:
        print(sys.argv[2] + " is neither vac or air.")
        vacorair = "AIR"

    print("Wavelength system: %s" % vacorair)

    cur.execute(
        "select y.pipelineID, y.FrameNum, x.telluricPath, y.pipelinever, y.path, x.autoflag, x.pipelineIDtel, x.advanced, y.obsdate from telluriccorrection as x join datareduction as y on x.pipelineIDobj=y.pipelineID where x.telluricID='%s';" %
        telluricID)
    rows = cur.fetchall()
    if rows == []:
        print(telluricID + " could not be found.")
        sys.exit()
    elif len(rows) == 1:
        ppid = rows[0][0]
        fnum = int(rows[0][1])
        telpath = rows[0][2]
        pver = rows[0][3]
        path = rows[0][4]
        autoflag = True if rows[0][5] == 1 else False
        ppidtel = rows[0][6]
        advflag = True if rows[0][7] == 1 else False
        obsdate = rows[0][8]
    else:
        print("Multiple records were found for " + telluricID)
        sys.exit()

    if not autoflag:
        print(telluricID + "is the manually reduced data.")
        sys.exit()

    v_helio = read_rvfile_ID(telluricID)

    pdffile = telpath + telluricID + ".pdf"
    if os.path.exists(pdffile):
        print("PDF was already made.")
        sys.exit()
    pp = PdfPages(pdffile)

    m, frameset, spf = obtainObjTelFileList(telpath, vacorair=vacorair, fluxornorm="flux")
    telfiles = obtainWscorTelluricFileList(path, ppidtel, advflag, vacorair=vacorair)
    colorsp = ["k" for i in spf]
    colorsp[-1] = "b"

    fig = plt.figure(figsize=(20, 5 + 2 * fnum))

    for i in range(len(m)):
        spxlist = []
        spylist = []
        lammin, lammax, ymin, ymax, yrange, ymed, ystd = [], [], [], [], [], [], []
        scale, scaleerr, shift, shifterr = [], [], [], []

        telx, tely, _, _, _ = openspecfits(telfiles[i])

        for j in range(len(frameset)):
            cur.execute(
                "select scale, scaleerr, shift, shifterr from telluricresult where telluricID='%s' and frame='%s' and echelleorder=%d;" %
                (telluricID, frameset[j], m[i]))
            rows = cur.fetchall()
            scale.append(float(rows[0][0]))
            scaleerr.append(float(rows[0][1]))
            shift.append(float(rows[0][2]))
            shifterr.append(float(rows[0][3]))

            spx, spy, _, _, _ = openspecfits(spf[frameset[j]][m[i]])
            # ymed1 = numpy.median(spy)
            # ystd1 = numpy.std(spy)
            # yclip = numpy.absolute(spy - ymed1) < plotsig * ystd1
            # # print("start", m[i], ymed1, ystd1, numpy.sum(yclip) / yclip.size)
            # for k in range(ite):
            #     ystdprev = copy.copy(ystd1)
            #     ymed1 = numpy.median(spy[yclip])
            #     ystd1 = numpy.std(spy[yclip])
            #     if ystdprev == ystd1:
            #         # print("end",k, m[i], ymed1, ystd1, numpy.sum(yclip) / yclip.size)
            #         break
            #     yclip = numpy.absolute(spy - ymed1) < plotsig * ystd1
            #     # if k == (ite - 1):
            #     #     print("end",k, m[i], ymed1, ystd1, numpy.sum(yclip) / yclip.size)
            #
            # ymed.append(ymed1)
            # ystd.append(ystd1 / ymed1)
            spylist.append(spy)
            spxlist.append(spx)
            # lammin.append(numpy.amin(spxlist[j]))
            # lammax.append(numpy.amax(spxlist[j]))
            # ymin.append(numpy.amin(spy / ymed1))
            # ymax.append(numpy.amax(spy / ymed1))
            # yrange.append(ymax[j] - ymin[j])

        ax1 = plt.axes([0.04, 0.2, 0.85, 0.7])
        texts = [" %s\n  (scale=%.3f$\pm$%.3f)\n  (shift=%.3f$\pm$%.3f)" % (
            frameset[j], scale[j], scaleerr[j], shift[j], shifterr[j]) for j in range(len(spf))]
        MultiSpecPlotter(spxlist, spylist, pp, colorsp, texts,
                         "%s (m=%d, auto=%s, advanced=%s)" % (telluricID, m[i], boolstr[autoflag], boolstr[advflag]),
                         fs=6., obsdates=obsdate, order=m[i], v_helio=v_helio, telx=telx, tely=tely, telflag=True)

        # lamshort = min(lammin)
        # lamlong = max(lammax)
        #
        # DIBreturn = GetDIBList(lamshort, lamlong, vacorair=vacorair)
        # STEreturn = GetLineList(lamshort, lamlong, sptype="OB", vacorair=vacorair)
        # AFreturn = GetAFList(lamshort, lamlong, obsdate, m[i], v_helio, vacorair=vacorair)
        #
        # ystdave = numpy.average(ystd)
        # ystep = ystdave * plotsig
        # ax1 = plt.axes([0.04, 0.2, 0.85, 0.7])
        # spectrum_plot(ax1, spxlist, spylist, [lamshort, lamlong],
        #               [1.0 - ystep, 1.0 + ystep * (fnum + 4)], colors=colorsp,
        #               yshift=[ystep * (fnum - k) for k in range(fnum + 1)],
        #               yaxis_label=r"Normalized flux + offset")
        # plt.title("%s (m=%d, auto=%s, advanced=%s)" % (telluricID, m[i], boolstr[autoflag], boolstr[advflag]))
        # for j in range(len(spf)):
        #     plt.text(lamlong, 1. + ystep * (fnum - j), " %s\n  (scale=%.3f$\pm$%.3f)\n  (shift=%.3f$\pm$%.3f)" % (
        #         frameset[j], scale[j], scaleerr[j], shift[j], shifterr[j]), color=colorsp[j])
        # if AFreturn != None:
        #     [AFstart, AFend, slitposition] = AFreturn
        #     ShadeAFregions(ax1, AFstart, AFend, slitposition, 1.0 - ystep, 1.0 + ystep * (fnum + 4))
        # if DIBreturn != None:
        #     [DIBID, wavDIB, refDIB, cateDIB, _, _, _] = DIBreturn
        #     DIBtags = numpy.array(["DIB%d" % (wavDIB[k]) for k in range(len(DIBID))])
        #     linelist_marker(ax1, wavDIB, DIBtags, 0., lamshort, lamlong, 1.0 + ystep * (fnum + 0.2),
        #                     1.0 + ystep * (fnum + 0.6), 1.0 + ystep * (fnum + 1.5))
        # if STEreturn != None:
        #     [StellarID, wavSte, atomSte, ionSte, sptypeSte] = STEreturn
        #     STEtags = numpy.array(["%s%s %.1f" % (atomSte[k], ionSte[k], wavSte[k]) for k in range(len(StellarID))])
        #     linelist_marker(ax1, wavSte, STEtags, 0., lamshort, lamlong, 1.0 + ystep * (fnum + 1.5),
        #                     1.0 + ystep * (fnum + 1.9), 1.0 + ystep * (fnum + 2.8), itemax=1000)
        # C2y = 1.0 + ystep * (fnum + 2.8)
        # for k in range(3):
        #     if lamshort < C2range[k][1] and lamlong > C2range[k][0]:
        #         C2marks(ax1, bandnames[k], 0., 0.,
        #                 ycoord=[C2y + 0.85 * ystep, C2y + 0.55 * ystep, C2y + 0.25 * ystep, C2y])
        #
        # ax2 = plt.axes([0.04, 0.1, 0.85, 0.1])
        # spectrum_plot(ax2, [telx], [tely], [lamshort, lamlong], [0., 1.2],
        #               xaxis_label=r"Wavelength %s ($\AA$)" % vacorair.lower(),
        #               yaxis_label="Transmittance")
        #
        # plt.savefig(pp, format="pdf")
        # plt.clf()

    pp.close()
    conn.close()
