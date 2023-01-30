# -*- coding:utf-8 -*-

import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from open_mysql_project import openproject
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import math
from add_lineDIB_mysql import GetDIBdict, GetDIBLine
from spectra_plotter import spectrum_plot, linelist_marker
from DIBanalysis import openDIBresult
from Spec1Dtools import openspecfits
import scipy.constants
from DIBsummary_forWeakDIB import telStatus

if __name__ == '__main__':
    c = scipy.constants.c * 1.e-3

    h_imnum = 3
    v_max = 3
    v_baseratio = 1-0.13*(6-v_max) - 0.02*(6-v_max)
    h_baseratio = 1-0.12*(6-h_imnum) - 0.03*(6-h_imnum)
    h_size = 0.12 / h_baseratio
    h_interval = 0.03 / h_baseratio
    h_orig = 0.1 / h_baseratio

    outputpdf = sys.argv[1]

    objectID = [149, 66]
    DIBIDs = [33,34,35,37,38]
    objstr = ''
    for i in objectID:
        objstr += "%d," % i
    objstr = objstr.rstrip(",")
    DIBIDstr = ''
    for i in DIBIDs:
        DIBIDstr += "%d," % i
    DIBIDstr = DIBIDstr.rstrip(",")

    conn, cur = openproject()


    [combineID, measurementID, DIBID, echelleorder, DIBspecpath, primaryflag, autonormalizeflag, automeasurementflag,
     EW, EWerr, centerlam_air, helio_velocity, FWHM, SNR, integration_start, integration_end, comment, depth,
     objIDref] = [np.array([]) for i in range(19)]

    cur.execute(
        "select objectID,E_BV,sptype,registeredname from object join objectdict using(objectID) "
        "where objectID IN (%s) and priority=1 order by E_BV;" % objstr)
    rows = cur.fetchall()
    objectIDdb = np.array([int(i[0]) for i in rows])
    ebv = np.array([float(i[1]) for i in rows])
    sptype = np.array([i[2] for i in rows])
    registeredname = np.array([i[3] for i in rows])

    cur.execute("SELECT DIBID, wavelength_air, category from DIBlist where category NOT IN ('fake', 'Cs I', 'weak') "
                "and wavelength_air between 9000. and 13500. and DIBID IN ({}) order by wavelength_air;".format(DIBIDstr))
    rows = cur.fetchall()
    DIBIDall = np.array([i[0] for i in rows])
    wavDIB = np.array([i[1] for i in rows])
    category = np.array([i[2] for i in rows])
    DIBtags = []
    for i in range(len(DIBIDall)):
        if category[i] == "candidate":
            DIBtags.append("DIB?\n {:.0f}".format(math.floor(round(wavDIB[i],1))))
        elif category[i] != "candidate":
            DIBtags.append("DIB\n {:.0f}".format(math.floor(round(wavDIB[i],1))))
    DIBtags = np.array(DIBtags)

    for o in objectID:
        cur.execute("SELECT x.combineID, x.measurementID, x.DIBID, x.echelleorder, x.DIBspecpath, "
                    "x.primaryflag, x.autonormalizeflag, x.automeasurementflag, x.EW, x.EWerr, x.centerlam_air, "
                    "x.helio_velocity, x.FWHM, x.SNR, x.integration_start, x.integration_end, x.comment, x.depth "
                    "from DIBmeasurement as x join combinesummary as y using(combineID) join DIBlist as z using(DIBID) "
                    "where y.objectID=%d and z.category NOT IN ('fake', 'Cs I', 'weak', 'candidate') and z.DIBID IN (%s);" % (o, DIBIDstr))
        rows = cur.fetchall()

        if rows == []:
            print("DIBs are not measured for Object ID = %d" % objectID)
            sys.exit()
        else:
            combineID = np.append(combineID, np.array([i[0] for i in rows]))
            measurementID = np.append(measurementID, np.array([i[1] for i in rows]))
            DIBID = np.append(DIBID, np.array([i[2] for i in rows]))
            echelleorder = np.append(echelleorder, np.array([i[3] for i in rows]))
            DIBspecpath = np.append(DIBspecpath, np.array([i[4] for i in rows]))
            primaryflag = np.append(primaryflag, np.array([i[5] for i in rows]))
            autonormalizeflag = np.append(autonormalizeflag, np.array([i[6] for i in rows]))
            automeasurementflag = np.append(automeasurementflag, np.array([i[7] for i in rows]))
            EW = np.append(EW, np.array([i[8] for i in rows]))
            EWerr = np.append(EWerr, np.array([i[9] for i in rows]))
            centerlam_air = np.append(centerlam_air, np.array([i[10] for i in rows]))
            helio_velocity = np.append(helio_velocity, np.array([i[11] for i in rows]))
            FWHM = np.append(FWHM, np.array([i[12] for i in rows]))
            SNR = np.append(SNR, np.array([i[13] for i in rows]))
            integration_start = np.append(integration_start, np.array([i[14] for i in rows]))
            integration_end = np.append(integration_end, np.array([i[15] for i in rows]))
            comment = np.append(comment, np.array([i[16] for i in rows]))
            depth = np.append(depth, np.array([i[17] for i in rows]))
            objIDref = np.append(objIDref, np.array([o for i in rows]))



    DIBinfo = GetDIBdict()
    [DIBIDlist, wav_air, reference, category, DIBfwhm, wavenumber, DIBcomment] = DIBinfo

    DIBset = list(set(DIBID))
    wav_set = [wav_air[i] for i in DIBset]
    DIBzip = zip(wav_set, DIBset)
    DIBzip_sorted = sorted(DIBzip)
    wav_set, DIBset = zip(*DIBzip_sorted)

    dibnum = len(DIBIDs)

    pdfcounter = 1
    dibcounter = 0
    firstflag = True

    print("Number of DIBs: ", dibnum)
    pdfnumflag = True if dibnum < h_imnum * v_max else False

    while dibcounter < dibnum - 1:
        if pdfnumflag:
            pp = PdfPages(outputpdf)
        else:
            pp = PdfPages(outputpdf.rstrip("pdf").rstrip(".") + "_{}.pdf".format(pdfcounter))

        v_imnum = min(v_max, (dibnum - dibcounter + h_imnum - 1)//h_imnum)
        print("Yeah", (dibnum - dibcounter), (dibnum - dibcounter + h_imnum - 1)//h_imnum)
        v_ratio = 1-0.13*(v_max-v_imnum) - 0.02*(v_max-v_imnum)
        v_size = 0.13 / v_ratio / v_baseratio#/v_imnum*v_max
        v_interval = 0.02 / v_ratio / v_baseratio#/v_imnum*v_max
        v_orig = 0.08 / v_ratio / v_baseratio#/v_imnum*v_max
        plt.figure(figsize=(18 * h_baseratio, 22 * v_ratio * v_baseratio))

        print(pdfcounter, v_imnum)

        for j in range(v_imnum):
            for i in range(h_imnum):
                print(DIBIDs[dibcounter])
                ax = plt.axes([h_orig + (h_size + h_interval) * i, v_orig + (v_size + v_interval) * (v_imnum - j - 1), h_size, v_size])
                req = DIBID == DIBIDs[dibcounter]
                num = np.sum(req)

                cID = combineID[req]
                tstatus = np.array([telStatus(cur, c) for c in cID])
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
                spxlist = []
                spvlist = []
                spylist = []
                lammin = []
                lammax = []
                velmin = []
                velmax = []
                yshift = []
                yshift_factor = []
                intlist = []
                pflist = []
                colorlist = []

                for n in range(len(objectID)):
                    if 1 in pf[oIDref == objectID[n]]:
                        repID = (oIDref == objectID[n]) & (pf == 1)
                        pflist.append(True)
                        colorlist += ["0.5", "0.5", "k"]
                    else:
                        print(SNRset)
                        SNRmax = np.amax(SNRset[(oIDref == objectID[n])])# & (tstatus == "AT")])
                        repID = (oIDref == objectID[n]) & (tstatus == "AT") & (SNRset == SNRmax)
                        pflist.append(False)
                        colorlist += ["0.25", "0.25", "0.5"]

                    repIDdict[objectID[n]] = repID

                    DIBdir = os.path.dirname(sp[repID][0]) + "/"
                    resultnpz = sp[repID][0].rstrip("fits") + "npz"
                    results = openDIBresult(resultnpz)
                    [mask, spx, spy_interp, sn, error, error_interp, error_cont] = results
                    intrange = (spx > ints[repID][0] - 1.e-6) & (spx < inte[repID][0] + 1.e-6)
                    onespec = np.ones(spx.shape)
                    _, spy, _, _, _ = openspecfits(sp[repID][0])
                    if depthset[repID][0] != 0.:
                        yshift.append(max(depthset[repID][0] * 1.5, 1. / SNRset[repID][0] * 10.))
                    else:
                        yshift.append(1. / SNRset[repID][0] * 10.)

                    spylist += [onespec, spy, spy_interp]
                    spxlist += [spx, spx, spx]
                    spv = (spx - wav_air[DIBset[dibcounter]]) / wav_air[DIBset[dibcounter]] * c
                    spvlist += [spv, spv, spv]
                    yshift_factor += [n, n, n]
                    lammin.append(np.amin(spx))
                    lammax.append(np.amax(spx))
                    velmin.append(np.amin(spv))
                    velmax.append(np.amax(spv))
                    intlist.append(intrange)

                yshift_max = np.max(yshift)
                ylow = 1. - yshift_max
                yupp = 1. + yshift_max * (len(objectID) + 1)
                spyshifted = [spylist[k] + yshift_max * yshift_factor[k] for k in range(len(yshift_factor))]
                ylevel = [1. + yshift_max * n for n in range(len(objectID))]

                for n in range(len(objectID)):
                    repID = repIDdict[objectID[n]]  # repIDlist[n]
                    spx = spxlist[n * 3]
                    spv = spvlist[n * 3]
                    onespec = np.ones(spx.shape)
                    spy_interp = spylist[2 + 3 * n]
                    if EWset[repID][0] != 0. and pflist[n]:
                        ax.fill_between(spv[intlist[n]], onespec[intlist[n]] + yshift_max * n,
                                           spy_interp[intlist[n]] + yshift_max * n, facecolor="y", alpha=0.5)

                wavDIBclip = wavDIB[DIBIDall != DIBset[dibcounter]]
                DIBtagsclip = DIBtags[DIBIDall != DIBset[dibcounter]]

                velDIB = (wavDIBclip - wav_air[DIBset[dibcounter]]) / wav_air[DIBset[dibcounter]] * c
                spectrum_plot(ax, spvlist, spyshifted, [-300, 300], [ylow, yupp],
                              colors=colorlist, yaxis_label="",
                              lines=["--", "-", "-"], linew=[1., 1., 2.])
                linelist_marker(ax, velDIB, DIBtagsclip, 0., -300., 300., 1. + yshift_max * (len(objectID) - 0.8),
                                1. + yshift_max * (len(objectID) - 0.5), 1. + yshift_max * (len(objectID) - 0.5), linew=1.5)

                ax.plot([0., 0.], [ylow, yupp], "lightskyblue", linestyle="--")
                ax.text(-290, yupp - (yupp - ylow)*0.02, r"$\lambda${:.0f}".format(math.floor(round(wav_air[DIBset[dibcounter]],1))), ha='left', va='top', fontsize=20)
                if firstflag:
                    for n in range(len(objectID)):
                        req = objectIDdb == objectID[n]
                        print(objectID[n],req)
                        print(registeredname, ebv, sptype)
                        ax.text(-290, 1. + yshift_max * 0.2 +  yshift_max * yshift_factor[n * 3] * 1.05, '{} \n(E(B-V)={})'.format(registeredname[req][0], ebv[req][0]), fontsize=8, ha='left', va='bottom')
                        # ax.text(-290, 1. + yshift_max * yshift_factor[n * 3] * 1.05, '{} \n(E(B-V)={}, {})'.format(registeredname[req][0], ebv[req][0], sptype[req][0].split()[0]), fontsize=8, ha='left', va='bottom')
                    firstflag = False

                dibcounter += 1

                if dibcounter > dibnum - 1: break
            if dibcounter > dibnum - 1: break

        plt.gcf().text(0.55, 0.03, "Velocity (km s$^{-1}$)", ha="center", fontsize=30)
        plt.gcf().text(0.03, 0.55, "Normalized flux", va="center", fontsize=30, rotation=90)

        plt.savefig(pp, format="pdf")
        plt.clf()
        firstflag = True

        pp.close()
        pdfcounter += 1

    conn.close()
