import sys, os, datetime
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == '__main__':
    rf = open(sys.argv[1], "r")
    rl = rf.readlines()
    rf.close()

    [pipelineID, order, frame, totalexp, totalSNR, mode, pipelinever, linewave, centerlams, centerpix, lambdac,
     linewave0, centerlams0, centerpix0, lambac0] = [np.array([i.split()[j] for i in rl]) for j in range(15)]

    years = np.array([int(p[0:4]) for p in pipelineID])
    months = np.array([int(p[5:7]) for p in pipelineID])
    days = np.array([int(p[8:10]) for p in pipelineID])

    order = order.astype(np.int32)
    totalexp = totalexp.astype(np.float64)
    totalSNR = totalSNR.astype(np.float64)
    # wavelength corrected
    linewave = linewave.astype(np.float64) #reference wavelength, true wavelength
    centerlams = centerlams.astype(np.float64) #shift from reference wavelength, should be 0 if the wavelength is correct.
    centerpix = centerpix.astype(np.float64) #corresponding pixel coordinate
    lambdac = lambdac.astype(np.float64) #no meaning?
    # pipeline reduced
    linewave0 = linewave0.astype(np.float64) #reference wavelength, true wavelength
    centerlams0 = centerlams0.astype(np.float64) #shift from reference wavelength, should be 0 if the wavelength is correct.
    centerpix0 = centerpix0.astype(np.float64) #corresponding pixel coordinate
    lambac0 = lambac0.astype(np.float64) #no meaning?

    pIDset = list(set(list(pipelineID)))
    print(len(pIDset))

    req_w = np.logical_not(np.logical_and(linewave > 10000., linewave < 10500.))

    pp = PdfPages(sys.argv[2])
    plt.figure()

    pipelineData = []
    waveshiftData = []
    pipelineDataS = []
    waveshiftDataS = []
    dates = []

    sig = 3.
    ite = 5

    for i in range(len(pIDset)):
        req_p = pipelineID == pIDset[i]
        frameset = list(set(list(frame[req_p])))
        for f in frameset:
            req_f = f == frame
            req = np.logical_and(np.logical_and(req_f, req_p), req_w)

            pipeave = np.average(centerlams0[req])
            wsave = np.average(centerlams[req])
            pipestd = np.std(centerlams0[req])
            wsstd = np.std(centerlams[req])

            for _ in range(ite):
                pipeclip = np.absolute(centerlams0[req] - pipeave) < pipestd * sig
                pipeave = np.average(centerlams0[req][pipeclip])
                pipestd = np.std(centerlams0[req][pipeclip])

            for _ in range(ite):
                wsclip = np.absolute(centerlams[req] - wsave) < wsstd * sig
                wsave = np.average(centerlams[req][wsclip])
                wsstd = np.std(centerlams[req][wsclip])

            # plt.scatter(linewave[req][wsclip], centerlams[req][wsclip], marker="x", label="Wavelength corrected spectra")
            # plt.scatter(linewave0[req][pipeclip], centerlams0[req][pipeclip], label="Pipeline reduced spectra")
            # plt.legend()
            # plt.title("{} ({})".format(pIDset[i], f))
            # plt.savefig(pp, format="pdf")
            # plt.clf()

            pipelineData.append(pipeave)
            waveshiftData.append(wsave)
            pipelineDataS.append(pipestd)
            waveshiftDataS.append(wsstd)
            dates.append(np.average(years[req] + (months[req]-1)/12. + days[req]/30./12.))

    pipelineData = np.array(pipelineData)
    waveshiftData = np.array(waveshiftData)
    pipelineDataS = np.zeros(pipelineData.shape)#np.array(pipelineDataS)
    waveshiftDataS = np.zeros(waveshiftData.shape)#np.array(waveshiftDataS)
    dates = np.array(dates)

    print(dates.size)

    plt.errorbar(dates, waveshiftData, yerr=waveshiftDataS, label="Wavelength corrected spectra", linestyle="None", marker=".")
    plt.errorbar(dates, pipelineData, yerr=pipelineDataS, label="Pipeline reduced spectra", linestyle="None", marker="x")
    plt.plot([min(dates), max(dates)], [0.2, 0.2], "k--", label="1 pix")
    plt.plot([min(dates), max(dates)], [0., 0.], "0.5")
    plt.plot([min(dates), max(dates)], [-0.2, -0.2], "k--")
    plt.legend(loc="lower left")
    plt.xlabel("Year")
    plt.ylabel("Diff. from true wavelength (angstrom)")
    plt.savefig(pp, format="pdf")
    plt.clf()

    plt.errorbar(dates, waveshiftData, yerr=waveshiftDataS, label="Wavelength corrected spectra", linestyle="None", marker=".")
    plt.errorbar(dates, pipelineData, yerr=pipelineDataS, label="Pipeline reduced spectra", linestyle="None", marker="x")
    plt.plot([min(dates), max(dates)], [0.2, 0.2], "k--", label="1 pix")
    plt.plot([min(dates), max(dates)], [0., 0.], "0.5")
    plt.plot([min(dates), max(dates)], [-0.2, -0.2], "k--")
    plt.legend(loc="lower right")
    plt.xlabel("Year")
    plt.ylabel("Diff. from true wavelength (angstrom)")
    plt.ylim(-0.25, 0.25)
    plt.savefig(pp, format="pdf")
    plt.clf()

    plt.figure(figsize=(8,8))
    plt.errorbar(waveshiftData, pipelineData, xerr=waveshiftDataS, yerr=pipelineDataS, linestyle="None", marker=".")
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5,0.5)
    plt.plot([-0.5, 0.5], [0.,0.], "0.5")
    plt.plot([0.,0.], [-0.5,0.5], "0.5")
    plt.xlabel("Wavelength corrected spectra")
    plt.ylabel("Pipeline reduced spectra")
    plt.savefig(pp, format="pdf")


    plt.figure(figsize=(8,8))
    plt.errorbar(waveshiftData[dates>2017], pipelineData[dates>2017], xerr=waveshiftDataS[dates>2017], yerr=pipelineDataS[dates>2017], linestyle="None", marker=".")
    plt.xlim(-0.05, 0.05)
    plt.ylim(-0.05,0.05)
    plt.plot([-0.05, 0.05], [0.,0.], "0.5")
    plt.plot([0.,0.], [-0.05,0.05], "0.5")
    plt.xlabel("Wavelength corrected spectra")
    plt.ylabel("Pipeline reduced spectra")
    plt.title("NTT17c")
    plt.savefig(pp, format="pdf")

    print(np.average(np.absolute(waveshiftData[dates>2017])))
    print(np.average(np.absolute(pipelineData[dates>2017])))

    pp.close()
