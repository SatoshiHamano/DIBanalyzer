import sys, os, datetime
import numpy, math
from pyraf import iraf
import scipy.optimize
from Spec1Dtools import openspecfits
from open_mysql_project import openproject
from obtain_filelist import obtainWscorFileList
import glob

iraf.noao()
iraf.onedspec()


# Description:
#   This "Waveshift_main.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to measure the wavelength shift.
#
#   The change of the wavelength is measured by fitting Gaussian functions to the
#   telluric absorption lines, which are enough strong and not blended with other lines.
#   Because the change of the wavelength is currently known to linearly depend on the
#   wavelength, a linear function is fitted to the plot of wavelength - wavelength change.
#   During the fitting, the outlier, which would be originated from the failure of Gaussian
#   fitting, is clipped.
#
# Usage:
#
#   $ python Waveshift_fit.py <speclist> <output> <telluric linelist> <resolution>
#
#   speclist -- the list of fits files to be measured
#   output -- the result of the fitting
#   telluric linelist -- the list of telluric lines to be used to measure the wavelength shift
#   resolution -- the spectral resolving power (R= lambda / delta lambda) to be used for Gaussian fitting
#
# Output:
#
#   The function of wavelength shift as a function of wavelength.
#
# Updates:
#
#   (ver.2)
#
#   Updated by S.Hamano in 2016.05.09
#
#   A bug in the clipping procedure is fixed.
#   Clipping procedure becomes efficient a little bit.
#
##

def gaussianfunc(p, x, c, w):
    y = 1.0 - p[0] * numpy.exp(-((x - p[1] - c) / w) ** 2)
    return y


def residue(p, y, x, c, w):
    res = (y - gaussianfunc(p, x, c, w))
    return (res)


def linearfunc(p, x):
    y = p[0] * x + p[1]
    return y


def residue2(p, y, x):
    res = (y - linearfunc(p, x))
    return (res)


# def openspecfits(fitsfile):
#     spfits = pyfits.open(fitsfile)
#     splength = spfits[0].header["NAXIS1"]
#     spdata = spfits[0].data
#
#     rcrval1 = float(spfits[0].header["CRVAL1"])
#     rcdelt1 = float(spfits[0].header["CDELT1"])
#     rcrpix1 = float(spfits[0].header["CRPIX1"])
#
#     lamx = numpy.array([rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)])
#     spfits.close()
#
#     return lamx, spdata, rcrval1, rcdelt1, rcrpix1

def Waveshift(speclist, orderlist, linecenterlist, resolution, center_wave, width=15., sigma=3., iterate=3):
    peak = []
    centerlams = []
    linewave = []
    centerpix = []
    inc = []
    order = []
    lambdac = []
    offset = []
    for i in range(len(speclist)):
        spf = speclist[i].split()[0]
        lamx, spdata, rcrval1, rcdelt1, rcrpix1 = openspecfits(spf)
        minlam = numpy.amin(lamx)
        maxlam = numpy.amax(lamx)

        for j in range(len(linecenterlist)):
            if linecenterlist[j] - width > minlam and linecenterlist[j] + width < maxlam:
                iraf.scopy(spf, "tmp%s" % spf.split("/")[-1], w1="%.2f" % (linecenterlist[j] - width),
                           w2="%.2f" % (linecenterlist[j] + width))
                iraf.continuum("tmp%s" % spf.split("/")[-1], "normtmp%s" % spf.split("/")[-1], interactive="no",
                               low_rej=2., high_rej=3., grow=1., order=5, func="legendre")

                normlamx, normdata, rcrval2, rcdelt2, rcrpix2 = openspecfits("./normtmp%s" % spf.split("/")[-1])

                os.remove("./normtmp%s" % spf.split("/")[-1])
                os.remove("./tmp%s" % spf.split("/")[-1])

                gwidth = linecenterlist[j] / resolution / 2.35 * 1.414

                p0 = [0.5, 0.1]
                param_output = scipy.optimize.leastsq(residue, p0, args=(normdata, normlamx, linecenterlist[j], gwidth),
                                                      full_output=True)

                peak.append(param_output[0][0])
                centerlams.append(param_output[0][1])
                linewave.append(linecenterlist[j])
                centerpix.append((linecenterlist[j] - rcrval1) / rcdelt1 + rcrpix1 - 1.)
                lambdac.append(center_wave[i])
                order.append(orderlist[i])

    return peak, centerlams, linewave, centerpix, lambdac, order

    # for i in range(len(centerlams)):
    #     wshiftfile.write("%.5f\t%.5f\t%.5f\t%.0f\n" % (linewave[i], centerlams[i], centerpix[i], lambdac[i]))


if __name__ == "__main__":

    filename = sys.argv[1:]

    conn, cur = openproject()

    cur.execute(
        "SELECT pipelineID, FrameNum, totalexp, totalSNR, mode, obsdate, path, pipelinever from datareduction;")# where obsdate between '2016-06-01' and '2017-12-31';")
    rows = cur.fetchall()
    pipelineID = [i[0] for i in rows]
    frameNum = [i[1] for i in rows]
    totalexp = [i[2] for i in rows]
    totalSNR = [i[3] for i in rows]
    mode = [i[4] for i in rows]
    obsdate = [i[5] for i in rows]
    path = [i[6] for i in rows]
    pipelinever = [i[7] for i in rows]

    print("{} data was found in DIB database.".format(len(pipelineID)))

    tellinefile = open("Waveshift_main_ver4/telluric_single_list_ordered_selected3.dat", "r")
    tellinelines = tellinefile.readlines()
    tellinefile.close()
    linecenterlist = [float(i.split()[0]) for i in tellinelines]

    center_wave = [13349., 13042., 12749., 12465., 12195., 11936., 11687., 11451., 11222., 11001., 10791., 10588.,
                   10392., 10204., 10023., 9846., 9675., 9512., 9353., 9200.]

    output = open(sys.argv[1], "a")
    R = 28000
    wsstr1 = "waveshift_measure/"
    vacorair = "VAC"
    fsr = "fsr1.30"

    for i in range(len(pipelineID)):
        if os.path.exists(path[i] + "waveshift_correct/"):
            print(pipelineID[i])
            targetframeglob = glob.glob(path[i] + wsstr1 + "*")
            targetframeglob.sort()
            targetframe = [k.split("/")[-1] for k in targetframeglob]
            targetnormdir = [glob.glob("%s*_%s/%s_norm/%s/" % (path[i], k, vacorair, fsr)) for k in targetframe]
            targetfluxdir = [glob.glob("%s*_%s/%s_flux/%s/" % (path[i], k, vacorair, fsr)) for k in targetframe]
            targetspec = [glob.glob("%s*_%s/%s_norm/%s/*_norm.fits" % (path[i], k, vacorair, fsr)) for k in
                          targetframe]
            targetspecflux = [glob.glob("%s*_%s/%s_flux/%s/*%s.fits" % (path[i], k, vacorair, fsr, vacorair)) for k
                              in targetframe]
            for k in targetspec: k.sort()
            for k in targetspecflux: k.sort()

            orderlist, frameset, filelist = obtainWscorFileList(path[i], fluxornorm="flux", vacorair=vacorair,
                                                                helio=False)
            for f in frameset:
                speclist = []
                for m in orderlist:
                    speclist.append(filelist[f][m])
                peak, centerlams, linewave, centerpix, lambdac, order = Waveshift(speclist, orderlist, linecenterlist,
                                                                                  R, center_wave)
                peak0, centerlams0, linewave0, centerpix0, lambdac0, order0 = Waveshift(
                    targetspecflux[targetframe.index(f)], orderlist, linecenterlist, R, center_wave)
                for n in range(len(linewave)):
                    output.write("%s\t%d\t%s\t%.1f\t%.1f\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.0f\t%.5f\t%.5f\t%.5f\t%.0f\n" % (
                        pipelineID[i], order[n], f, totalexp[i], totalSNR[i], mode[i],
                        pipelinever[i], linewave[n], centerlams[n], centerpix[n], lambdac[n],
                        linewave0[n], centerlams0[n], centerpix0[n], lambdac0[n]))

    # print "The shift as a function of wavelength is determined as follows:"
    # print "SHIFT = %.5e * WAVELENGTH + %.5e" % (inc, offset)
    # print "(Inclination: %.5e)" % inc
    # print "(Offset: %.5e)" % offset

    output.close()
    conn.close()
