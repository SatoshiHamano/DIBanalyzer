# -*- coding:utf-8 -*-

__version__ = "1.2"

import os, glob, shutil
import numpy
import scipy.optimize
import astropy.io.fits as fits
import scipy.constants
from specutils import Spectrum1D
from astropy import units as u

# Description:
#   This script contains some basic and useful functions to manipulate the spectrum files.
#
# Updates:
#   ver1.0 made by Hamano in 2017.01.25
#
#   ver1.1 updated by Hamano in 2018.05.24
#       - pyfits --> astropy.io.fits
#
#   ver1.2 updated by Hamano in 2019.08.24
#       version was moved from file name to __version__ parameter.
#       apall.py was integrated.
#       The name was changed to Spec1Dtools.
#       FSR reading function was integrated.
#       pix2wave was integrated.
#       open_ec_specfiles was integrated from Auto_ecidentify.py.
#


def openspecfits(fitsfile):
    if fitsfile.find("fits") == -1:
        fitsfile += ".fits"

    spfits = fits.open(fitsfile)
    splength = spfits[0].header["NAXIS1"]
    spdata = spfits[0].data

    rcrval1 = float(spfits[0].header["CRVAL1"])
    rcdelt1 = float(spfits[0].header["CDELT1"])
    rcrpix1 = float(spfits[0].header["CRPIX1"])

    lamx = numpy.array([rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)])
    spfits.close()

    return lamx, spdata, rcrval1, rcdelt1, rcrpix1

def savespecfits(spec, inputfits, outputfits):
    spf = fits.open(inputfits)
    spf[0].data = spec.flux.value
    newcrval1 = numpy.amin(spec.wavelength.value)
    deltalam = spec.wavelength.value - numpy.roll(spec.wavelength.value, 1)
    newcdelt1 = numpy.median(deltalam[1:])
    spf[0].header["CRVAL1"] = newcrval1
    spf[0].header["CDELT1"] = newcdelt1
    spf[0].header["CD1_1"] = newcdelt1
    spf[0].header["CRPIX1"] = 1.
    spf.writeto(outputfits)
    spf.close()

def open_ec_specfiles(specfile):
    # open echelle format spectra fits file
    # return the list of x-coordinate, y-coordinate and echelle orders

    spf = fits.open(specfile)
    spdata = spf[0].data
    apnum = len(spdata)

    rcrval1 = float(spf[0].header["CRVAL1"])
    rcdelt1 = float(spf[0].header["CDELT1"])
    rcrpix1 = float(spf[0].header["CRPIX1"])
    splength = spf[0].header["NAXIS1"]

    commonx = numpy.array([[rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)] for i in range(apnum)])

    aporders = []
    for i in range(apnum):
        aporders.append(int(spf[0].header["APNUM%d" % (i + 1)].split()[0]))

    spf.close()

    return commonx, spdata, numpy.array(aporders), rcrval1, rcdelt1, rcrpix1

def FSR_angstrom(select_date="latest"):
    datapath = os.path.dirname(os.path.abspath(__file__))
    basename = "FSR/FSR_winered_"
    extention = "txt"
    datalist = glob.glob("%s/%s*%s" % (datapath, basename, extention))
    dates = [int(i.split(basename)[-1].rstrip(extention).rstrip(".")) for i in datalist]
    latest_index = dates.index(max(dates))

    selected_data = datalist[latest_index]
    if select_date != "latest":
        try:
            for i in range(len(dates)):
                if int(select_date) == dates[i]:
                    selected_data = datalist[i]
                    break
        except:
            print("Warning: %s cannot be converted to integer." % select_date)
        if selected_data == datalist[latest_index]:
            print("Your input \"%s\" is not found. Existing data list:" % select_date)
            for i in range(len(dates)):
                print("\t%s" % datalist[i])
            print("\nLatest data \"%s\" is shown.\n" % dates[latest_index])
    else:
        selected_data = datalist[latest_index]

    print(selected_data)

    datafile = open(selected_data, "r")
    datalines = datafile.readlines()
    datafile.close()

    for i in range(len(datalines)):
        if datalines[i][0] != "#" or datalines[i] != "":
            if datalines[i].find("WIDE") != -1:
                wideline = i
            elif datalines[i].find("HIRES-Y") != -1:
                hiresyline = i
            elif datalines[i].find("HIRES-J") != -1:
                hiresjline = i

    lines = [wideline, hiresyline, hiresjline]
    lines.append(len(datalines) - 1)
    lines.sort()

    fsr_angstrom = {}
    for i in range(len(lines) - 1):
        for j in range(lines[i], lines[i + 1] + 1):
            if len(datalines[j].split()) >= 3:
                linecomp = datalines[j].split()
                order = int(linecomp[0])
                lowlim = float(linecomp[1])
                upplim = float(linecomp[2])
                fsr_angstrom[order] = [lowlim, upplim]

    return fsr_angstrom

def binning_spec(lambdax, fluxy, binning_size):
    lambdax_bin = []
    fluxy_bin = []

    for i in range(len(lambdax) // binning_size):
        tmp_x = 0.
        tmp_y = 0.

        for j in range(binning_size):
            tmp_x += lambdax[binning_size * i + j]
            tmp_y += fluxy[binning_size * i + j]

        lambdax_bin.append(tmp_x / binning_size)
        fluxy_bin.append(tmp_y / binning_size)

    return numpy.array(lambdax_bin), numpy.array(fluxy_bin)



def model_opener(modeltxt, lambdacolumn=0, fluxcolumn=2):
    rf = open(modeltxt, "r")
    rl = rf.readlines()
    rf.close()

    modelx, modely = [], []
    for i in rl:
        if i[0] != "#":
            rl1 = i.split()
            modelx.append(float(rl1[lambdacolumn]))
            modely.append(float(rl1[fluxcolumn]))

    modelx = numpy.array(modelx)
    modely = numpy.array(modely)

    return modelx, modely


def cutSpectrum(inputfits, outputfits, w1, w2):
    inputspec = Spectrum1D.read(inputfits)
    slicesp = inputspec[(w1-1.e-7)*u.AA:(w2+1.e-7)*u.AA]
    savespecfits(slicesp, inputfits, outputfits)

def IRAF_pysarith(spin, ope, factor, spout):
    curdir = os.getcwd()
    outputdir = os.path.dirname(spout)
    outputfile = os.path.basename(spout)
    os.chdir(outputdir)

    counter = 1
    inputfile = os.path.basename(spin)
    copyinputfile = inputfile.replace(".fits", "-cp%d.fits" % counter)
    if not os.path.exists(outputdir + "/" + copyinputfile):
        shutil.copy(spin, outputdir + "/" + copyinputfile)
    factorfile = os.path.basename(factor)
    copyfactor = factorfile.replace(".fits", "-cp%d.fits" % counter)
    if not os.path.exists(outputdir + "/" + copyfactor):
        shutil.copy(factor, outputdir + "/" + copyfactor)


    iraf.sarith(copyinputfile, ope, copyfactor, outputfile)

    os.remove("./%s" % os.path.basename(copyinputfile))
    os.remove("./%s" % os.path.basename(copyfactor))

    os.chdir(curdir)