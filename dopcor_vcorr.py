import sys, os, datetime
import astropy.io.fits as fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
# from pyraf import iraf
from open_mysql_project import openproject
from specutils import Spectrum1D
import scipy.constants
from Spec1Dtools import savespecfits

# iraf.images()
# iraf.imutil()
# iraf.imgeom()


def IRAF_dopcor_vcorr(inputsp, outputsp, v):
    # inputcopy = "./work_place/" + inputsp.split("/")[-1].replace(".fits", "-cp.fits")
    # outputcopy = "./work_place/" + outputsp.split("/")[-1].replace(".fits", "-cp.fits")
    # shutil.copyfile(inputsp, inputcopy)
    #
    # iraf.dopcor(inputcopy, outputcopy, redshift="%.2f" % v, isvelocity="yes", add="yes")
    #
    # shutil.copyfile(outputcopy, outputsp)
    # os.remove(inputcopy)
    # os.remove(outputcopy)

    outputdir = os.path.dirname(outputsp)
    tempname = outputdir + "/tmp_dopcor_vcorr.fits"

    iraf.dopcor(inputsp, tempname, redshift="%.2f" % v, isvelocity="yes", add="yes")

    os.rename(tempname, outputsp)

def dopcor_vcorr(inputsp, outputsp, v):
    inputspec = Spectrum1D.read(inputsp)
    shiftedspec = Spectrum1D(spectral_axis=inputspec.wavelength * (1.-v/(scipy.constants.c*1.e-3)), flux=inputspec.flux)
    savespecfits(shiftedspec, inputsp, outputsp)


def read_rvfile(rvfile, mode="helio"):
    rvf = open(rvfile, "r")
    rvl = rvf.readlines()
    rvf.close()

    v = "INDEF"

    for i in rvl:
        if i[0] != "#":
            if mode == "helio":
                v = float(i.split()[2])
            elif mode == "lsr":
                v = float(i.split()[3])

    return v

def read_rvfile_ID(ID):
    rvfile = "rvcorrect.txt"
    conn, cur = openproject()
    cur.execute(
        "select y.path from telluriccorrection as x JOIN datareduction as y on x.pipelineIDobj=y.pipelineID where x.telluricID='%s';" % ID)
    rows = cur.fetchall()
    if rows == []:
        cur.execute(
            "select path from datareduction where pipelineID='%s';" % ID)
        rows = cur.fetchall()
        if rows == []:
            print("%s is not registered." % ID)
            sys.exit()
        elif len(rows) == 1:
            path = rows[0][0]
            return read_rvfile(path + rvfile)
    else:
        path = rows[0][0]
        return read_rvfile(path + rvfile)

def dopcor_rvcorrect(inputsp, outputsp, rvfile, mode="helio"):
    v = read_rvfile(rvfile, mode=mode)
    if v == "INDEF":
        print("Something is wrong!")
    else:
        dopcor_vcorr(inputsp, outputsp, -v)


def IRAF_rvcorrect(inputfits, output, observatory="INDEF"):

    thres_date = datetime.date(2017, 1, 1)

    spfits = fits.open(inputfits)
    tmpdec = spfits[0].header["DEC"].split(":")
    tmpra = spfits[0].header["RA"].split(":")
    tmput = spfits[0].header["ACQTIME1"].split("-")
    spfits.close()

    decpm = float(tmpdec[0]) / math.fabs(float(tmpdec[0]))

    dec = float(tmpdec[0]) + float(tmpdec[1]) / 60. * decpm + float(tmpdec[2]) / 3600. * decpm
    # dec = float(tmpdec[0]) + float(tmpdec[1]) / 60. + float(tmpdec[2]) / 3600.
    ra = float(tmpra[0]) + float(tmpra[1]) / 60. + float(tmpra[2]) / 3600.

    year = tmput[0]
    month = tmput[1]
    day = tmput[2]
    hour = float(tmput[3].split(":")[0]) + float(tmput[3].split(":")[1]) / 60. + float(tmput[3].split(":")[2]) / 3600.

    obsdate = datetime.date(int(year), int(month), int(day))

    if observatory == "INDEF":
        if obsdate < thres_date:
            observatory = "KAO"
        else:
            observatory = "lasilla"

    fo = open(output, "w")
    sys.stdout = fo

    iraf.rvcorrect(observatory=observatory, year=year, month=month, day=day, ut="%.2f" % hour, ra="%.2f" % ra,
                   dec="%.2f" % dec)

    fo.close()
    sys.stdout = sys.__stdout__


def rvcorrect(inputfits, output, observatory="INDEF"):
    thres_date = Time("2017-1-1")

    spfits = fits.open(inputfits)
    sc = SkyCoord(ra=spfits[0].header["RA"], dec=spfits[0].header["DEC"], unit=(u.hourangle, u.deg))
    acqtime1 = spfits[0].header["ACQTIME1"]
    acqtime = Time(acqtime1[:10] + "T" + acqtime1[11:])
    spfits.close()

    if observatory == "INDEF":
        if acqtime < thres_date:
            observatory = "KAO"
        else:
            observatory = "lasilla"

    kao = EarthLocation.from_geodetic(lat=35.0703*u.deg, lon=135.758*u.deg, height=136*u.m)
    lasilla = EarthLocation.from_geodetic(lat=-29.25666667*u.deg, lon=-70.73*u.deg, height=2347*u.m)

    fo = open(output, "w")
    fo.write("# Calculated using astropy.coordinates.SkyCoords.radial_velocity_correction\n")

    if observatory == "KAO":
        helcorr = sc.radial_velocity_correction("heliocentric", obstime=acqtime, location=kao)
        fo.write("# \tlatitude = {}\n# \tlongitude = {}\n# \taltitude = {}\n".format(kao.lat, kao.lon, kao.height))
        fo.write("## HJD\t\tVOBS\tVHELIO\tVLSR\n")
        fo.write("{:.5f}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(acqtime.jd, 0.0, helcorr.to(u.km/u.s).value, 999.9))
    elif observatory == "lasilla":
        helcorr = sc.radial_velocity_correction("heliocentric", obstime=acqtime, location=lasilla)
        fo.write("# \tlatitude = {}\n# \tlongitude = {}\n# \taltitude = {}\n".format(lasilla.lat, lasilla.lon, lasilla.height))
        fo.write("## HJD\t\tVOBS\tVHELIO\tVLSR\n")
        fo.write("{:.5f}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(acqtime.jd, 0.0, helcorr.to(u.km/u.s).value, 999.9))

    fo.close()


if __name__ == "__main__":
    inputsp = sys.argv[1]
    outputsp = sys.argv[2]
    outputf = sys.argv[3]

    rvcorrect(inputsp, outputf)
    dopcor_rvcorrect(inputsp, outputsp, outputf)


#
# filename = sys.argv[1:]
#
# fitslist = open(filename[0], "r")
# fitslines = fitslist.readlines()
#
# spfiles = []
# vhelio = []
# for i in range(len(fitslines)):
#     spfiles.append(fitslines[i].split()[0])
#
# starlist = open(filename[1], "r")
# starlines = starlist.readlines()
#
# star = []
# for i in range(len(starlines)):
#     star.append(starlines[i].split()[0])
#     vhelio.append(float(starlines[i].split()[4]) - float(starlines[i].split()[3]))
#
# fitslist.close()
#
# for i in range(len(spfiles)):
#     for j in range(len(star)):
#         if spfiles[i].find(star[j]) != -1:
#             print
#             spfiles[i], star[j], vhelio[j]
#             iraf.dopcor(spfiles[i], spfiles[i].rstrip("fits").rstrip(".") + "v.fits", redshift="%.2f" % (vhelio[j]),
#                         isvelocity="yes", add="yes")
