import sys,os,datetime
import numpy,math,glob,shutil
import pyfits
from datadownload_merlot_ver3 import *
from Waveshift_fit_ver2 import *
from Waveshift_check import *
from PySpecshift import *
from PyScombine import *
from PyContinuum import *

filename = sys.argv[1:]

#Description:
#   This "Waveshift_main.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to automatically download the 1D WINERED spectrum data from
#   merlot server, correct the wavelength shift of each frame caused by the change
#   of the ambient temperature during the observation, combine the shift-corrected
#   fits files of all frames, and continuum fit with arbitraly parameters.
#
#   The correction of wavelength shift is the major purpose of this script.
#   The change of the wavelength is measured by fitting Gaussian functions to the
#   telluric absorption lines, which are enough strong and not blended with other lines.
#   Because the change of the wavelength is currently known to linearly depend on the
#   wavelength, a linear function is fitted to the plot of wavelength - wavelength change.
#   During the fitting, the outlier, which would be originated from the failure of Gaussian
#   fitting, is clipped.
#
#   This script can be applied only to the WIDE mode data of WINERED spectrograph.
#
#Usage:
#
#   $ python Waveshift_main.py <input> <telluric line> <continuum parameter>
#
#   input -- the list of observational information made with the following format
#       ===starting format===
#       starname    year    month   day     observational_run   cut_number
#       .
#       .
#       .
#       ===ending format===
#
#       starname: the name of the star you want to reduce
#       year, month, day: observation date (UT).
#                         year must be given by 4 digits.
#                         month and day must be given by 2 or 1 digits. (e.g., 4, 04, 28)
#       observational_run: specify the observational run of the observation
#                          following expressions are allowed: (2nd, 3rd, 4th, 5th).
#       cut_number: specifty the cut number you want to use.
#                   following expressions are allowed: (cut1, cut2, cut3, cut4, cut5).
#
#   telluric line -- the list of wavelengths of telluric absorption lines used for the
#                    measurement of wavelength shift.
#
#   continuum parameter -- the list of parameters of continuum fitting with the following format
#       ===starting format===
#       echelle_order    low_rej    high_rej   fitting_order     fitting_fuction   not_sampled_region
#       .
#       .
#       .
#       ===ending format===
#
#       echelle_order: the echelle order of spectrum. from m=42 to 61 for WINERED WIDE mode.
#       low_rej: low_rej parameter for IRAF/continuum.
#       high_rej: high_rej parameter for IRAF/continuum.
#       fitting_order: order parameter for IRAF/continuum.
#       fitting_fuction: fitting function for IRAF/continuum.
#       not_sampled_region: specify the region (in angstrom) not included in the continuum fitting.
#                           (e.g.,) 10010:10200,10440:10560
#                           Then, two regions, 10010-10200A and 10440-10560A, are NOT included in fitting.
#                           If you want to include whole range of the spectrum, please type "no" instead.
#
#   After running this script, you will be asked a PASSWORD to access the merlot server.
#
#
#Output:
#   Following directory are created under the current working directory.
#
#   <star1>_<obsdate1>
#       - pipeline_data: contain downloaded data.
#       - specshift: contain wavelength-shift-corrected fits files and fitting result.
#       - scombine: contain combined fits files
#       - normalize: contain normalized fits files
#   <star2>_<obsdate2>
#   .
#   .
#   <starN>_<obsdateN>
#
#
#Updates:
#
#   (ver.2)
#
#   Updated by S.Hamano in 2016.05.09
#
#   Waveshift_check is newly introduced.
#   The result is stored into "status_waveshift_summary.txt", which is created in currrent directory.
#   If "status_waveshift_summary.txt" is exist in current directory, the file is updated.
#
#   Adopting the updated Waveshift_fit (ver2)
#
#   (ver.3)
#
#   Updated by S.Hamano in 2016.05.31
#
#   This script become able to be used in Linux environments by changing the way to use
#   glob.glob() function, which behave differently between Linux and Mac.
#
#   (ver.4)
#
#   Updated by S.Hamano in 2016.09.05
#
#   1) Waveshift become able to be used for the objects whose name is longer than 11 characters.
#   2) datadownload_merlot_ver3.py is installed.
#
##

if __name__ == "__main__":

    filename = sys.argv[1:]
    
    #Parameters

    user = "WINERED"
    server = "merlot.kyoto-su.ac.jp"
    center_wave = [13349., 13042., 12749., 12465., 12195., 11936., 11687., 11451., 11222., 11001., 10791., 10588., 10392., 10204., 10023., 9846., 9675., 9512., 9353., 9200.]
    ordernumber = 20
    orders = range(42,62)
    
    ### Password flag (added in ver.4)
    ### 0(default): you need to type password to access merlot server.
    ### non-zero number: you do not need to type password.
    
    passflag = 0
    #passflag = 1
    
    ###

    #Read the wavelengths of telluric absorption lines

    tellinefile = open(filename[1],"r")
    tellinelines = tellinefile.readlines()
    tellinefile.close()
    linecenterlist = [float(tellinelines[i].split()[0]) for i in range(len(tellinelines))]
    
    #Open the output text file to store the status of Waveshift
    #(added in 2016.05.09 by S.Hamano)
    
    status_waveshift = open("status_waveshift_summary.txt","a")
    status_waveshift.write("Time: %s\n\n" % datetime.datetime.today())
    
    #Ask the Password to access the server

    print "Accessing %s:%s..." % (user, server)
    if passflag == 0:
        password = command_question("Password: ")
    else:
        password = "hoge"
    
    #Obtain required information (e.g., number of frames) by accessing the merlot server.
    
    print "Obtaining required information from %s..." % server
    starname, starnum, obsdate, cut, run, pipelinedirec, framenumber = readstardirec(filename[0], user, server, password, passflag)
    
    #Specify the name of the directories, which will be created.
    
    main_direc = ["%s_%s" % (starname[i], obsdate[i].replace("_","")) for i in range(starnum)]
    target_direc_main = [["%s_%s/pipeline_data%s/frame_NO%d" % (starname[i], obsdate[i].replace("_",""), pipelinedirec[i].split("pipeline")[-1], j+1) for j in range(framenumber[i])] for i in range(starnum)]
    target_direc_sub = [["pipeline_data%s/frame_NO%d" % (pipelinedirec[i].split("pipeline")[-1], j+1) for j in range(framenumber[i])] for i in range(starnum)]
    shift_direc = [["specshift/frame_NO%d" % (j+1) for j in range(framenumber[i])] for i in range(starnum)]
    scombine_direc = ["scombine" for i in range(starnum)]
    continuum_direc = ["normalize" for i in range(starnum)]

    #Specify the spectral resolving power from the infomation of observational run number.

    resolution = []
    for i in range(starnum):
        if run[i] == "4th":
            resolution.append(20000)
        else:
            resolution.append(28000)

    #Download the data

    Data_Download(starname, starnum, target_direc_main, cut, pipelinedirec, framenumber, user, server, password, passflag)

    #Create directories

    for i in range(starnum):
        os.chdir(main_direc[i])
        for j in range(framenumber[i]):
            tmp_outputdir = "specshift/output_data"
            if not os.path.exists(shift_direc[i][j]):
                os.makedirs(shift_direc[i][j])
            if not os.path.exists(tmp_outputdir):
                os.makedirs(tmp_outputdir)
        if not os.path.exists(scombine_direc[i]):
            os.makedirs(scombine_direc[i])
        if not os.path.exists(continuum_direc[i]):
            os.makedirs(continuum_direc[i])
        os.chdir("../")

        print "%s_%s directory and sub directories are created." % (starname[i], obsdate[i].replace("_",""))

    #Measure the wavelength shifts as a function of wavelength and adopt them to the spectrum.

    for i in range(starnum):
        
        os.chdir(main_direc[i])
        
        print "Now correcting the wavelengths..."
        tmp_outputfilelist = [] # added in 2016/05/09 by S.Hamano
        tmp_statusfile = "specshift/status_waveshift.txt"
        for j in range(framenumber[i]):
            tmp_pipelinedatalist = []
            for k in orders:
                tmp_pipelinedatalist.append(glob.glob("%s/*_m%d_*fits" % (target_direc_sub[i][j],k))[0])
            tmp_outputdir = "specshift/output_data"
            tmp_outputfile = "%s/%s_%s_%d.dat" % (tmp_outputdir, starname[i], obsdate[i].replace("_",""), j+1)
            
            tmp_outputfilelist.append(tmp_outputfile)

            inc, offset = Waveshift(tmp_pipelinedatalist, tmp_outputfile, linecenterlist, resolution[i], center_wave)

            for k in range(ordernumber):
                PySpecshift(tmp_pipelinedatalist[k], "%s/%s_%s_No%d_m%d_shift.fits" % (shift_direc[i][j], starname[i], obsdate[i].replace("_",""), j+1, orders[k]), (inc * center_wave[k] + offset) * -1)

        print "Correction of wavelength is done for %s_%s." % (starname[i], obsdate[i].replace("_",""))

    #Check the quality of Waveshift
    #(added in 2016.05.09 by S.Hamano)

        tmp_status_waveshift = Waveshift_check(tmp_outputfilelist, tmp_statusfile)
        status_waveshift.write("%s_%s: %s\n" % (starname[i], obsdate[i].replace("_",""), tmp_status_waveshift))

    #Combine and Continuum fitting

        minlam = []
        maxlam = []
        scombinefits = []
        continuumfits = []

        for k in range(ordernumber):
            tmp_specshift = []
            tmp_specshift_file = []
            for n in range(framenumber[i]):
                tmp_specshift.append(glob.glob("specshift/frame_NO%d/*m%d_shift.fits" % (n+1, orders[k]))[0])
                tmp_specshift_file.append(glob.glob("specshift/frame_NO%d/*m%d_shift.fits" % (n+1, orders[k]))[0].split("/")[-1])
            
            for n in range(len(tmp_specshift)):
                shutil.copy(tmp_specshift[n],scombine_direc[i])

            specshiftdatalist = "%s/%s_%s_m%d_shift.list" % (scombine_direc[i], starname[i], obsdate[i].replace("_",""), orders[k])
            tmp_shiftlist = open(specshiftdatalist, "w")
            scombinefits.append("%s/%s_%s_m%d.fits" % (scombine_direc[i], starname[i], obsdate[i].replace("_",""), orders[k]))
            continuumfits.append("%s/%s_%s_m%dn.fits" % (continuum_direc[i], starname[i], obsdate[i].replace("_",""), orders[k]))

            for j in range(framenumber[i]):
                tmp_shiftlist.write("%s/%s\n" % (scombine_direc[i], tmp_specshift_file[j]))

            tmp_shiftlist.close()

            PyScombine(specshiftdatalist, scombinefits[k])

            for j in range(framenumber[i]):
                os.remove("%s/%s" % (scombine_direc[i], tmp_specshift_file[j]))

            lamx, spdata, rcrval1, rcdelt1, rcrpix1 = openspecfits(scombinefits[k])
            minlam.append(numpy.amin(lamx))
            maxlam.append(numpy.amax(lamx))

        continuum_lowrej, continuum_highrej, continuum_order, continuum_func, continuum_sample = read_continuum_parameters("../%s" % filename[2], minlam, maxlam)

        for k in range(ordernumber):
            PyContinuum(scombinefits[k], continuumfits[k], continuum_lowrej[k], continuum_highrej[k], continuum_order[k], continuum_func[k], continuum_sample[k])

        print "Combine and continuum fitting are done for %s_%s." % (starname[i], obsdate[i].replace("_",""))

        os.chdir("../")


    #Close the output text file
    #(added in 2016.05.09 by S.Hamano)
    status_waveshift.write("\n\n")
    status_waveshift.close()

