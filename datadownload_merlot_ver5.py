#!/usr/bin/env python
# -*- coding: utf-8 -*-
## modules
import os, time, glob
import pexpect
from open_mysql_project import openproject


# Description:
#   This "datadownload_merlot.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to automatically download the 1D WINERED spectrum data from
#   merlot server.
#
# Usage:
#
#   $ python datadownload_merlot.py <input>
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
#   After running this script, you will be asked a PASSWORD to access the merlot server.
#
#
# Output:
#   Following directory are created under the current working directory.
#
#   <star1>_<obsdate1>
#       - pipeline_data: contain downloaded data.
#
# Updates:
#
#   (ver.2)
#
#   The version number of pipeline is added to the end of directory name.
#   The "Data_Download_sum" function is newly added.
#
#   (ver.3)
#
#   Minor bugs fixed.
#
##


def command_question(question):
    answer = input(question)
    return answer


def Obtain_ls_Result(user, server, password, directory_path, lsinput, waittime=2., lsoption="-1d", pflag=0):
    p = pexpect.spawn("ssh %s@%s" % (user, server))
    time.sleep(waittime)
    if pflag == 0:
        p.expect("%s@%s's password:" % (user, server))
        p.sendline(password)
    p.expect_exact("$")
    p.sendline("cd %s" % directory_path)
    p.expect_exact("$")
    p.sendline("ls %s --color=no %s" % (lsoption, lsinput))
    time.sleep(waittime)
    p.expect_exact("$")
    resultls = p.before
    p.sendline("exit")
    p.close()

    return resultls.split("\r\n")[1:-1]


def Pipeline_info(user, server, password, waittime=2., pflag=0):
    p = pexpect.spawn("ssh %s@%s" % (user, server))
    time.sleep(waittime)
    if pflag == 0:
        p.expect("%s@%s's password:" % (user, server))
        p.sendline(password)
    p.expect_exact("$")
    p.sendline("python /media/WD_ext/temporary/hamano/pipeline_info_hamano_ver2.py")
    p.expect_exact("$")
    p.sendline("exit")

    p = pexpect.spawn(
        "scp %s@%s:/media/WD_ext/temporary/hamano/pipeline_info_hamano.dat /Users/hamano/DIB_analysis/script/asciidata/" % (
            user, server))
    if pflag == 0:
        p.expect("%s@%s's password:" % (user, server))
        p.sendline(password)

    p.interact()
    p.close()


def Obtain_scp_Files(user, server, password, download_data, target_direc, waittime=2., pflag=0):
    p = pexpect.spawn("scp %s@%s:%s %s" % (user, server, download_data, target_direc))
    time.sleep(waittime)
    if pflag == 0:
        p.expect("%s@%s's password:" % (user, server))
        p.sendline(password)
    p.interact()
    p.close()


def Obtain_scp_Directory(user, server, password, download_data, target_direc, waittime=2., pflag=0):
    p = pexpect.spawn("scp -r %s@%s:%s %s" % (user, server, download_data, target_direc))
    time.sleep(waittime)
    if pflag == 0:
        p.expect("%s@%s's password:" % (user, server))
        p.sendline(password)
    p.interact()
    p.close()


def Data_Download(starname, starnum, target_direc, cut, pipelinedirec, framenumber, user, server, password, passflag):
    for i in range(starnum):
        for j in range(framenumber[i]):
            download_data = "%s/%s_NO%d/onedspec/object/flux/%s/*fits" % (pipelinedirec[i], starname[i], j + 1, cut[i])
            if not os.path.exists(target_direc[i][j]):
                os.makedirs(target_direc[i][j])
            Obtain_scp_Files(user, server, password, download_data, target_direc[i][j], pflag=passflag)


def Data_Download_sum(starname, starnum, target_direc, cut, pipelinedirec, user, server, password, passflag):
    for i in range(starnum):
        download_data = "%s/%s_sum_data/onedspec/object/flux/%s/*fits" % (pipelinedirec[i], starname[i], cut[i])
        if not os.path.exists(target_direc[i]):
            os.makedirs(target_direc[i])
        Obtain_scp_Files(user, server, password, download_data, target_direc[i], pflag=passflag)


def collect_info(pipeline_direc):
    objectfits = []
    skyfits = []
    bgregion = []

    print(pipeline_direc)
    print(glob.glob(pipeline_direc + "*.txt"))

    if len(glob.glob(pipeline_direc + "*.txt")) == 1:
        listf = open(glob.glob(pipeline_direc + "*.txt")[0], "r")
        listl = listf.readlines()
        listf.close()
    elif len(glob.glob(pipeline_direc + "*.txt")) > 1:
        if len(glob.glob(pipeline_direc + "*list.txt")) == 1:
            listf = open(glob.glob(pipeline_direc + "*list.txt")[0], "r")
            listl = listf.readlines()
            listf.close()
        elif len(glob.glob(pipeline_direc + "*list_edit.txt")) == 1:
            listf = open(glob.glob(pipeline_direc + "*list_edit.txt")[0], "r")
            listl = listf.readlines()
            listf.close()
        else:
            print("FAILED: List is not found")
            return False
    else:
        print("FAILED: List is not found")
        return False

    for i in range(len(listl)):
        if len(listl[i]) > 3:
            listl1line = listl[i].split()
            objectfits.append(listl1line[0])
            skyfits.append(listl1line[1])
            if len(listl1line) > 2:
                bgregion.append(listl1line[-1])
            else:
                bgregion.append("INDEF")

    aperture = []
    badpix = []
    waveshift_ave = []
    waveshift_stv = []
    waveshift_num = []
    UTstart = []
    UTend = []
    RA = []
    Dec = []
    AirmassStart = []
    AirmassEnd = []
    pipelineparameter = {}
    calibdata = {}

    if os.path.exists(pipeline_direc + "reduction_log/"):
        if os.path.exists(pipeline_direc + "reduction_log/aperture_log.txt"):
            aplogf = open(pipeline_direc + "reduction_log/aperture_log.txt", "r")
            aplogl = aplogf.readlines()
            aplogf.close()
            for i in aplogl:
                if i.find("m=45") != -1:
                    aplogl1line = i.split()
                    for j in range(1, len(aplogl1line)):
                        aperture.append(aplogl1line[j])
                elif i.find("m=142") != -1:
                    aplogl1line = i.split()
                    for j in range(1, len(aplogl1line)):
                        aperture.append(aplogl1line[j])
                elif i.find("m=171") != -1:
                    aplogl1line = i.split()
                    for j in range(1, len(aplogl1line)):
                        aperture.append(aplogl1line[j])
        else:
            print("FAILED: Aperture log file is not found.")
            return False

        if os.path.exists(pipeline_direc + "reduction_log/badpix_log.txt"):
            bplogf = open(pipeline_direc + "reduction_log/badpix_log.txt", "r")
            bplogl = bplogf.readlines()
            bplogf.close()
            for i in bplogl:
                if i.find("No.") != -1:
                    bplogl1line = i.split()
                    if len(bplogl1line) == 3:
                        badpix.append(int(bplogl1line[1]))

        if len(objectfits) == 1:
            waveshift_ave.append(0.)
            waveshift_stv.append(0.)
            waveshift_num.append(0)
        elif os.path.exists(pipeline_direc + "reduction_log/waveshift_log.txt"):
            wslogf = open(pipeline_direc + "reduction_log/waveshift_log.txt", "r")
            wslogl = wslogf.readlines()
            wslogf.close()
            for i in wslogl:
                if i.find("Average") != -1:
                    wslogl1line = i.split(":")[1].split()
                    for j in range(len(wslogl1line)):
                        waveshift_ave.append(float(wslogl1line[j]))
                if i.find("Std dev") != -1:
                    wslogl1line = i.split(":")[1].split()
                    for j in range(len(wslogl1line)):
                        waveshift_stv.append(float(wslogl1line[j]))
                if i.find("Number of orders") != -1:
                    wslogl1line = i.split(":")[1].split()
                    for j in range(len(wslogl1line)):
                        waveshift_num.append(int(wslogl1line[j]))
        else:
            print("FAILED: Waveshift log file is not found.")
            return False

        if len(glob.glob(pipeline_direc + "reduction_log/*.tex")) > 0:
            texlogf = open(glob.glob(pipeline_direc + "reduction_log/*.tex")[0], "r")
            texlogl = texlogf.readlines()
            texlogf.close()
            reductiondataflag = False
            pipelineflag = False
            calibflag = False
            for i in texlogl:
                if i.find("section*{Reduced data}") != -1:
                    reductiondataflag = True
                if reductiondataflag:
                    for j in objectfits:
                        if i.replace("\\_", "_").find(j) != -1:
                            textable1line = i.split("&")
                            if len(textable1line) > 8:
                                UTstart.append(textable1line[3].split()[0])
                                UTend.append(textable1line[4].split()[0])
                                RA.append(textable1line[5].split()[0])
                                Dec.append(textable1line[6].split()[0])
                                AirmassStart.append(textable1line[7].split()[0])
                                AirmassEnd.append(textable1line[8].split()[0])
                if i.find("section*{Pipeline information}") != -1:
                    reductiondataflag = False
                    pipelineflag = True
                if pipelineflag:
                    if i.find("ver.") != -1:
                        pipelineparameter["ver"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Scattered light subtraction") != -1:
                        pipelineparameter["ScatteredLight"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Manual aperture range setting") != -1:
                        pipelineparameter["Manual"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Background subtraction mode") != -1:
                        pipelineparameter["Background"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Hot pix correction") != -1:
                        pipelineparameter["HotPix"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Sky emission spectra") != -1:
                        pipelineparameter["SkyEmission"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Transform dy") != -1:
                        pipelineparameter["Transformdy"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Transform flux") != -1:
                        pipelineparameter["Transformflux"] = i.split("&")[1].rstrip("\\\n")
                    if i.find("Cut ranges (x FSRs)") != -1:
                        pipelineparameter["CutRange"] = i.split("&")[1].rstrip("\\\n")
                if i.find("section*{Calibration data}") != -1:
                    pipelineflag = False
                    calibflag = True
                if calibflag:
                    if i.find("Flat fielding image") != -1:
                        calibdata["Flat"] = i.split("&")[1].rstrip("\\\n").replace("\\_", "_")
                    if i.find("Bad pixel mask") != -1:
                        calibdata["BPMask"] = i.split("&")[1].rstrip("\\\n").replace("\\_", "_")
                    if i.find("Comparison image") != -1:
                        calibdata["Comparison"] = i.split("&")[1].rstrip("\\\n").replace("\\_", "_")
                    if i.find("Aperture trace image") != -1:
                        calibdata["Trace"] = i.split("&")[1].rstrip("\\\n").replace("\\_", "_")
                    if i.find("Aperture mask for apscatter") != -1:
                        calibdata["ApscatterMask"] = i.split("&")[1].rstrip("\\\n").replace("\\_", "_")
                if i.find("section*{Aperture range and Linear shift of spectra}") != -1:
                    calibflag = False
        else:
            print("FAILED: tex file is not found.")
            return False
    else:
        print("FAILED: reduction_log/ was not found.")
        return False

    snrflag = True
    if len(objectfits) == 1:
        snratio = 0.
    elif os.path.exists(pipeline_direc + "SNR_dat/fsr1.05/SNratio_m45_fsr1.05.dat"):
        snrf = open(pipeline_direc + "SNR_dat/fsr1.05/SNratio_m45_fsr1.05.dat", "r")
    elif os.path.exists(pipeline_direc + "SNR_dat/fsr1.05/SNratio_m142_fsr1.05.dat"):
        snrf = open(pipeline_direc + "SNR_dat/fsr1.05/SNratio_m142_fsr1.05.dat", "r")
    elif os.path.exists(pipeline_direc + "SNR_dat/fsr1.05/SNratio_m171_fsr1.05.dat"):
        snrf = open(pipeline_direc + "SNR_dat/fsr1.05/SNratio_m171_fsr1.05.dat", "r")
    else:
        snrflag = False
        print("FAILED: SNR dat file is not found.")
        snratio = 0

    if len(objectfits) > 1 and snrflag:
        snrl = snrf.readlines()
        snrf.close()

        for i in snrl:
            if i.find("S/N (total):") != -1:
                snratio = float(i.split(":")[1].split()[1])

    # for i in range(len(objectfits)):
    #     framedirec = glob.glob(pipeline_direc+"*_NO%d" % (i+1))[0]
    #     specfitsfiles = glob.glob(framedirec+"/onedspec/VAC_flux/fsr1.05/*_m45*fits")
    #     if len(specfitsfiles) > 0:
    #         _, spy, _,_,_ = openspecfits(specfitsfiles[0])
    #     else:
    #         specfitsfiles = glob.glob(framedirec + "onedspec/VAC_flux/fsr1.05/*_m142*fits")
    #         if len(specfitsfiles) > 0:
    #             _, spy, _,_,_ = openspecfits(specfitsfiles[0])
    #         else:
    #             specfitsfiles = glob.glob(framedirec + "onedspec/VAC_flux/fsr1.05/*_m171*fits")
    #             if len(specfitsfiles) > 0:
    #                 _, spy, _, _, _ = openspecfits(specfitsfiles[0])
    #             else:
    #                 return False
    #
    #     specaverage.append(numpy.average(spy))

    return [objectfits, skyfits, bgregion, aperture, badpix, waveshift_ave, waveshift_stv, waveshift_num, UTstart,
            UTend, RA, Dec, AirmassStart, AirmassEnd, pipelineparameter, calibdata, snratio]


if __name__ == "__main__":
    user = "WINERED"
    server = "merlot.kyoto-su.ac.jp"

    print("Accessing %s:%s..." % (user, server))

    ### Password flag (added in ver.3)
    ### 0(default): you need to type password to access merlot server.
    ### non-zero number: you do not need to type password.

    passflag = 0

    ###
    if passflag == 0:
        password = command_question("Password: ")
    else:
        password = "hoge"

    Pipeline_info(user, server, password)

    rf = open("/Users/hamano/DIB_analysis/script/asciidata/pipeline_info_hamano.dat", "r")
    rl = rf.readlines()
    rf.close()

    wf = open("/Users/hamano/DIB_analysis/script/asciidata/pipeline_info_error.dat", "w")

    redIDfull = [rl[i].split(" || ")[0] for i in range(len(rl))]
    objframefull = [rl[i].split(" || ")[2] for i in range(len(rl))]
    skyframefull = [rl[i].split(" || ")[3] for i in range(len(rl))]
    redID = []
    objframe = {}
    skyframe = {}
    fnumlist = {}
    for i in range(len(redIDfull)):
        if not redIDfull[i] in redID:
            redID.append(redIDfull[i])
            objframe[redIDfull[i]] = [objframefull[i]]
            skyframe[redIDfull[i]] = [skyframefull[i]]
            fnumlist[redIDfull[i]] = 1
        else:
            objframe[redIDfull[i]].append(objframefull[i])
            skyframe[redIDfull[i]].append(skyframefull[i])
            fnumlist[redIDfull[i]] += 1

    redDate = [redID[i][0:10] for i in range(len(redID))]
    redYear = [redDate[i][0:4] for i in range(len(redID))]
    redDir = [redID[i][11:] for i in range(len(redID))]
    redVer = [redDir[i].split("_ver")[1] for i in range(len(redID))]
    redObj = [redDir[i].split("_pipeline_ver")[0] for i in range(len(redID))]

    conn, cur = openproject()

    cur.execute("select frame,objectid,type,exptime,instmode from observation;")
    rows = cur.fetchall()
    if rows == []:
        framelist = []
        objidobslist = []
        typelist = []
        explist = []
        modelist = []
    else:
        framelist = [rows[i][0] for i in range(len(rows))]
        objidobslist = [rows[i][1] for i in range(len(rows))]
        typelist = [rows[i][2] for i in range(len(rows))]
        explist = [float(rows[i][3]) for i in range(len(rows))]
        modelist = [rows[i][4] for i in range(len(rows))]

    cur.execute("select objectid,objectname from object;")
    rows = cur.fetchall()
    if rows == []:
        objidlist = []
        namelist = []
    else:
        objidlist = [rows[i][0] for i in range(len(rows))]
        namelist = [
            rows[i][1].replace("*", "").lstrip(" ").replace("   ", " ").replace("  ", " ").replace(" ", "_").replace(
                "[", "").replace("]", "") for i in range(len(rows))]

    cur.execute("select pipelineID from datareduction;")
    rows = cur.fetchall()
    if rows == []:
        pipidlist = []
    else:
        pipidlist = [rows[i][0] for i in range(len(rows))]

    objpathcore = "/Users/hamano/DIB_analysis/DIB_pipeline_dir"
    stdpathcore = "/Users/hamano/DIB_analysis/Standard_pipeline_dir"
    pippath = "/media/WD_ext/PIPELINE_DATA"

    calibkeys = ["Flat", "BPMask", "Comparison", "Trace", "ApscatterMask"]
    pipelinekeys = ["ver", "ScatteredLight", "Manual", "Background", "HotPix", "SkyEmission", "Transformdy",
                    "Transformflux", "CutRange"]

    conv = {"yes": "true", "no": "false"}

    for i in range(len(redID)):
        if not redID[i] in pipidlist:
            objids = set()
            totalexp = 0.
            for j in objframe[redID[i]]:
                if j in framelist:
                    objids.add(objidobslist[framelist.index(j)])
                    totalexp += explist[framelist.index(j)]
            # print(redID[i],objids,totalexp)
            if len(objids) > 1:
                idstr = ""
                for k in objids:
                    idstr += str(k)
                    idstr += ", "
                idstr.rstrip(", ")
                wf.write(idstr + " is included in the dataset %s.\n" % redID[i])
            elif len(objids) == 1:
                objidPL = list(objids)[0]
                objnamePL = namelist[objidlist.index(objidPL)]
                if not objframe[redID[i]][0] in framelist:
                    continue
                if typelist[framelist.index(objframe[redID[i]][0])] == "OBJECT":
                    targetdirec = "%s/%s/%s/" % (objpathcore, objnamePL, redID[i])
                    if not os.path.exists("%s/%s" % (objpathcore, objnamePL)):
                        os.makedirs("%s/%s" % (objpathcore, objnamePL))
                else:
                    targetdirec = "%s/%s/%s/" % (stdpathcore, redDate[i], redID[i])
                    if not os.path.exists("%s/%s" % (stdpathcore, redDate[i])):
                        os.makedirs("%s/%s" % (stdpathcore, redDate[i]))
                if not os.path.exists(targetdirec):
                    Obtain_scp_Directory(user, server, password, "%s/%s/%s/%s_summary/%s_small/" % (
                        pippath, redYear[i], redDate[i], redDate[i], redDir[i]), targetdirec)
                if not os.path.exists(targetdirec + "SNR_dat/"):
                    Obtain_scp_Directory(user, server, password, "%s/%s/%s/%s/SNR_dat/" % (
                        pippath, redYear[i], redDate[i], redDir[i]), targetdirec)
                if os.path.exists(targetdirec + "framelist.par"):
                    rflist = open(targetdirec + "framelist.par", "r")
                    rllist = rflist.readlines()
                    rflist.close()
                    if len(rllist) != 0:
                        frameflag = True
                    else:
                        frameflag = False
                else:
                    frameflag = False
                if not frameflag:
                    Obtain_scp_Files(user, server, password, "%s/%s/%s/%s/*.txt" % (
                        pippath, redYear[i], redDate[i], redDir[i]), targetdirec)
                    if not os.path.exists(targetdirec):
                        continue
                    wflist = open(targetdirec + "framelist.par", "w")
                    wflist.write("The list file is taken from server.")
                    wflist.close()

                info = collect_info(targetdirec)
                if info != False:
                    [objectfits, skyfits, bgregion, aperture, badpix, waveshift_ave, waveshift_stv, waveshift_num,
                     UTstart,
                     UTend, RA, Dec, AirmassStart, AirmassEnd, plpara, calibdata, snratio] = info
                    calibtext = ""
                    plparatext = ""
                    for k in calibkeys:
                        calibtext += "'%s'," % calibdata[k]
                    for k in pipelinekeys:
                        if k == "Background":
                            plparatext += "'%s'," % plpara[k]
                        elif k == "Transformdy":
                            plparatext += "%s," % plpara[k]
                        elif k == "CutRange":
                            plparatext += "'%s'" % plpara[k].replace(",", " ")
                        elif k == "ver":
                            plparatext += "'%s'," % plpara[k]
                        else:
                            plparatext += "%s," % conv[plpara[k].lstrip().rstrip()]

                    valuetext = "'%s',%s,%s,%s,%s,'%s',cast('%s' as datetime),'%s',%s%s" % (
                        redID[i], objidPL, fnumlist[redID[i]], totalexp, snratio,
                        modelist[framelist.index(objframe[redID[i]][0])], redDate[i], targetdirec, calibtext,
                        plparatext.rstrip(",")
                    )

                    cur.execute(
                        "INSERT IGNORE INTO datareduction (pipelineID,objectID,FrameNum,totalexp,totalSNR,mode,obsdate,path,flat,bpmask,comparison,trace,apsmask,pipelinever,scatteredlight,manual,background,hotpix,skyemission,transformdy,transformflux,cutrange) VALUES (%s);" % valuetext)
                    for m in range(len(objectfits)):
                        if objectfits[m].find("WINA") != -1:
                            objectfits[m] = objectfits[m].lstrip("CDS_WINA0")
                        if skyfits[m].find("WINA") != -1:
                            skyfits[m] = skyfits[m].lstrip("CDS_WINA0")
                        if objectfits[m] in framelist:
                            frametext = "'%s','%s','%s',%s,'%s','%s',%s,%s,%s,%s" % (
                                redID[i], objectfits[m], skyfits[m], explist[framelist.index(objectfits[m])],
                                aperture[m],
                                bgregion[m], badpix[m], waveshift_ave[m], waveshift_stv[m], waveshift_num[m])
                        else:
                            frametext = "'%s','%s','%s',%s,'%s','%s',%s,%s,%s,%s" % (
                                redID[i], objectfits[m], skyfits[m], str(0.), aperture[m],
                                bgregion[m], badpix[m], waveshift_ave[m], waveshift_stv[m], waveshift_num[m])

                        cur.execute(
                            "INSERT IGNORE INTO reducedframe (pipelineID,objectframe,skyframe,exptime,aperture,background,badpix,wsave,wsstd,wsnum) VALUES(%s);" % frametext)
                else:
                    wf.write("Warning: Could not retrieve pipeline info from %s.\n" % redID[i])

            conn.commit()

    conn.close()

    wf.close()

    # targetDir = filename[1]
    # for i in range(len(rl)):
    #     if os.path.exists("%s/%s/%s" % (targetDir, redDate[i], redDir[i])):
    #         print "You have %s in %s!" % (redDir[i], redDate[i])
    #     else:
    #         if not os.path.exists("%s/%s" % (targetDir, redDate[i])):
    #             os.makedirs("%s/%s" % (targetDir, redDate[i]))
    #         Obtain_scp_Directroy(user, server, password, "%s/%s/%s/%s_summary/%s_small/" % (
    #         pippath, redYear[i], redDate[i], redDate[i], redDir[i]),
    #                              "%s/%s" % (targetDir, redDate[i]))
    #
