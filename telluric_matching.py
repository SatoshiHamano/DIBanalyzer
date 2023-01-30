# -*- coding: utf-8 -*-

import datetime
import sys
import numpy
import math
import mysql.connector
from urllib.parse import urlparse
from open_mysql_project import openproject


def TelluricMatch(objSN, objAM, objUT, telID, telSN, telAM, telUT):
    SNratio = telSN / objSN
    SNtotal = SNratio / (1. + SNratio ** 2) ** 0.5
    SNscore = (SNtotal - 0.5) * 200.

    AMdif = numpy.absolute(objAM - telAM)
    AMscore = 100. - AMdif * 100.

    UTdif = []
    UTscore = []
    for i in telUT:
        td = math.fabs((objUT - i).total_seconds()) / 60. / 60.
        UTdif.append(td)
        UTscore.append(-20. * td + 100.)
    UTscore = numpy.array(UTscore)

    totalScore = SNscore + AMscore + UTscore
    tsid = numpy.argmax(totalScore)
    ScoreDiff = totalScore[tsid] - totalScore
    ScoreDiffMinList = []
    ScoreDiffMinIDList = []
    for i in range(len(ScoreDiff)):
        if ScoreDiff[i] != 0.:
            ScoreDiffMinList.append(ScoreDiff[i])
            ScoreDiffMinIDList.append(telID[i])
    if ScoreDiffMinList != []:
        ScoreDiffMin = ScoreDiffMinList[numpy.argmin(ScoreDiffMinList)]
        ScoreDiffMinID = ScoreDiffMinIDList[numpy.argmin(ScoreDiffMinList)]
    else:
        ScoreDiffMin = 0
        ScoreDiffMinID = None

    # for i in range(len(telID)):
    #     print(telID[i], telSN[i], telAM[i], telUT[i],
    #           "SN%.1f + AM%.1f + UT%.1f = %.1f" % (SNscore[i], AMscore[i], UTscore[i], totalScore[i]))

    return telID[tsid], ScoreDiffMin, ScoreDiffMinID


def TelluricIDList(objidlist):
    conn, cur = openproject()

    ppid = []
    targetSNR = []
    targetmode = []
    targetobsdate = []
    targetslit = []
    for i in range(len(objidlist)):
        cur.execute(
            "SELECT pipelineID,totalSNR,obsdate,mode from datareduction WHERE pipelineID='%s';" % objidlist[i])
        rows = cur.fetchall()

        if len(rows) == 1:
            ppid.append(rows[0][0])
            targetSNR.append(float(rows[0][1]))
            targetobsdate.append(rows[0][2])
            targetmode.append(rows[0][3])
        else:
            print(rows)
            sys.exit()

        cur.execute(
            "SELECT obs.slit from reducedframe as fr JOIN observation as obs ON fr.objectframe=obs.frame WHERE fr.pipelineID='%s';" %
            ppid[i])
        rows = cur.fetchall()

        slitlist = [rows[j][0] for j in range(len(rows))]
        slitset = set(slitlist)
        slits = list(slitset)
        slits.sort()
        slitstr = ""
        for j in slits:
            slitstr += str(j)
        targetslit.append(slitstr)

    BestID_list = []
    Dif1st2nd_list = []
    SecondBestID_list = []
    for i in range(len(ppid)):
        cur.execute(
            "SELECT red.pipelineID,red.totalSNR,red.mode from datareduction as red JOIN object as obj USING(objectID) WHERE red.obsdate='%s' AND obj.type='STANDARD';" %
            targetobsdate[i])
        rows = cur.fetchall()

        if rows != []:
            stppid = [rows[n][0] for n in range(len(rows))]
            tel_snr_all = numpy.array([float(rows[n][1]) for n in range(len(rows))])
            tel_mode = [rows[n][2] for n in range(len(rows))]

            cur.execute(
                "SELECT obs.airmass,obs.obsdate from reducedframe AS fr JOIN observation AS obs ON fr.objectframe=obs.frame WHERE fr.pipelineID='%s';" %
                ppid[i])
            rows = cur.fetchall()
            medid = int(len(rows) / 2) if len(rows) % 2 == 0 else int((len(rows) - 1) / 2)
            airmass = [float(rows[n][0]) for n in range(len(rows))]
            obsdate = [rows[n][1] for n in range(len(rows))]

            target_am = numpy.average(airmass)
            target_obstime = obsdate[medid]

            tel_id = []
            tel_am = []
            tel_obstime = []
            tel_snr = []
            for j in range(len(stppid)):
                cur.execute(
                    "SELECT obs.airmass,obs.obsdate,obs.slit from reducedframe AS fr JOIN observation AS obs ON fr.objectframe=obs.frame WHERE fr.pipelineID='%s';" %
                    stppid[j])
                rows = cur.fetchall()
                airmass = [float(rows[n][0]) for n in range(len(rows))]
                obsdate = [rows[n][1] for n in range(len(rows))]
                telslit = [rows[n][2] for n in range(len(rows))]
                telslit_set = set(telslit)
                telslit_list = list(telslit_set)
                telslit_list.sort()
                telslitstr = ""
                for n in telslit_list:
                    telslitstr += str(n)

                if targetmode[i] == tel_mode[j] and targetslit[i] == telslitstr:
                    am_ave = numpy.average(airmass)
                    medid = int(len(rows) / 2) if len(rows) % 2 == 0 else int((len(rows) - 1) / 2)
                    obsdate_med = obsdate[medid]
                    tel_id.append(stppid[j])
                    tel_am.append(am_ave)
                    tel_obstime.append(obsdate_med)
                    tel_snr.append(tel_snr_all[j])

            tel_am = numpy.array(tel_am)
            tel_obstime = numpy.array(tel_obstime)
            tel_snr = numpy.array(tel_snr)

            BestID, Dif1st2nd, SecondBestID = TelluricMatch(targetSNR[i], target_am, target_obstime, tel_id, tel_snr,
                                                            tel_am, tel_obstime)

            BestID_list.append(BestID)
            Dif1st2nd_list.append(Dif1st2nd)
            SecondBestID_list.append(SecondBestID)
            # print(ppid[i], targetSNR[i], target_am, target_obstime, targetmode[i], targetslit[i])
            # print(TelluricMatch(targetSNR[i], target_am, target_obstime, tel_id, tel_snr, tel_am, tel_obstime))
            # print("\n")
        else:
            BestID_list.append(None)
            Dif1st2nd_list.append(0)
            SecondBestID_list.append(None)

    return BestID_list, Dif1st2nd_list, SecondBestID_list

if __name__ == "__main__":
    filename = sys.argv[1:]
    rf = open(filename[0], "r")
    rl = rf.readlines()
    rf.close()

    objidlist = [i.split()[0] for i in rl]
    BestID_list, Dif1st2nd_list, SecondBestID_list = TelluricIDList(objidlist)
    for i in range(len(objidlist)):
        print(objidlist[i], BestID_list[i], Dif1st2nd_list[i], SecondBestID_list[i])

    telidlist = []
    for i in range(len(objidlist)):
        if not BestID_list[i] is None:
            telidlist.append(BestID_list[i])
        if not SecondBestID_list[i] is None:
            if Dif1st2nd_list[i] < 30:
                telidlist.append(SecondBestID_list[i])

    telidlist = list(set(telidlist))
    telidlist.sort()

    for i in telidlist:
        print(i)