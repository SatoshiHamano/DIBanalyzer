import sys,os,time
import pexpect

#Description:
#   This "datadownload_merlot.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to automatically download the 1D WINERED spectrum data from
#   merlot server.
#
#Usage:
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
#Output:
#   Following directory are created under the current working directory.
#
#   <star1>_<obsdate1>
#       - pipeline_data: contain downloaded data.
#
#Updates:
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
    answer = raw_input(question)
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

def Obtain_scp_Files(user, server, password, download_data, target_direc, waittime=2., pflag=0):
    p = pexpect.spawn("scp %s@%s:%s %s" % (user, server, download_data, target_direc))
    time.sleep(waittime)
    if pflag == 0:
        p.expect("%s@%s's password:" % (user, server))
        p.sendline(password)
    p.interact()
    p.close()


def readstardirec(starlistfile, user, server, password, passflag):

    starlist = open(starlistfile,"r")
    starlines = starlist.readlines()
    starnum = len(starlines)
    starname = [starlines[i].split()[0] for i in range(starnum)]
    year = [starlines[i].split()[1] for i in range(starnum)]
    month_tmp = [starlines[i].split()[2] for i in range(starnum)]
    day_tmp = [starlines[i].split()[3] for i in range(starnum)]
    run = [starlines[i].split()[4] for i in range(starnum)]
    cut = [starlines[i].split()[5] for i in range(starnum)]
    
    day = []
    month = []
    for i in range(starnum):
        if len(day_tmp[i]) == 1:
            day.append("0%s" % day_tmp[i])
        elif len(day_tmp[i]) == 2:
            day.append(day_tmp[i])
        else:
            print "Error: illegal input for day."
            print "Input day \"%s\": Too long." % day_tmp[i]
            sys.exit()
        if len(month_tmp[i]) == 1:
            month.append("0%s" % month_tmp[i])
        elif len(month_tmp[i]) == 2:
            month.append(month_tmp[i])
        else:
            print "Error: illegal input for month."
            print "Input month \"%s\": Too long." % month_tmp[i]
            sys.exit()


    obsdate = ["%s_%s_%s" % (year[i], month[i], day[i]) for i in range(starnum)]
    
    direcpath = []
    for i in range(len(run)):
        if run[i] == "2nd":
            direcpath.append("/data/2013_02-03_reduction_data_revised")
        elif run[i] == "3rd":
            direcpath.append("/data/2013_11-12_reduction_data_revised")
        elif run[i] == "4th":
            direcpath.append("/data/2014_08-09_reduction_data")
        elif run[i] == "5th":
            direcpath.append("/media/WD_ext_1/5th_run_reduction_data")
        else:
            print "Error: illegal input for run number."
            print "Input run number \"%s\": Such Observational Run never exists." % run[i]
            sys.exit()

    pipelinedirec = ["%s/%s" % (direcpath[i],obsdate[i]) for i in range(starnum)]

    for i in range(starnum):
        tmplist_pipelinedirec = Obtain_ls_Result(user, server, password, pipelinedirec[i], "*_ver???", pflag=passflag)
        flag_pipelinedirec = 0
        tmplist_match = []
        for j in range(len(tmplist_pipelinedirec)):
            if tmplist_pipelinedirec[j].find("%s_pipeline_ver" % starname[i]) != -1:
                tmplist_match.append(j)
                flag_pipelinedirec += 1
        if flag_pipelinedirec == 0:
            print "Error: No reduced data directory is found."
            print "Cannot identify the pipeline directory for \"%s\"." % starname[i]
            sys.exit()
        else:
            pipelinedirec[i] += "/%s" % tmplist_pipelinedirec[max(tmplist_match)]


    framenumber = []
    for i in range(starnum):

        framenumber.append(len(Obtain_ls_Result(user, server, password, "%s" % pipelinedirec[i], "%s_NO*" % starname[i], pflag=passflag)))

    return starname, starnum, obsdate, cut, run, pipelinedirec, framenumber


def Data_Download(starname, starnum, target_direc, cut, pipelinedirec, framenumber, user, server, password, passflag):

    for i in range(starnum):
        for j in range(framenumber[i]):
            download_data = "%s/%s_NO%d/onedspec/object/flux/%s/*fits" % (pipelinedirec[i], starname[i], j+1, cut[i])
            if not os.path.exists(target_direc[i][j]):
                os.makedirs(target_direc[i][j])
            Obtain_scp_Files(user, server, password, download_data, target_direc[i][j], pflag=passflag)



def Data_Download_sum(starname, starnum, target_direc, cut, pipelinedirec, user, server, password, passflag):

    for i in range(starnum):
        download_data = "%s/%s_sum_data/onedspec/object/flux/%s/*fits" % (pipelinedirec[i], starname[i], cut[i])
        if not os.path.exists(target_direc[i]):
            os.makedirs(target_direc[i])
        Obtain_scp_Files(user, server, password, download_data, target_direc[i], pflag=passflag)




if __name__ == "__main__":

    filename = sys.argv[1:]

    user = "WINERED"
    server = "merlot.kyoto-su.ac.jp"

    print "Accessing %s:%s..." % (user, server)
    
    ### Password flag (added in ver.3)
    ### 0(default): you need to type password to access merlot server.
    ### non-zero number: you do not need to type password.
    
    passflag = 0
    
    ###
    if passflag == 0:
        password = command_question("Password: ")
    else:
        password = "hoge"
    
    starname, starnum, obsdate, cut, run, pipelinedirec, framenumber = readstardirec(filename[0], user, server, password, passflag)

    target_direc = [["%s_%s/pipeline_data%s/frame_NO%d" % (starname[i], obsdate[i].replace("_",""), pipelinedirec[i].split("pipeline")[-1], j+1) for j in range(framenumber[i])] for i in range(starnum)]

    target_direc_sum = ["%s_%s/pipeline_data%s/sum" % (starname[i], obsdate[i].replace("_",""), pipelinedirec[i].split("pipeline")[-1]) for i in range(starnum)]

    Data_Download(starname, starnum, target_direc, cut, pipelinedirec, framenumber, user, server, password, passflag)

    Data_Download_sum(starname, starnum, target_direc_sum, cut, pipelinedirec, user, server, password, passflag)


