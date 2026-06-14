import sys,os,datetime
import numpy,math

filename = sys.argv[1:]

#Description:
#   This "Waveshift_check.py" script was made by Satoshi Hamano in 2016/05/09.
#
#   This script enables you to check the quality of the Wavelength_fit measurement.
#
#   Checked points are the following:
#
#       1) Standard deviation of wavelength residuals after correcting the wavelength
#          shifts is lower than set threshold (default: 0.2 angstrom).
#       2) Number of clipped points is lower than set threshold (default: 5 points).
#       3) The relative differences of wavelengths of telluric lines become less by
#          correcting the wavelength shifts.
#          (This check is done for all pairs of frames.)
#
#Usage:
#
#   $ python Waveshift_check.py <output_list> <status_file>
#
#   output_list: the list of output files, which are created by Waveshift_fit script.
#   status_file: the file, in which the status of the waveshift_fit script is stored.
#
#Output:
#
#   The result of the three checks and the status of the Waveshift_fit script is written in status_file.
#
#Updates:
#
#
##



def Waveshift_check(outputlist,status_output,sd_pix=1.0,lambda_per_pix=0.2,clipsafe=5):

    numlist = len(outputlist)

    [inc, offset] = [[], []]
    [waveshift, centerlam, centerpix, order, lambdac, waveshift_cl, centerlam_cl, order_cl, lambdac_cl] = [ [[] for j in range(numlist)] for i in range(9)]

    for i in range(numlist):
        rf = open(outputlist[i],"r")
        rl = rf.readlines()
        rf.close()

        inc.append(float(rl[2].split()[5]))
        offset.append(float(rl[2].split()[9]))

        for j in range(4,len(rl)):
            if rl[j].find("begin") != -1:
                beginline = j
    
        for j in range(4,beginline-1):
            waveshift_cl[i].append(float(rl[j].split()[2]))
            centerlam_cl[i].append(float(rl[j].split()[1]))
            order_cl[i].append(float(rl[j].split()[3]))
            lambdac_cl[i].append(float(rl[j].split()[4]))
    
        for j in range(beginline+2,len(rl)):
            waveshift[i].append(float(rl[j].split()[1]))
            centerlam[i].append(float(rl[j].split()[0]))
            centerpix[i].append(float(rl[j].split()[2]))
            order[i].append(float(rl[j].split()[3]))
            lambdac[i].append(float(rl[j].split()[4]))

    list_of_list = [waveshift, centerlam, centerpix, order, lambdac, waveshift_cl, centerlam_cl, order_cl, lambdac_cl]
    [waveshift, centerlam, centerpix, order, lambdac, waveshift_cl, centerlam_cl, order_cl, lambdac_cl] = [ [numpy.array(list_of_list[i][j]) for j in range(len(list_of_list[i]))] for i in range(len(list_of_list)) ]

    wf = open(status_output,"w")

    flag_sd_eachframe = [0 for i in range(numlist)]
    sd_eachframe = []

    for i in range(numlist):
        sd_eachframe.append(numpy.std(waveshift[i] - inc[i] * lambdac[i] - offset[i]))
        if sd_eachframe[i] > sd_pix * lambda_per_pix:
            flag_sd_eachframe[i] = 1

    flag_clip = [0 for i in range(numlist)]
    num_clip = []
    for i in range(numlist):
        num_clip.append(len(waveshift_cl[i]))
        if num_clip[i] >= clipsafe:
            flag_clip[i] = 1

    flag_relative_shift = [[0 for i in range(numlist)] for j in range(numlist)]
    shift_sd_before = [[] for i in range(numlist)]
    shift_sd_after = [[] for i in range(numlist)]
    for i in range(numlist):
        for j in range(numlist):
            if i != j:
                tmp_shift_before = []
                tmp_shift_after = []
                for k in range(len(waveshift[i])):
                    for n in range(len(waveshift[j])):
                        if [centerlam[i][k],order[i][k]] == [centerlam[j][n],order[j][n]]:
                            tmp_shift_before.append(waveshift[i][k] - waveshift[j][n])
                            tmp_shift_after.append((waveshift[i][k] - inc[i] * lambdac[i][k] - offset[i]) - (waveshift[j][n] - inc[j] * lambdac[j][n] - offset[j]))
                shift_sd_before[i].append(numpy.median(numpy.array(tmp_shift_before)))
                shift_sd_after[i].append(numpy.median(numpy.array(tmp_shift_after)))
        
            if i == j:
                shift_sd_before[i].append(0)
                shift_sd_after[i].append(0)

    for i in range(numlist):
        for j in range(numlist):
            if math.fabs(shift_sd_before[i][j]) < math.fabs(shift_sd_after[i][j]) and shift_sd_after[i][j] > numpy.std(numpy.array(shift_sd_after))*5.:
                flag_relative_shift[i][j] = 1

    wf.write("Time: %s\n" % datetime.datetime.today())
    wf.write("Checked Waveshift outputs:\n")
    for i in range(numlist):
        wf.write("\t%s\n" % outputlist[i])
    flags_all = 0
    for i in range(numlist):
        flags_all += flag_sd_eachframe[i]
        flags_all += flag_clip[i]
        for j in range(numlist):
            flags_all += flag_relative_shift[i][j]
    if flags_all == 0:
        wf.write("\nStatus: SUCCESS\n")
        check_status = "success"
    else:
        wf.write("\nStatus: FAILURE\n")
        check_status = "failure"

    wf.write("\nS.D. of each frame (<%.2f angstrom: OK):\n" % (sd_pix * lambda_per_pix))
    for i in range(numlist):
        if flag_sd_eachframe[i] == 0:
            wf.write("\tFrame %d: %.4f --- %s\n" % (i+1, sd_eachframe[i], "OK"))
        else:
            wf.write("\tFrame %d: %.4f --- %s\n" % (i+1, sd_eachframe[i], "NG"))

    wf.write("\nNumber of clipped lines (<%d lines: OK):\n" % clipsafe)
    for i in range(numlist):
        if flag_clip[i] == 0:
            wf.write("\tFrame %d: %d --- %s\n" % (i+1, num_clip[i], "OK"))
        else:
            wf.write("\tFrame %d: %d --- %s\n" % (i+1, num_clip[i], "NG"))

    wf.write("\nImprovement of wavelength shifts (S.D.: %.2e)\n" % numpy.std(numpy.array(shift_sd_after)))
    for i in range(numlist):
        for j in range(numlist):
            if i != j:
                if flag_relative_shift[i][j] == 0:
                    wf.write("\tFrame %d - %d: %.2e -> %.2e --- %s\n" % (i+1, j+1, shift_sd_before[i][j], shift_sd_after[i][j], "OK"))
                else:
                    wf.write("\tFrame %d - %d: %.2e -> %.2e --- %s\n" % (i+1, j+1, shift_sd_before[i][j], shift_sd_after[i][j], "NG"))

    wf.write("\n\nFin.")

    wf.close()

    return check_status

if __name__ == "__main__":

    filename = sys.argv[1:]
    rf = open(filename[0],"r")
    rl = rf.readlines()
    rf.close()
    outputlist = [rl[i].split()[0] for i in range(len(rl))]
    print "Status of Waveshift: %s" % Waveshift_check(outputlist, filename[1])


