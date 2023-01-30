import sys

if __name__ == "__main__":
    rf = open(sys.argv[1], "r")
    rl = rf.readlines()
    rf.close()

    rf2 = open("temporaly_files/advancedTelluric_20210107.txt", "r")
    rl2 = rf2.readlines()
    rf2.close()

    wp = "waveshift_parameter/waveshift_parameter_blank_wide.dat"
    #wp = sys.argv[3]

    advanced = [i.split()[0] for i in rl2]

    targetID = []
    telluricID = []
    advancedFlag = []
    for i in rl:
        rl1 = i.split()
        if len(rl1) > 1:
            targetID.append(rl1[0])
            telluricID.append(rl1[1])
            if rl1[1] in advanced:
                advancedFlag.append(True)
            else:
                advancedFlag.append(False)

    print(advancedFlag)

    wf = open(sys.argv[2], "w")

    for i in range(len(targetID)):
        wf.write("python telluric_auto.py %s %s\n" % (targetID[i], telluricID[i]))
        wf.write("python waveshift_correct.py %s--%s-T1 %s\n" % (targetID[i], telluricID[i][11:], wp))
        if advancedFlag[i]:
            wf.write("python telluric_auto.py %s %s -a\n" % (targetID[i], telluricID[i]))
            wf.write("python waveshift_correct.py %s--%s-T2 %s\n" % (targetID[i], telluricID[i][11:], wp))

    wf.close()