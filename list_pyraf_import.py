import sys
import glob

if __name__ == "__main__":
    searchfile = '*.py'

    files = glob.glob(searchfile)
    files.sort()

    for i in files:
        rf = open(i, "r")
        rl = rf.readlines()
        rf.close()
        pyrafflag = False
        for j in range(len(rl)):
            if rl[j].find("from pyraf import iraf") != -1 and rl[j][0] == "f":
                pyrafflag = True
        if pyrafflag:
            print(i)
            for j in range(len(rl)):
                if rl[j].find("iraf.") != -1:
                    print(rl[j].rstrip("\n"))
            print()