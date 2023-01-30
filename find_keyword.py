import sys
import glob

if __name__ == "__main__":
    searchkeyword = sys.argv[1]
    searchfile = '*.py'

    files = glob.glob(searchfile)

    for i in files:
        rf = open(i, "r")
        rl = rf.readlines()
        rf.close()
        for j in range(len(rl)):
            if rl[j].find(searchkeyword) != -1:
                print(i, j+1, "line")
