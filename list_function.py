import sys
import glob

if __name__ == "__main__":
    searchfile = '*.py'

    files = glob.glob(searchfile)
    files.sort()

    for i in files:
        print(i)
        rf = open(i, "r")
        rl = rf.readlines()
        rf.close()
        functionbody = []
        functionline = []
        returnstate = []
        for j in range(len(rl)):
            if rl[j].find("def ") != -1 and rl[j][0] == "d":
                functionbody.append(rl[j].lstrip("def ").rstrip("\n"))
                functionline.append(j)
        functionline.append(len(rl))
        for k in range(len(functionbody)):
            returnlist = []
            for j in range(functionline[k], functionline[k+1]):
                if rl[j].find("return ") != -1:
                    returnlist.append(rl[j].lstrip("return ").rstrip("\n"))
            if returnlist == []:
                print("\t", functionbody[k])
            else:
                print("\t", functionbody[k], "\n\t  -->", returnlist)

        print()