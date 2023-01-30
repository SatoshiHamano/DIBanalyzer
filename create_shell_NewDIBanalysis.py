import sys
from combine_MySQL import obtainCombinePath
from open_mysql_project import openproject

if __name__ == "__main__":
    conn, cur = openproject()

    #objectID = [138,50,75,23,22,24,33,43,30,123,31,42,124,90,246,32,27,63,45,59,41,66,40,39,14,10,12,13,15,9,132,144]
    objectID = [50,23,22,24,33,43,30,31,42,246,27,32,63,45,59,41,40,39,132,144]
    objIDstr = ""
    for i in objectID:
        objIDstr += str(i) + ","
    objIDstr = objIDstr.rstrip(",")
    #NewDIB = [679, 680, 681, 682]
    NewDIB = [683,684,685]

    cur.execute("select combineID,helio_velocity from DIBmeasurement join combinesummary using (combineID) where DIBID = 35 and objectID IN (%s) and primaryflag =1;" % objIDstr)
    rows = cur.fetchall()
    combineID = [i[0] for i in rows]
    helio_v = [i[1] for i in rows]
    conn.close()

    wf = open(sys.argv[1], "w")

    for i in range(len(combineID)):
        for j in range(len(NewDIB)):
            wf.write("python DIBcut.py %s -d %d -v %.2f\n" % (combineID[i], NewDIB[j], helio_v[i]))
            wf.write("python DIBanalysis.py -d %d -c %s -p\n" % (NewDIB[j], combineID[i]))

    wf.close()

