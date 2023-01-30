from open_mysql_project import openproject
import glob, sys, os

if __name__ == "__main__":
    conn, cur = openproject()

    checkid = sys.argv[1:]

    cur.execute("SELECT x.pipelineID,x.path from datareduction as x join object as y using(objectID) where y.type='STANDARD';")
    rows = cur.fetchall()
    ppid = [i[0] for i in rows]
    path = [i[1] for i in rows]

    ppidadv = []
    for i in range(len(ppid)):
        dirlist = glob.glob("%s*" % path[i])
        advdircand = []
        advdir = "NA"
        for j in dirlist:
            if os.path.isdir(j) and j.split(path[i])[1].find("_v") != -1:
                # print(ppid[i], j.split(path[i])[1])
                ppidadv.append(ppid[i])

    for i in checkid:
        if i in ppidadv:
            print(i, ": OK.")
        else:
            print(i, ": Yet.")

    conn.close()