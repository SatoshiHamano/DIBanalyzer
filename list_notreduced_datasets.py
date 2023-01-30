from open_mysql_project import openproject
import sys

if __name__ == "__main__":
    conn, cur = openproject()

    cur.execute("SELECT x.pipelineID from datareduction as x join object as y using(objectID) where y.type='OBJECT';")
    rows = cur.fetchall()
    ppid = [i[0] for i in rows]

    notreduced = []
    for i in range(len(ppid)):
        cur.execute("SELECT * from telluriccorrection where pipelineIDobj='%s';" % ppid[i])
        rows = cur.fetchall()
        if rows == []:
            notreduced.append(ppid[i])
            print(ppid[i])

    conn.close()