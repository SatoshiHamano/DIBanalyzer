import sys
import mysql.connector
from urllib.parse import urlparse
import glob
from Spec1Dtools import openspecfits

def telluric_log_reader(tellog):
    rf = open(tellog, "r")
    rl = rf.readlines()
    rf.close()

    outputfile = "NA"
    shift = "NA"
    scale = "NA"
    for i in rl:
        if i.find("Output:") != -1:
            outputfile = i.split()[1]
        if i.find("Tweak:") != -1:
            shift = float(i.split()[3].rstrip(","))
            scale = float(i.split()[6].rstrip(","))

    if "NA" in [outputfile, shift, scale]:
        print(tellog + " cannot be read!!!")
        print()
        for i in rl:
            print(rl)
        sys.exit()
    else:
        return outputfile, shift, scale


if __name__ == "__main__":
    telluricID = sys.argv[1]

    urlsql = urlparse('mysql://root:kwmjbqb9py@localhost:3306/DIBproject')

    conn = mysql.connector.connect(
        host=urlsql.hostname or 'localhost',
        port=urlsql.port or 3306,
        user=urlsql.username or 'root',
        password=urlsql.password or 'kwmjbqb9py',
        database=urlsql.path[1:],
    )

    cur = conn.cursor()

    cur.execute(
        "select autoflag, telluricpath from telluriccorrection where telluricID = '%s';" % telluricID)
    rows = cur.fetchall()
    if rows == []:
        print("%s is not found." % telluricID)
        sys.exit()
    else:
        autoflag = rows[0][0]
        telluricpath = rows[0][1]

    if autoflag == 1:
        print("%s is produced by auto mode." % telluricID)
        sys.exit()

    cur.execute(
        "select echelleorder, frame from telluricresult where telluricID = '%s';" % telluricID)
    rows = cur.fetchall()
    if rows == []:
        firstflag = True
    else:
        firstflag = False
        eorder_pre = [int(i[0]) for i in rows]
        frame_pre = [int(i[1]) for i in rows]

    frames = glob.glob(telluricpath + "*")
    frames.sort()

    for i in frames:
        tellog = glob.glob(i + "/" + "telluric_log*.txt")
        if tellog == []:
            continue
        else:
            for j in tellog:
                m = int(j.split("/")[-1].split("telluric_log_m")[1].rstrip(".txt"))
                if not firstflag:
                    if m in eorder_pre and i.split("/")[-1] in frame_pre:
                        continue
                else:
                    outputfile, shift, scale = telluric_log_reader(j)
                    spx, spy, _, _, _ = openspecfits(i + "/" + outputfile)
                    minlam = min(spx)
                    maxlam = max(spx)

                    valuetext = "'%s',%d,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,'%s','%s'" % (
                        telluricID, m, scale, -1, shift, -1, minlam, maxlam,
                        i.split("/")[-1], i + "/" + outputfile)
                    cur.execute(
                        "INSERT IGNORE INTO telluricresult (telluricID,echelleorder,scale,scaleerr,shift,shifterr,lambdamin,lambdamax,frame,telluricfilepath) VALUES (%s);" % valuetext)

                    print("The result of frame=%s and m=%d was recorded!" % (i.split("/")[-1], m))

    conn.commit()
    conn.close()