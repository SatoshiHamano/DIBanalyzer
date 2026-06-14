from matplotlib.backends.backend_pdf import PdfPages
import sys
import mysql.connector
from urllib.parse import urlparse
import numpy
from waveshift_measure import read_waveshift
import matplotlib.pyplot as plt


if __name__ == "__main__":
    urlsql = urlparse('mysql://root:kwmjbqb9py@localhost:3306/DIBproject')

    conn = mysql.connector.connect(
        host=urlsql.hostname or 'localhost',
        port=urlsql.port or 3306,
        user=urlsql.username or 'root',
        password=urlsql.password or 'kwmjbqb9py',
        database=urlsql.path[1:],
    )

    cur = conn.cursor()

    fsr = "fsr1.30"
    vacorair = "VAC"
    nite = 5
    lowsig = 2

    cur.execute(
        "select pipelineID, FrameNum, totalSNR, mode, obsdate, path from datareduction where obsdate between '2014-08-01 00:00:00' and '2014-11-01 00:00:00';")
    rows = cur.fetchall()
    ppid = [i[0] for i in rows]
    fnum = [int(i[1]) for i in rows]
    snr = [float(i[2]) for i in rows]
    mode = [i[3] for i in rows]
    obsdate = [i[4] for i in rows]
    path = [i[5] for i in rows]

    rf = open("waveshift_parameter_4thrun.dat", "r")
    rl = rf.readlines()
    rf.close()

    orderpa = [int(i.split()[0]) for i in rl]
    am = [float(i.split()[1]) for i in rl]
    shiftpa = [float(i.split()[3]) for i in rl]

    pp = PdfPages("waveshift_analysis_4th_3.pdf")
    plt.figure()

    for i in range(len(ppid)):
        frame = ["sum"]  # ["NO%d" % (j+1) for j in range(fnum[i])] + ["sum"]
        for j in range(len(frame)):
            wsfile = "%s%s%s/%s/%s_%s_%s_norm_%s.txt" % (
            path[i], "waveshift_measure/", frame[j], fsr, ppid[i], fsr, vacorair, frame[j])
            order, split, wavc, shift, r_edge, r_min, r_dif, r_ratio = read_waveshift(wsfile)
            m = list(set(order))
            m.sort()

            wavc_m = {}
            for k in m:
                wavc_m[k] = numpy.average(wavc[order == k])
            req1 = numpy.array([r_edge[k] > 0.5 for k in range(len(order))])

            clip1 = numpy.array([True for k in range(len(order))])
            for k in range(nite):
                r_ratio_av = numpy.average(r_ratio[numpy.logical_not(req1) & clip1])
                r_ratio_std = numpy.std(r_ratio[numpy.logical_not(req1) & clip1])
                clip1[(r_ratio - r_ratio_av) < - lowsig * r_ratio_std] = False

            req = numpy.logical_not(numpy.logical_or(req1, numpy.logical_not(clip1)))

            plt.scatter(wavc[req], shift[req], color="b")
            for k in range(len(m)):
                plt.scatter(wavc[numpy.logical_and(req, order == m[k])],
                            shift[numpy.logical_and(req, order == m[k])] - am[k] * (
                                        wavc[numpy.logical_and(req, order == m[k])] - wavc_m[m[k]]) - shiftpa[k],
                            color="r")
            plt.title(ppid[i])
            plt.savefig(pp, format="pdf")
            plt.clf()

    pp.close()
