import sys
import mysql.connector
from urllib.parse import urlparse
import numpy
from waveshift_measure import read_waveshift

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
        "select pipelineID, FrameNum, totalSNR, mode, obsdate, path, pipelinever from datareduction where obsdate between '2016-01-26 00:00:00' and '2016-03-26 00:00:00';")
    rows = cur.fetchall()
    ppid = [i[0] for i in rows]
    fnum = [int(i[1]) for i in rows]
    snr = [float(i[2]) for i in rows]
    mode = [i[3] for i in rows]
    obsdate = [i[4] for i in rows]
    path = [i[5] for i in rows]
    pver = [i[6] for i in rows]

    wf = open("waveshift_analysis_5th20160131_ver3.5.dat", "w")

    for i in range(len(ppid)):
        if "3.5" in pver[i]:
            frame = ["NO%d" % (j+1) for j in range(fnum[i])] + ["sum"]
            for j in range(len(frame)):
                wsfile = "%s%s%s/%s/%s_%s_%s_norm_%s.txt" % (path[i], "waveshift_measure/", frame[j], fsr, ppid[i], fsr, vacorair, frame[j])
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

                am = []
                sc = []
                for k in range(len(m)):
                    order_req = order[req]
                    wavc_req = wavc[req]
                    shift_req = shift[req]
                    if len(order_req[order_req == m[k]]) > 2:
                        a, b = numpy.polyfit(wavc_req[order_req == m[k]], shift_req[order_req == m[k]], 1)
                        am.append(a)
                        sc.append(b + a * wavc_m[m[k]])
                    else:
                        am.append(0.)
                        sc.append(0.)
                am = numpy.array(am)
                sc = numpy.array(sc)
                sc_ave = numpy.average(sc[(sc != 0.) & (am != 0.)])
                sc_cor = sc - sc_ave

                for k in range(len(m)):
                    wf.write("%s\t%s\t%d\t%.4e\t%.4e\t%.4e\n" % (ppid[i], frame[j], m[k], am[k], sc[k], sc_cor[k]))

    wf.close()
    conn.close()