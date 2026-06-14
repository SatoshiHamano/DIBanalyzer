import numpy
import matplotlib.pyplot as plt
import datetime
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == "__main__":
    rf = open("waveshift_analysis_4th_ver3.6.dat", "r")
    rl = rf.readlines()
    rf.close()

    ppid = numpy.array([i.split()[0] for i in rl])
    frame = numpy.array([i.split()[1] for i in rl])
    m = numpy.array([int(i.split()[2]) for i in rl])
    am = numpy.array([float(i.split()[3]) for i in rl])
    sc = numpy.array([float(i.split()[4]) for i in rl])
    sc_cor = numpy.array([float(i.split()[5]) for i in rl])

    fdate = [datetime.datetime(int(i[0:4]), int(i[5:7]), int(i[8:10])) for i in ppid]
    basedate = datetime.datetime(2014, 8, 12)
    datedif = [i - basedate for i in fdate]
    days = numpy.array([i.days for i in datedif])

    pp = PdfPages("waveshift_analysis_4th_2_ver3.6.pdf")

    plt.figure()
    mlist = list(range(42,62,1))

    amlist = []
    amliststd = []
    scclist = []
    sccliststd = []
    mlimlist = []

    for i in range(len(mlist)):
        req1 = (m == mlist[i]) & (frame == "sum")
        req2 = (am != 0.) & (sc != 0.)
        req = numpy.logical_and(req1, req2)
        if sum(req) != 0:
            plt.scatter(days[req], am[req], s=5, color="orange")
            plt.plot([min(days[req]), max(days[req])], [numpy.average(am[req]), numpy.average(am[req])])
            plt.title("m=%d" % mlist[i])
            plt.xlabel("days")
            plt.ylabel("a")
            plt.ylim(-0.004,0.004)
            plt.savefig(pp, format="pdf")
            plt.clf()

            amlist.append(numpy.average(am[req]))
            amliststd.append(numpy.std(am[req]))

            plt.scatter(days[req], sc[req], s=5., color="b")
            plt.plot([min(days[req]), max(days[req])], [numpy.average(sc[req]), numpy.average(sc[req])])
            plt.title("m=%d" % mlist[i])
            plt.xlabel("days")
            plt.ylabel("absolute shift from center")
            plt.ylim(-1.0, 1.0)
            plt.savefig(pp, format="pdf")
            plt.clf()

            plt.scatter(days[req], sc_cor[req], s=5., color="r")
            plt.plot([min(days[req]), max(days[req])], [numpy.average(sc_cor[req]), numpy.average(sc_cor[req])])
            plt.title("m=%d" % mlist[i])
            plt.xlabel("days")
            plt.ylabel("relative shift from center")
            plt.ylim(-0.3, 0.3)
            plt.savefig(pp, format="pdf")
            plt.clf()

            scclist.append(numpy.average(sc_cor[req]))
            sccliststd.append(numpy.std(sc_cor[req]))
            mlimlist.append(mlist[i])

    for i in range(len(mlimlist)):
        print(mlimlist[i], amlist[i], amliststd[i], scclist[i], sccliststd[i])

    plt.errorbar(mlimlist, amlist, yerr=amliststd)
    plt.xlabel("Echelle order")
    plt.ylabel("a")
    plt.savefig(pp, format="pdf")
    plt.clf()
    plt.errorbar(mlimlist, scclist, yerr=sccliststd)
    plt.xlabel("Echelle order")
    plt.ylabel("Relative shift")
    plt.savefig(pp, format="pdf")

    pp.close()