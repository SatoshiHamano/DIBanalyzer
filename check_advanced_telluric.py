from telluric_auto import GetListAdvancedFits
from open_mysql_project import openproject
import glob, sys
from spectra_plotter import MultiSpecPlotter, MultiSpecAxis
from matplotlib.backends.backend_pdf import PdfPages
from Spec1Dtools import openspecfits
import matplotlib.pyplot as plt
from atran_model_npz import modelnpzopen

if __name__ == "__main__":
    # ppid_tel -> path_tel -> GetAd
    conn, cur = openproject()

    # rf = open(sys.argv[1], "r")
    # rl = rf.readlines()
    # rf.close()

    # ppidtel = [i.split()[0] for i in rl]
    ppidtel = sys.argv[1]
    ppidobj = sys.argv[2]
    righttexts = ["   Pipeline spectra", "   Intrinsic line removed"]
    colors = ["r", "b"]
    telx, tely, _ = modelnpzopen()
    obsdate = ppidtel[0:10]

    # WIDE mode only
    modeltxt = {}
    for i in range(42, 62):
        modeltxt[i] = "A0V_template/template_m%d_cut5.txt" % i

    cur.execute(
        "select pipelinever, path from datareduction where pipelineID='%s';" %
        ppidtel)
    rows = cur.fetchall()
    if rows == []:
        print("%s is not found in the database table 'DIBproject.datareduction'." % ppidtel)
        sys.exit()

    pipelineve = rows[0][0]
    path = rows[0][1]

    spf_bef = glob.glob(path + "*_sum/VAC_norm/fsr1.30/*_norm.fits")
    spf_bef.sort()
    m_bef = [int(spf.split("_fsr")[-2].split("_m")[-1]) for spf in spf_bef]

    spf_aft = GetListAdvancedFits(path, "fsr?.??")
    spf_aft.sort()
    m_aft = [int(spf.split("_m")[-1].rstrip("fits").rstrip(".")) for spf in spf_aft]

    cur.execute("select telluricnumber from telluriccorrection where pipelineIDobj = '%s' and pipelineIDtel = '%s' and advanced=1;"  % (ppidobj, ppidtel))
    rows = cur.fetchall()
    if rows == []:
        print("No advanced telluric-corrected spectra.")
        sys.exit()
    else:
        advnumlist = [i[0] for i in rows]

    if len(advnumlist) == 1:
        advnum = advnumlist[0]
    else:
        print("List of advanced telluric number:", advnumlist)
        advnum = int(input("Enter telluric number:"))

    cur.execute(
        "select telluricnumber from telluriccorrection where pipelineIDobj = '%s' and pipelineIDtel = '%s' and advanced=0;" % (
        ppidobj, ppidtel))
    rows = cur.fetchall()
    if rows == []:
        print("No simply telluric-corrected spectra.")
        sys.exit()
    else:
        simplenumlist = [i[0] for i in rows]

    if len(simplenumlist) == 1:
        simnum = simplenumlist[0]
    else:
        print("List of simple telluric number:", simplenumlist)
        simnum = int(input("Enter telluric number:"))

    frame = "sum"

    advdir = spf_aft[0].split("/telluric/")[0].split("/")[-1]
    pp = PdfPages(path + "Check_%s.pdf" % advdir)
    plt.figure(figsize=(20, 20))

    for j in range(len(m_aft)):
        for k in range(len(m_bef)):
            if m_aft[j] == m_bef[k]:
                if m_aft[j] >= 42 and m_aft[j] <= 61:
                    modelkey = modeltxt[m_aft[j]]
                else:
                    modelkey = "INDEF"


                spx_b, spy_b, _, _, _ = openspecfits(spf_bef[k])
                spx_a, spy_a, _, _, _ = openspecfits(spf_aft[j])
                spxlist, spylist = [spx_b, spx_a], [spy_b, spy_a]

                cur.execute("select y.telluricfilepath from telluriccorrection as x join telluricresult as y using(telluricID) where x.pipelineIDobj = '%s' and x.pipelineIDtel = '%s' and y.frame='%s' and x.telluricNumber=%d and y.echelleorder=%d;" % (ppidobj, ppidtel, frame, advnum, m_aft[j]))
                rows = cur.fetchall()
                spx_adv, spy_adv, _, _, _ = openspecfits(rows[0][0])
                cur.execute("select y.telluricfilepath from telluriccorrection as x join telluricresult as y using(telluricID) where x.pipelineIDobj = '%s' and x.pipelineIDtel = '%s' and y.frame='%s' and x.telluricNumber=%d and y.echelleorder=%d;" % (ppidobj, ppidtel, frame, simnum, m_aft[j]))
                rows = cur.fetchall()
                spx_sim, spy_sim, _, _, _ = openspecfits(rows[0][0])
                spxlist_obj, spylist_obj = [spx_sim, spx_adv], [spy_sim, spy_adv]

                axtel = plt.axes([0.08, 0.1, 0.9, 0.4])

                MultiSpecAxis(axtel, spxlist, spylist, colors, xaxis_label="Wavelength", sptype="A", zerooffset=True, vacorair="VAC", C2plot=False, obsdates=obsdate, order=m_aft[j], model=modelkey, )

                axobj = plt.axes([0.08, 0.55, 0.9, 0.4])

                MultiSpecAxis(axobj, spxlist_obj, spylist_obj, colors, sptype="OB", zerooffset=False, vacorair="VAC", C2plot=True, obsdates=obsdate, order=m_aft[j])
                plt.title("%s (%s, %s, m=%d)" % (advdir, ppidtel, ppidobj, m_aft[j]))
                plt.savefig(pp, format="pdf")
                plt.clf()

                # MultiSpecPlotter(spxlist, spylist, pp, colors, righttexts, "%s (m=%d)" % (ppidtel[i], m_aft[j]),
                #                  sptype="A", telflag=True, vacorair="VAC", C2plot=False, zerooffset=True,
                #                  model=modelkey, telx=telx, tely=tely)

    pp.close()

    print(path + "Check_%s.pdf" % advdir)