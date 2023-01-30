from open_mysql_project import *
from Spec1Dtools import openspecfits
import numpy

if __name__ == "__main__":
    conn, cur = openproject()
    cur.execute("select * from telluricresult where telluricID='2014_10_17-CygOB2No.12_pipeline_ver3.5--2014_10_17-HR_196_pipeline_ver3.5-T1' and frame='sum';")
    rows = cur.fetchall()

    telid = [i[0] for i in rows]
    eorder = [i[1] for i in rows]
    scale = [i[2] for i in rows]
    scaleerr = [i[3] for i in rows]
    shift = [i[4] for i in rows]
    shifterr = [i[5] for i in rows]
    lammin = [i[6] for i in rows]
    lammax = [i[7] for i in rows]
    frame = [i[8] for i in rows]
    filepath = [i[9] for i in rows]

    telidt2 = "2014_10_17-CygOB2No.12_pipeline_ver3.5--2014_10_17-HR_196_pipeline_ver3.5-T2"
    filepatht2 = ["/Users/hamano/DIB_analysis/DIB_pipeline_dir/Schulte_12/2014_10_17-CygOB2No.12_pipeline_ver3.5/2014_10_17-CygOB2No.12_pipeline_ver3.5--2014_10_17-HR_196_pipeline_ver3.5-T2/sum/CygOB2No.12_sum_m%s_fsr1.30_VAC_tel.fits" % m for m in eorder]

    for i in range(20):
        spx,_,_,_,_ = openspecfits(filepatht2[i])
        lmin,lmax = numpy.amin(spx), numpy.amax(spx)
        valuetext = "'%s',%s,%s,%s,%s,%s,%.3f,%.3f,'%s','%s'" % (
                    telidt2, eorder[i], scale[i], scaleerr[i], shift[i], shifterr[i], lmin, lmax,
                    frame[i], filepatht2[i])
        cur.execute("INSERT IGNORE INTO telluricresult (telluricID,echelleorder,scale,scaleerr,shift,shifterr,lambdamin,lambdamax,frame,telluricfilepath) VALUES (%s);" % valuetext)

    conn.commit()

    conn.close()