import math

rf = open("WODB_std_all.txt", "r")
rl = rf.readlines()
rf.close()

pf = open("pipeline_info_table_4.dat", "r")
pl = pf.readlines()
pf.close()

wf = open("WODB_all_csv_revised.txt", "a")

airmass_pp = {}
object_pp = {}
idlist = []
NAlist = []

for i in pl:
    pl1 = i.split(" || ")
    frame = pl1[2]
    try:
        airmass_start = float(pl1[14])
        airmass_end = float(pl1[15])
    except:
        NAlist.append(frame)
        continue
    if not frame in idlist:
        if math.fabs(airmass_start - 1.0) > 1.5 or math.fabs(airmass_end - 1.0) > 1.5:
            print(frame + ": airmass (%.2f-%.2f) is out of range." % (airmass_start, airmass_end))
        else:
            idlist.append(frame)
            airmass_pp[frame] = [airmass_start, airmass_end, (airmass_start + airmass_end) / 2.]

for i in range(len(rl)):
    rl1 = rl[i].split(" || ")
    [frame, object, theme, type, obsdate, instmode, slit, exptime, position, airmass, flag, memo, a] = rl1
    try:
        airmass_v = float(airmass)
    except:
        airmass_v = 0.
    if airmass_v == 0.:
        if frame in idlist:
            airmass = str(airmass_pp[frame][2])
        else:
            airmass = "0.0"
    wf.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
        frame, object, theme, type, obsdate, instmode, slit, exptime, position, airmass, flag, memo.replace(",", " ")))

wf.close()
