import sys
from combine_MySQL import obtainCombinePath
from open_mysql_project import openproject

if __name__ == "__main__":
    conn, cur = openproject()

    cur.execute("select combineID, helio_velocity from DIBmeasurement where DIBID=37 and primaryflag=1;")
    rows = cur.fetchall()
    combineID = [i[0] for i in rows]
    heliov = [i[1] for i in rows]
    print(combineID)
    primaryFlag = ["-p" for i in range(len(combineID))]
    conn.close()

    # combineID = ["HD_149404_o246_c3T"]
    # primaryFlag = ["-p"]
    # heriov = [1.]

    # combineID = ["9_Gem_o124_c4T"]
    # primaryFlag = ["-p"]


    # combineID = ["HD_170740_o42_c3T"]
    # primaryFlag = ["-p"]

    # combineID = ["HD_169454_o41_c3T"]
    # primaryFlag = ["-p"]

    # combineID = ["HD_147889_o59_c3T", "b_Sco_o23_c3T", "HD_150898_o29_c3T", "HD_135591_o22_c3T",
    #              "ome_Sco_o24_c3T", "HD_151804_o30_c3T", "rho_Oph_D_o25_c3T", "chi_Oph_o26_c3T",
    #              "HD_152235_o63_c3T", "HD_185247_o45_c3T"]
    # primaryFlag = [" "] * 10
    # combineID = ["9_Gem_o124_c2T", "9_Gem_o124_c3T", "6_Cas_o90_c2T"]
    # primaryFlag = ["-p", " ", " "]

    # combineID = ["BD-16_4818_o8_c1T", "BD-16_4818_o8_c2T", "i_Sco_o60_c1T", "i_Sco_o60_c2T",
    #              "mu._Sgr_o65_c1T", "mu._Sgr_o65_c2T", "HD_148379_o27_c1T", "HD_148379_o27_c2T",
    #              "HD_152408_o31_c1T", "HD_152408_o31_c2T", "HD_154368_o32_c1T", "HD_154368_o32_c2T",
    #              "HD_155806_o33_c1T", "HD_155806_o33_c2T", "67_Oph_o34_c1T", "67_Oph_o34_c2T",
    #              "HD_164402_o35_c1T", "HD_164402_o35_c2T", "15_Sgr_o37_c1T", "15_Sgr_o37_c2T",
    #              "HD_168607_o39_c1T", "HD_168607_o39_c2T", "HD_168625_o40_c1T", "HD_168625_o40_c2T",
    #              "HD_169454_o41_c1T", "HD_169454_o41_c2T", "HD_170740_o42_c1T", "HD_170740_o42_c2T",
    #              "20_Aql_o43_c1T", "20_Aql_o43_c2T", "kap_Aql_o44_c1T", "kap_Aql_o44_c2T",
    #              "35_Aqr_o50_c1T", "35_Aqr_o50_c2T", "HD_214080_o52_c1T", "HD_214080_o52_c2T",
    #              "Hen_3-1250_o132_c1T", "Hen_3-1250_o132_c2T"
    #              ]
    # primaryFlag = [" ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p", " ", "-p",
    #                " ", "-p"
    #                ]

    # combineID = ["Wd1_W33_o144_c8T", "Kleinmann_o135_c5T", "Herschel_36_o133_c5T"]
    # primaryFlag = [" ", " ", " "]

    # combineID = ["Wd1_W33_o144_c3T", "Wd1_W33_o144_c4T", "Kleinmann_o135_c3T", "Kleinmann_o135_c4T", "Herschel_36_o133_c3T", "Herschel_36_o133_c4T"]
    # primaryFlag = [" ", "-p", " ", "-p", " ", "-p"]

    # # combineID = ["bet_Ori_o138_c3T", "bet_Ori_o138_c6R", "bet_Ori_o138_c7T"]
    # primaryFlag = [" ", " ", "-p"]

    # combineID = ["Schulte_12_o9_c1T", "Schulte_12_o9_c2T",
    #              "HD_183143_o66_c1T", "HD_183143_o66_c2T", "HD_20041_o75_c1T", "HD_20041_o75_c3T", "HD_12953_o55_c1T",
    #              "HD_21389_o86_c1T", "HD_50064_o128_c1T", "HD_50064_o128_c2T", "zet_Ori_A_o120_c1T",
    #              "chi02_Ori_o123_c1T", "chi02_Ori_o123_c2T", "6_Cas_o90_c1T"]
    # primaryFlag = [" ", "-p",
    #                "-p", " ", "-p", " ", "-p",
    #                "-p", "-p", " ", "-p",
    #                " ", "-p", "-p"]

    combpath = [obtainCombinePath(i) + "DIBs/" for i in combineID]

    DIBID = [686, 687, 688, 689]

    # DIBID = [54, 55, 25, 56, 57, 58, 43, 59, 60, 61, 62, 26, 63.44, 64, 65, 66, 27, 67, 68, 69, 70, 71, 72, 73, 74, 28,
    #          75, 29, 76, 30, 77, 78, 79, 80, 81, 82, 45, 46, 83, 84, 31, 47, 48, 85, 86, 32, 87, 50, 88, 51, 89, 90,
    #          91]  # selected
    # DIBID = [33, 34, 53, 24, 39, 54, 55, 25, 56, 57, 58, 40, 41, 42, 43, 59, 60, 61, 62, 26, 63, 44, 64, 65, 66, 35, 36,
    #          27, 67, 68, 69, 70, 71, 72, 73, 74, 28, 75, 29, 76, 37, 30, 77, 78, 79, 80, 81, 82, 45, 46, 83, 84, 31, 47,
    #          48, 85, 49, 86, 32, 87, 50, 88, 51, 89, 90, 52, 91, 38]  # all DIBs
    removeDIB = [24, 25, 59, 26, 63, 27, 29, 84]
    reduceDIB = []
    for i in DIBID:
        if not i in removeDIB:
            reduceDIB.append(i)

    wf = open(sys.argv[1], "w")

    for i in range(len(combineID)):
        for j in range(len(reduceDIB)):
            wf.write("python DIBcut.py {} -d {} -v {}\n".format(combineID[i], reduceDIB[j], heliov[i]))
            wf.write("python DIBanalysis.py -d {} -c {} {}\n".format(reduceDIB[j], combineID[i], primaryFlag[i]))

    wf.close()
