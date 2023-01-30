#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
import matplotlib.pyplot as plt

def sp2color(TYPE, SUBTYPE, CLASS):
    path = "./"
    csv = path + "Sp_B-V.csv"
    df = pd.read_csv(csv)
    return df.query("type==@TYPE & subtype==@SUBTYPE").iat[0,CLASS+1]

if __name__ == '__main__':

    a = sys.argv[1]
    b = sys.argv[2]
    c = int(sys.argv[3])
    d = float(sys.argv[4])

    col = sp2color(a,b,c)

    print(d - col)

    # type = ["O", "B", "A"]
    # typenum = [0, 10, 20]
    # subtype = [0.5 * i for i in range(20)]
    # lclass = [1,2,3,4,5]
    # colors = {1:"b", 2:"g", 3:"cyan", 4:"k", 5:"r"}
    #
    # plt.figure()
    # for t in range(len(type)):
    #     for s in subtype:
    #         for l in lclass:
    #             if sp2color(type[t],s,l) != "nan":
    #                 plt.scatter(typenum[t]+s, sp2color(type[t],s,l), color=colors[l])
    #
    # plt.savefig("Sp_B-V.png")
