# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import numpy as np

root_path = os.path.dirname(os.path.abspath(__file__))

atomfile = root_path + '/linelist.dat'

linelist = np.loadtxt(atomfile, dtype=[('trans', 'S13'),
                                       ('ion', 'S6'),
                                       ('l0', 'f4'),
                                       ('f', 'f8'),
                                       ('gam', 'f8'),
                                       ('mass', 'f4')]
                      )

molecname = set(["C2", "12C13C", "CN"])  # , "12C15N", "13C14N"])
molecID = {"C2": "C2J",
           "12C13C": "C2iJ",
           "CN": "CN"
           }

bandID = {"C2": ["AX00", "AX10", "AX20"],
          "12C13C": ["AX00", "AX10"],
          "CN": ["AX00", "AX10"]
          }

bandname_short = {"C2":
                      {"AX00": "AX (0,0)",
                       "AX10": "AX (1,0)",
                       "AX20": "AX (2,0)"
                       },
                  "12C13C":
                      {"AX00": "AX (0,0)",
                       "AX10": "AX (1,0)"
                       },
                  "CN":
                      {"AX00": "AX (0,0)",
                       "AX10": "AX (0,0)"
                       }
                  }

bandname_long = {"C2":
                     {"AX00": "A^1\Pi_u-X^1\Sigma_g^+ (0,0)",
                      "AX10": "A^1\Pi_u-X^1\Sigma_g^+ (1,0)",
                      "AX20": "A^1\Pi_u-X^1\Sigma_g^+ (2,0)"
                      },
                 "12C13C":
                     {"AX00": "A^1\Pi_u-X^1\Sigma_g^+ (0,0)",
                      "AX10": "A^1\Pi_u-X^1\Sigma_g^+ (1,0)"
                      },
                 "CN":
                     {"AX00": "A^2\Pi_u-X^2\Sigma^+ (0,0)",
                      "AX10": "A^2\Pi_u-X^2\Sigma^+ (0,0)"
                      }
                 }

band_lamrange = {"C2":
                     {"AX00": [12050., 12500.],
                      "AX10": [10120., 10600.],
                      "AX20": [8700., 9100.]
                      },
                 "12C13C":
                     {"AX00": [12050., 12500.],
                      "AX10": [10120., 10600.]
                      },
                 "CN":
                     {"AX00": [10900., 11100.],
                      "AX10": [9100., 9200.]
                      }
                 }

firsttrans_PQR = {"C2":
                      {"AX00": ["C2J2_12105.5", "C2J2_12096.0", "C2J0_12089.6"],
                       "AX10": ["C2J2_10157.7", "C2J2_10151.1", "C2J0_10146.5"],
                       "AX20": ["C2J2_8768.4", "C2J2_8763.6", "C2J0_8760.1"]
                       },
                  "12C13C":
                      {"AX00": ["C2iJ1_12096.0", "C2iJ0_12090.9", "C2iJ0_12086.4"],
                       "AX10": ["C2iJ1_10182.4", "C2iJ0_10178.8", "C2iJ0_10175.7"]
                       },
                  }


def vactoair(lambda_vac):
    s = pow(10, 4) / lambda_vac
    n = 1.0 + 5.792105e-2 / (238.0185 - s ** 2) + 1.67917e-3 / (57.362 - s ** 2)

    return lambda_vac / n

def airtovac(lambda_air):
    s = pow(10, 4) / lambda_air
    n = 1.0 + 5.792105e-2 / (238.0185 - s ** 2) + 1.67917e-3 / (57.362 - s ** 2)

    return lambda_air * n


class lineparam:
    def __init__(self, molec, band, linetag, branch, rotJ, lamvac, fvalue):
        self.molec = molec
        self.band = band
        self.linetag = linetag
        self.branch = branch
        self.rotJ = rotJ
        self.transition = "%s(%d)" % (self.branch, self.rotJ)
        self.lamvac = lamvac
        self.wavenumber = 1.e+8 / self.lamvac
        self.lamair = vactoair(self.lamvac)
        self.fvalue = fvalue
        self.level = self.linetag.split("_")[0]


class molecband:
    def __init__(self, molec, band):
        self.molec = molec
        self.band = band
        self.bandshort = bandname_short[self.molec][self.band]
        self.bandlong = bandname_long[self.molec][self.band]
        self.lines = []
        if self.molec == "C2" or molec == "12C13C":
            branch = ["P", "Q", "R"]
            branchid = -1
        else:
            branch = [" "]
            branchid = 0

        for line in linelist:
            transstr = line["trans"].astype(str)
            if transstr.find(molecID[self.molec]) != -1:
                if band_lamrange[self.molec][self.band][0] < line["l0"] < band_lamrange[self.molec][self.band][1]:
                    if transstr == firsttrans_PQR[self.molec][self.band][0]:
                        branchid = 0
                    elif transstr == firsttrans_PQR[self.molec][self.band][1]:
                        branchid = 1
                    elif transstr == firsttrans_PQR[self.molec][self.band][2]:
                        branchid = 2

                    self.lines.append(lineparam(self.molec, self.band, transstr, branch[branchid],
                                                int(line["ion"].astype(str).replace(molecID[self.molec], "")),
                                                line["l0"], line["f"]))

    def branchlist(self, branch):
        branchlist = []
        for line in self.lines:
            if line.branch == branch:
                branchlist.append(line)

        return branchlist

    def minimumlambda(self, unit="air"):
        if unit == "air":
            minlam = self.lines[0].lamair
            for line in self.lines:
                minlam = min(minlam, line.lamair)
        elif unit == "vac":
            minlam = self.lines[0].lamvac
            for line in self.lines:
                minlam = min(minlam, line.lamvac)
        else:
            print("Unit is vac or air! Not %s!" % unit)
            minlam = "INDEF"
        return minlam

    def maximumlambda(self, unit="air"):
        if unit == "air":
            maxlam = self.lines[0].lamair
            for line in self.lines:
                maxlam = max(maxlam, line.lamair)
        elif unit == "vac":
            maxlam = self.lines[0].lamvac
            for line in self.lines:
                maxlam = max(maxlam, line.lamvac)
        else:
            print("Unit is vac or air! Not %s!" % unit)
            maxlam = "INDEF"
        return maxlam


class moleclinelist:
    def __init__(self, molec):
        self.molec = molec
        self.bandlist = {}
        for i in bandID[self.molec]:
            self.bandlist[i] = molecband(self.molec, i)
