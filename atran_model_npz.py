# -*- coding:utf-8 -*-

import numpy, math
import matplotlib.pyplot as plt
import sys
from vac2air_spec import vac2air
import scipy.constants

def modelopen(modeldat):
    rf1 = open(modeldat, "r")
    rl = rf1.readlines()
    rf1.close()

    id = numpy.array([int(i.split()[0]) for i in rl])
    wav = numpy.array([float(i.split()[1]) * 1.e+4 for i in rl])
    flux = numpy.array([float(i.split()[2]) for i in rl])

    return id, wav, flux

def modelnpzopen(vacorair="VAC"):
    modelnpz = "atran.smo.11513_R28000_0ft.npz"

    atran_npz = numpy.load(modelnpz)
    wav = atran_npz["wav"]
    flux = atran_npz["flux"]
    id = atran_npz["id"]
    if vacorair == "AIR":
        wav_air = vac2air(wav)
        return wav_air, flux, id
    else:
        return wav, flux, id


# def modelnpzhelioSave(modeldirpath, vhelio):
#     telmodelhelio = "telluric_model_AIR_helio.npz"
#     vel_light = scipy.constants.c * 1.e-3
#
#     wav, flux, id = modelnpzopen(vacorair="AIR")
#     wav_shifted = wav * (1. + vhelio / vel_light)
#
#     numpy.savez(modeldirpath + telmodelhelio, id=id, wav=wav_shifted, flux=flux)


def modelnpzopenHelio(vhelio, vacorair="AIR"):
    wav, flux, id = modelnpzopen(vacorair=vacorair)
    vel_light = scipy.constants.c * 1.e-3

    return wav * (1. + vhelio / vel_light), flux, id


# def modelnpzhelioopen(modeldirpath):
#     telmodelhelio = "telluric_model_AIR_helio.npz"
#
#     atran_npz = numpy.load(modeldirpath + telmodelhelio)
#     wav = atran_npz["wav"]
#     flux = atran_npz["flux"]
#     id = atran_npz["id"]
#
#     return wav, flux, id


if __name__ == "__main__":
    # id1, wav1, flux1 = modelopen("atran.smo.1412.dat")
    # id2, wav2, flux2 = modelopen("atran.smo.10958_R28000_1700ft.dat")
    # id3, wav3, flux3 = modelopen("atran.smo.11513_R28000_0ft.dat")
    #
    # plt.figure()
    # plt.plot(wav1,flux1)
    # plt.plot(wav2,flux2-1.0)
    # plt.plot(wav2,flux3-2.0)
    # plt.grid()
    #
    # plt.savefig("modelcomp.png")

    id, wav, flux = modelopen(sys.argv[1])
    dwav = (wav[-1]-wav[0]) / wav.size
    wav_new = numpy.array([wav[0] + i * dwav for i in range(wav.size)])

    numpy.savez(sys.argv[2], id=id, wav=wav_new, flux=flux)

