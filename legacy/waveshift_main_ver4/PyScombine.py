import sys,os,datetime
from pyraf import iraf

iraf.noao()
iraf.onedspec()


#Description:
#   This "Waveshift_main.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to combine the fits files.
#
#Usage:
#
#   $ python PyScombine.py <input> <output>
#
#   input -- the list of fits files
#   output -- the outptu fits file
#

def PyScombine(input,output):
    iraf.scombine("@%s" % input, output, combine="average" )

if __name__ == "__main__":

    filename = sys.argv[1:]
    iraf.scombine("@%s" % filename[0], filename[1], combine="average")

