import sys,os,datetime
from pyraf import iraf

iraf.noao()
iraf.onedspec()

#Description:
#   This "PySpecshift.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to shift the spectrum.
#
#Usage:
#
#   $ python PySpecshift.py <input> <output> <shift>
#
#   input -- the fits file to be shifted
#   output -- the shifted fits file
#   shift -- the value of shift
#

def PySpecshift(input, output, shift):
    iraf.scopy(input, output)
    iraf.specshift(output, shift)

if __name__ == "__main__":

    filename = sys.argv[1:]
    PySpecshift(filename[0], filename[1], float(filename[2]))