#!/opt/caastro/ext/anaconda/bin/python


import astropy.io.fits as pyfits
import numpy as np
import sys
import os
import errno
import getopt
import math

# global parameters :
debug=0
fitsname="file.fits"
fitsname2="fits2.fits"
check_diagonal=True
limit=0.0000000000000001

def usage():
   print("comparefits.py FITS_FILE1 FITS_FILE2 CHECK_DIAGONAL[default %d] limit[default %.12f]" % (check_diagonal,limit))
   print "\n"

if len(sys.argv) > 1:
   fitsname = sys.argv[1]

if len(sys.argv) > 2:
   fitsname2 = sys.argv[2]

if len(sys.argv) > 3:
   check_diagonal = (int(sys.argv[3]) > 0 )
   
if len(sys.argv) > 4:
   limit = float(sys.argv[4])
   


print("####################################################")
print("PARAMTERS :")
print("####################################################")
print("fitsname        = %s" % fitsname)
print("fitsname2       = %s" % fitsname2)
print("check_diagonal  = %d" % check_diagonal)
print("limit           = %.12f" % (limit))
print("####################################################")

fits = pyfits.open(fitsname)
x_size=fits[0].header['NAXIS1']
# channels=100
y_size=fits[0].header['NAXIS2']
print 'Read fits file %s' % fitsname
print 'FITS size = %d x %d' % (x_size,y_size)

fits2 = pyfits.open(fitsname2)
x_size2=fits2[0].header['NAXIS1']
# channels=100
y_size2=fits2[0].header['NAXIS2']
print 'Read fits file 2 %s' % fitsname
print 'FITS size 2 = %d x %d' % (x_size,y_size)

if x_size!=x_size2 or y_size!=y_size2 :
   print "ERROR : cannot execute operation %s on files of different sizes (%d,%d) != (%d,%d)" % (oper,x_size,y_size,x_size2,y_size2)
   exit;

data1=fits[0].data
data2=fits2[0].data

# print 'BEFORE (%d,%d) = %.2f' % (x_size/2,y_size/2,data[y_size/2][x_size/2])

hdu_out = pyfits.PrimaryHDU()
hdu_out.data = np.random.random((x_size,y_size))
data_out=hdu_out.data

b_equal=True
max_diff=-1e20
x_max_diff=-1
y_max_diff=-1
for y in range(y_size) :
   for x in range(x_size) :   
      if x == y :
         if not check_diagonal :
            continue
            
      diff = data1[y][x] - data2[y][x]     
      
      if math.fabs(diff) > limit :
         b_equal = False
      
      if math.fabs(diff) > max_diff :
         max_diff = math.fabs(diff)
         x_max_diff = x
         y_max_diff = y


if b_equal :
   print("COMPARISON RESULT : Images are EQUAL !!!\n")
else :
   print("COMPARISON RESULT : Images differ by maximum value = %.12f at pixel (%d,%d)" % (max_diff,x_max_diff,y_max_diff))   
        