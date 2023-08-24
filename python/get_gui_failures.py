#!/usr/bin/env python

#from astropy.io import fits
from astropy.table import Table
import os
import sys

# read in combined table
tabname = sys.argv[1]

tab = Table.read(tabname)

# find objects with ELLIP_RA of zero
guifailed = tab['ELLIP_XCENTROID'] == 0

# get the set of images

dirname = tab['POINTING'][guifailed]

outfile = 'gui_failures_fullpath.txt'
if os.path.exists(outfile):
    os.remove(outfile)
    
# grep full file list and redirect to new output file
for d in dirname:
    os.system(f'grep {d} virgo-coadds-fullpath.txt >> gui_failures_fullpath.txt')
