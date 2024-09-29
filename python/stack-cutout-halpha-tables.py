#!/usr/bin/env python

'''
GOAL:
* this is an updated Virgo/programs/stack-gui-halpha-tables.py
* I am updating to run on the fits files that are in the cutouts/VFID* directories
* the difference is that these were run AFTER the color-based continuum subtraction with the hand-tuned scale factor,
so the sequencing is different from just running after running the gui in auto mode.


* combine tables that are created from get_gr_cont_phot.py
* astropy.table vstack is really slow, so I am going to try to do it more manually

DATA:
* the most recent directory is /data-pool/Halpha/halphagui-output-20240522/cutouts/
* run from this directory

PROCEDURE:
* count total lines in input files
* create a new empty array with the right number of rows and dataype of the individual tables
* write individual tables into the combined table


OUTPUT:
* combined table
'''

from datetime import date
import glob
import numpy as np
import os
from astropy.io import fits
from astropy.table import Table

# get input files
alldirs = open('virgo-cutouts.txt','r')
alldirlist = alldirs.readlines()

for d in alldirlist:
    print(f"{d}/halpha-csgr-rfinn-2024-Sep-29.fits")
    if os.path.exists(f"{d}/halpha-csgr-rfinn-2024-Sep-29.fits"):
        continue
    else:
        print(f"Missing halpha-csgr file for {d}")

# now gather list of phot output
flist = glob.glob("VFID*/halpha-csgr-rfinn*.fits")
flist.sort()
print(f"Found {len(flist)} files to stack")
print()
# get total number of lines
print('getting total number of lines')
nlines=0
for f in flist:
    t = fits.getdata(f)
    nlines += len(t)

# create output table
print(f"nlines = {nlines}")
outtab = np.zeros(nlines,dtype=t.dtype)

# loop through input tables and write the rows into the output table
print('filling in output table')
startindex=0
for f in flist:
    t = fits.getdata(f)
    n = len(t)
    outtab[startindex:startindex+n] = t
    startindex += n

# write the output table
print('writing output table')
today = date.today()
str_date_today = today.strftime('%Y-%b-%d')
outtab = Table(outtab)
outtab_name = 'csgr-output-combined-{}.fits'.format(str_date_today)

outtab.write(outtab_name,overwrite=True)
