#!/usr/bin/env

import os
from astropy.table import Table
import glob
# read in virgo-coadds-fullpath.txt
infile = open('virgo-coadds-fullpath.txt')
filenames = []
for line in infile:
    filenames.append(line.strip)

print(filenames[0])

for f in filenames:
    prefix = os.path.basename(rootname).replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    matchfile = glob.glob(prefix+"*")
    if len(matchfile) < 1:
        print("no match for ",f)
