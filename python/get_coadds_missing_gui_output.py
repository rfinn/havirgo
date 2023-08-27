#!/usr/bin/env

import os
from astropy.table import Table
import glob

print()
print("Checking for coadds in virgo-coadds-fullpath.txt with no output file")
print()
# read in virgo-coadds-fullpath.txt
infile = open('virgo-coadds-fullpath.txt')
filenames = []
for line in infile:
    filenames.append(line.strip())

#print(filenames[0])

for f in filenames:
    #print(f)
    prefix = os.path.basename(f).replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    matchfile = glob.glob(prefix+"*")
    if len(matchfile) < 1:
        print("no match for ",f)
