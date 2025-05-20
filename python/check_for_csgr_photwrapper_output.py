from datetime import date
import glob
import numpy as np
import os
from astropy.io import fits
from astropy.table import Table
import sys

# get input files
alldirs = open('virgo-cutouts.txt','r')
alldirlist = alldirs.readlines()

if len(sys.argv) > 1:
    newdate = sys.argv[1]
    print(f"got a new date to use: {newdate}")
else:
    newdate = '2024-Sep-29'
for d in alldirlist:
    
    dirname = d.rstrip()
    #print(f"{d}/halpha-csgr-rfinn-2024-Sep-29.fits")
    #print(dirname)
    if os.path.exists(f"{dirname}/halpha-csgr-rfinn-{newdate}.fits"):
        continue
    else:
        print(f"Missing halpha-csgr file for {dirname}")
        if len(sys.argv) > 2:
            os.system(f"python ~/github/havirgo/python/get_gr_cont_phot.py {dirname}")
