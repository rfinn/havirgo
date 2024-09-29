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
    dirname = d.replace("\n","")
    if os.path.exists(f"{dirname}/halpha-csgr-rfinn-2024-Sep-29.fits"):
        continue
    else:
        print(f"Missing halpha-csgr file for {d}")

