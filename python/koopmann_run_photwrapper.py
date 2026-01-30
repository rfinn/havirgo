#!/usr/bin/env python
"""
GOAL:
measure photometry of the continuum subtracted images that are created using the g-r color image


PROCEDURE:
* this will run in each cutout directory
* can be run in parallel

USAGE:
* this takes the directory name as input, like

python ~/github/havirgo/python/get_gr_cont_phot.py VFID6352-NGC5806-BOK-20210418-VFID6406


* to run in parallel

parallel --eta  python ~/github/havirgo/python/get_gr_cont_phot.py :::: virgo-cutouts.txt 


"""
import sys
import os
from datetime import date
import numpy as np

from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table, Column
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

homedir = os.getenv("HOME")
sys.path.append(homedir+"/github/halphagui/")
#from halphamain import create_output_table

from photwrapper import ellipse
from fit_profile import profile, dualprofile, rprofile, haprofile, ratio_error

sys.path.append(homedir+"/github/havirgo/python/")

from get_gr_cont_phot import output_table


if __name__ == '__main__':
    # loop through directory structures
    
