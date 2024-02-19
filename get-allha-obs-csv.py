#!/usr/bin/env python

"""
* Create a csv file that lists all of the galaxies observed with halpha
* CSV will have the galname and filter correction
* we will use this to review the halpha galaxies to identify those with nearby bright stars or galaxies
* We will use this list to create a clean halpha sample

"""

import numpy as np
from astropy.table import Table
import os
homedir = os.getenv("HOME")
# read in table that contains full sample

infile = os.path.join(homedir,'research','Virgo','halpha-tables','halphagui-output-combined-2023-Aug-27.fits')

htab = Table.read(infile)

outfile = os.path.join(homedir,'research','Virgo','halpha-tables','hagalaxies-including-duplicates.csv')

out1 = open(outfile,'w')

# sort by VFID

sorted_indices = np.argsort(htab['VFID'])
# write out the header
s = f"galid,filter_cor,bright star,bright galaxy,mask prob,interesting ha,comment\n"
out1.write(s)
for i in sorted_indices:
    # get telescope, date, pointing
    t = htab['POINTING'][i].split('-')
    suffix = t[-3]+'-'+t[-2]+'-'+t[-1]

    galobsname = htab['prefix'][i]+'-'+suffix
    

    s = f"{galobsname},{htab['FILT_COR'][i]:.3f},0,0,0,0, \n"
    out1.write(s)

out1.close()
