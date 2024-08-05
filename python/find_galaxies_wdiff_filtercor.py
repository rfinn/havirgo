#!/usr/bin/env

"""
go through all galaxies in halphagui-output-combined-2024-Jul-08.fits 

compare filter_cor in halphagui-output vs filter_cor in vf_v2_halpha.fits

if different, print prefix, tel, date obs, correct filter_cor, filter_cor in vf_v2_halpha.fits

"""

from astropy.table import Table
import os
import math
import sys
import numpy as np
homedir = os.getenv("HOME")


################################
# READ IN HALPHAGUI TABLE
################################
hagui_name = os.path.join(homedir,'research/Virgo/halpha-tables/halphagui-output-combined-2024-Jul-08.fits')
hagui_tab = Table.read(hagui_name)


################################
# READ IN MERGED HALPHA TABLE
################################
vfha_name = os.path.join(homedir,'research/Virgo/tables-north/v2/vf_v2_halpha.fits')
vfha_tab = Table.read(vfha_name)


################################
# LOOP OVER HALPHAGUI TABLE
################################
nprob = 0
print("prefix                              tel date     corr_fc fc_used")
print('----------------------------------------------------------------')

sorted_indices = np.argsort(hagui_tab['VFID'])
for i in sorted_indices:
    j =  vfha_tab['VFID'] == hagui_tab['VFID'][i] 
    #print(hagui_tab['VFID'][i],vfha_tab['VFID'][j])

    vfid = hagui_tab['VFID'][i]
    fc_gui = float(hagui_tab['FILT_COR'][i])
    try:
        fc_vf = float(vfha_tab['FILT_COR'][j][0])
    except:
        print("ERROR converting to float: ",vfha_tab['FILT_COR'][j])
        sys.exit()

    if fc_gui > 0:
        #print(f"{fc_gui:.3f}, {fc_vf:.3f}")
        if not math.isclose(fc_gui,fc_vf,abs_tol=.001):

            s = f"{hagui_tab['prefix'][i]:35} {hagui_tab['TEL'][i]} {hagui_tab['DATE-OBS'][i]} {fc_gui:.4f} {fc_vf:.4f} "
            print(s)
            nprob += 1


print()
print("number of problem galaxies = ",nprob)
