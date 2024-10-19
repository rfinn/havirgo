#!/usr/bin/env python

"""

first argument = stacked file from halpha gui

second argument = file from running photwrapper on color-based continuum subtracted (csgr) images


USAGE:

Download halpha-galaxies-including-duplicate as csv (https://docs.google.com/spreadsheets/d/12eDw8tZsiyt8atk3yj07o5h5HK7RtzxnODJaXkFzPYk/edit?usp=sharing)

  cd ~/research/Virgo/halpha-tables
  mv ~/Downloads/halpha-galaxies-including-duplicates\ -\ hagalaxies-including-duplicates.csv .
  mv halpha-galaxies-including-duplicates\ -\ hagalaxies-including-duplicates.csv hagalaxies-including-duplicates.csv

from laptop, ~/research/Virgo/halpha-tables

  scp draco:/data-pool/Halpha/halphagui-output-20240522/cutouts/csgr-output-combined-2024-Sep-29.fits .

  scp draco:/data-pool/Halpha/halphagui-output-20240522/*-2024-Sep-29.fits .

  python ~/github/havirgo/python/combine_csgr_gui_outputs.py halphagui-output-combined-2024-Sep-29.fits from-cutouts/csgr-output-combined-2024-Sep-29.fits 

This creates

hgui_csgrphot_combined.fits

This table has all of the same columns as the halphagui-output-combined files, but the photwrapper columns 
have been replace by the columnts in csgr-output-combined-2024-Sep-29.fits b/c these were run on the 
continuum subtracted images that used the gr image and the manually adjusted scale factor. 

"""

import sys
from astropy.table import Table, Column
import numpy as np
from datetime import date

########################################
# get filenames
########################################
hgui_filename = sys.argv[1]
csgr_filename = sys.argv[2]

########################################
# read in tables
########################################
hgui = Table.read(hgui_filename)
csgr = Table.read(csgr_filename)

########################################
# read google spreadsheet
########################################
gsheet_filename = 'hagalaxies-including-duplicates.csv'
gtab = Table.read(gsheet_filename)
# get rid of first row, which is totals
gtab = gtab[1:]
print("len of gtab after slicing = ",len(gtab))

########################################
# only keep vfids in csgr file
########################################
matchflag = np.zeros(len(hgui),'bool')
matchindex = np.zeros(len(hgui),'i')
gmatchindex = np.zeros(len(hgui),'i')
csindices = np.arange(len(csgr))
for i,v in enumerate(hgui['VFID']):
    test = (v == csgr['VFID']) & (hgui['TEL'][i]==csgr['TEL']) & (hgui['DATE-OBS'][i]==csgr['DATE-OBS']) 
    if np.sum(test) > 0:
        matchflag[i] = True
        matchindex[i] = csindices[test][0]
        #print(i,matchindex[i])
        
####################################################
# cut halpha gui table to keep only
# those galaxies that made it through analysis
####################################################
hgui = hgui[matchflag]
matchindex = matchindex[matchflag]

# check that tables have the same length
if len(hgui) == len(csgr):
    print("fear not, all is well...")
else:
    print("tables are not the same length!!!")
    sys.exit()


#################################################
# remove the SMORPH_FLAG and SMORPH_HFLAG so they can be replaced with integer arrays
#################################################
hgui.remove_column('SMORPH_FLAG')
hgui.remove_column('SMORPH_HFLAG')


#################################################
# overwrite hgui w/good columns from csgr table
#################################################
# column names to keep from csgr table

# column names to keep from csgr table
csgr_keep_columns = ['ELLIP_XCENTROID',
                     'ELLIP_YCENTROID',
                     'ELLIP_EPS',
                     'ELLIP_THETA',
                     'ELLIP_GINI',
                     'ELLIP_HGINI',
                     'ELLIP_M20',
                     'ELLIP_HM20',
                     'ELLIP_UNMASKED_AREA',
                     'ELLIP_TOTAL_AREA',
                     'ELLIP_SUM',
                     'ELLIP_SUM_MAG',
                     'ELLIP_ASYM',
                     'ELLIP_ASYM_ERR',
                     'ELLIP_HSUM',
                     'ELLIP_HSUM_MAG',
                     'ELLIP_HASYM',
                     'ELLIP_HASYM_ERR',
                     'R_SKYNOISE',
                     'H_SKYNOISE',
                     'R_SKY',
                     'H_SKY',
                     'ELLIP_R30',
                     'ELLIP_R50',
                     'ELLIP_R90',
                     'ELLIP_HR30',
                     'ELLIP_HR50',
                     'ELLIP_HR90',
                     'R24',
                     'R24_ERR',
                     'R25',
                     'R25_ERR',
                     'R26',
                     'R26_ERR',
                     'R_F25',
                     'R_F25_ERR',
                     'R24V',
                     'R24V_ERR',
                     'R25V',
                     'R25V_ERR',
                     'R_F50',
                     'R_F50_ERR',
                     'R_F75',
                     'R_F75_ERR',
                     'M24',
                     'M24_ERR',
                     'M25',
                     'M25_ERR',
                     'M26',
                     'M26_ERR',
                     'F_30R24',
                     'F_30R24_ERR',
                     'F_R24',
                     'F_R24_ERR',
                     'C30',
                     'C30_ERR',
                     'PETRO_R',
                     'PETRO_R_ERR',
                     'PETRO_FLUX',
                     'PETRO_FLUX_ERR',
                     'PETRO_R50',
                     'PETRO_R50_ERR',
                     'PETRO_R90',
                     'PETRO_R90_ERR',
                     'PETRO_CON',
                     'PETRO_CON_ERR',
                     'PETRO_MAG',
                     'PETRO_MAG_ERR',
                     'HR16',
                     'HR16_ERR',
                     'HR17',
                     'HR17_ERR',
                     'HR_F25',
                     'HR_F25_ERR',
                     'HR_F50',
                     'HR_F50_ERR',
                     'HR_F75',
                     'HR_F75_ERR',
                     'HM16',
                     'HM16_ERR',
                     'HM17',
                     'HM17_ERR',
                     'HF_30R24',
                     'HF_30R24_ERR',
                     'HF_R24',
                     'HF_R24_ERR',
                     'HC30',
                     'HC30_ERR',
                     'HR_F95R24',
                     'HR_F95R24_ERR',
                     'HF_TOT',
                     'HF_TOT_ERR',
                     'HPETRO_R',
                     'HPETRO_R_ERR',
                     'HPETRO_FLUX',
                     'HPETRO_FLUX_ERR',
                     'HPETRO_R50',
                     'HPETRO_R50_ERR',
                     'HPETRO_R90',
                     'HPETRO_R90_ERR',
                     'HPETRO_CON',
                     'HPETRO_CON_ERR',
                     'HPETRO_MAG',
                     'HPETRO_MAG_ERR',
                     'LOG_SFR_HA',
                     'LOG_SFR_HA_ERR',
                     'LOG_SFR_HA_FLAG',
                     'SSFR_IN',
                     'SSFR_IN_ERR',
                     'SSFR_OUT',
                     'SSFR_OUT_ERR',
                     'SMORPH_XCENTROID',
                     'SMORPH_YCENTROID',
                     'SMORPH_RPETRO_CIRC',
                     'SMORPH_RPETRO_ELLIP',
                     'SMORPH_RHALF_ELLIP',
                     'SMORPH_R20',
                     'SMORPH_R80',
                     'SMORPH_GINI',
                     'SMORPH_M20',
                     'SMORPH_F_GM20',
                     'SMORPH_S_GM20',
                     'SMORPH_C',
                     'SMORPH_A',
                     'SMORPH_S',
                     'SMORPH_FLAG',
                     'SMORPH_HXCENTROID',
                     'SMORPH_HYCENTROID',
                     'SMORPH_HRPETRO_CIRC',
                     'SMORPH_HRPETRO_ELLIP',
                     'SMORPH_HRHALF_ELLIP',
                     'SMORPH_HR20',
                     'SMORPH_HR80',
                     'SMORPH_HGINI',
                     'SMORPH_HM20',
                     'SMORPH_HF_GM20',
                     'SMORPH_HS_GM20',
                     'SMORPH_HC',
                     'SMORPH_HA',
                     'SMORPH_HS',
                     'SMORPH_HFLAG',
                     'ELLIP_RA',
                     'ELLIP_DEC']

csgr_keep_columns = ['ELLIP_XCENTROID',
                     'ELLIP_YCENTROID',
                     'ELLIP_EPS',
                     'ELLIP_THETA',
                     'ELLIP_GINI',
                     'ELLIP_HGINI',
                     'ELLIP_M20',
                     'ELLIP_HM20',
                     'ELLIP_UNMASKED_AREA',
                     'ELLIP_TOTAL_AREA',
                     'ELLIP_SUM',
                     'ELLIP_SUM_MAG',
                     'ELLIP_ASYM',
                     'ELLIP_ASYM_ERR',
                     'ELLIP_HSUM',
                     'ELLIP_HSUM_MAG',
                     'ELLIP_HASYM',
                     'ELLIP_HASYM_ERR',
                     'R_SKYNOISE',
                     'H_SKYNOISE',
                     'R_SKY',
                     'H_SKY',
                     'ELLIP_R30',
                     'ELLIP_R50',
                     'ELLIP_R90',
                     'ELLIP_HR30',
                     'ELLIP_HR50',
                     'ELLIP_HR90',
                     'R24',
                     'R24_ERR',
                     'R25',
                     'R25_ERR',
                     'R26',
                     'R26_ERR',
                     'R_F25',
                     'R_F25_ERR',
                     'R24V',
                     'R24V_ERR',
                     'R25V',
                     'R25V_ERR',
                     'R_F50',
                     'R_F50_ERR',
                     'R_F75',
                     'R_F75_ERR',
                     'M24',
                     'M24_ERR',
                     'M25',
                     'M25_ERR',
                     'M26',
                     'M26_ERR',
                     'F_30R24',
                     'F_30R24_ERR',
                     'F_R24',
                     'F_R24_ERR',
                     'C30',
                     'C30_ERR',
                     'PETRO_R',
                     'PETRO_R_ERR',
                     'PETRO_FLUX',
                     'PETRO_FLUX_ERR',
                     'PETRO_R50',
                     'PETRO_R50_ERR',
                     'PETRO_R90',
                     'PETRO_R90_ERR',
                     'PETRO_CON',
                     'PETRO_CON_ERR',
                     'PETRO_MAG',
                     'PETRO_MAG_ERR',
                     'HR16',
                     'HR16_ERR',
                     'HR17',
                     'HR17_ERR',
                     'HR_F25',
                     'HR_F25_ERR',
                     'HR_F50',
                     'HR_F50_ERR',
                     'HR_F75',
                     'HR_F75_ERR',
                     'HM16',
                     'HM16_ERR',
                     'HM17',
                     'HM17_ERR',
                     'HF_30R24',
                     'HF_30R24_ERR',
                     'HF_R24',
                     'HF_R24_ERR',
                     'HC30',
                     'HC30_ERR',
                     'HR_F95R24',
                     'HR_F95R24_ERR',
                     'HF_TOT',
                     'HF_TOT_ERR',
                     'HPETRO_R',
                     'HPETRO_R_ERR',
                     'HPETRO_FLUX',
                     'HPETRO_FLUX_ERR',
                     'HPETRO_R50',
                     'HPETRO_R50_ERR',
                     'HPETRO_R90',
                     'HPETRO_R90_ERR',
                     'HPETRO_CON',
                     'HPETRO_CON_ERR',
                     'HPETRO_MAG',
                     'HPETRO_MAG_ERR',
                     'LOG_SFR_HA',
                     'LOG_SFR_HA_ERR',
                     'LOG_SFR_HA_FLAG',
                     'SSFR_IN',
                     'SSFR_IN_ERR',
                     'SSFR_OUT',
                     'SSFR_OUT_ERR',
                     'SMORPH_XCENTROID',
                     'SMORPH_YCENTROID',
                     'SMORPH_RPETRO_CIRC',
                     'SMORPH_RPETRO_ELLIP',
                     'SMORPH_R20',
                     'SMORPH_R50',                     
                     'SMORPH_R80',
                     'SMORPH_FLUX_CIRC',                     
                     'SMORPH_RHALF_CIRC',
                     'SMORPH_RMAX_CIRC',                     
                     'SMORPH_FLUX_ELLIP',                     
                     'SMORPH_RHALF_ELLIP',
                     'SMORPH_RMAX_ELLIP',                     
                     'SMORPH_GINI',
                     'SMORPH_M20',
                     'SMORPH_F_GM20',
                     'SMORPH_S_GM20',
                     'SMORPH_C',
                     'SMORPH_A',
                     'SMORPH_S',
                     'SMORPH_SERSIC_AMP',
                     'SMORPH_SERSIC_RHALF',
                     'SMORPH_SERSIC_N',
                     'SMORPH_SERSIC_XC',
                     'SMORPH_SERSIC_YC',
                     'SMORPH_SERSIC_ELLIP',                     
                     'SMORPH_SERSIC_THETA',
                     'SMORPH_SERSIC_CHISQ',
                     'SMORPH_SERSIC_FLAG',
                     'SMORPH_SKY_MEAN',
                     'SMORPH_SKY_MED',
                     'SMORPH_SKY_STD',                     
                     'SMORPH_FLAG',
                     'SMORPH_HXCENTROID',
                     'SMORPH_HYCENTROID',
                     'SMORPH_HRPETRO_CIRC',
                     'SMORPH_HRPETRO_ELLIP',
                     'SMORPH_HR20',
                     'SMORPH_HR50',                     
                     'SMORPH_HR80',
                     'SMORPH_HFLUX_CIRC',                     
                     'SMORPH_HRHALF_CIRC',
                     'SMORPH_HRMAX_CIRC',                     
                     'SMORPH_HFLUX_ELLIP',                     
                     'SMORPH_HRHALF_ELLIP',
                     'SMORPH_HRMAX_ELLIP',                     
                     'SMORPH_HGINI',
                     'SMORPH_HM20',
                     'SMORPH_HF_GM20',
                     'SMORPH_HS_GM20',
                     'SMORPH_HC',
                     'SMORPH_HA',
                     'SMORPH_HS',
                     #'SMORPH_HSERSIC_AMP',
                     #'SMORPH_HSERSIC_RHALF',
                     #'SMORPH_HSERSIC_N',
                     #'SMORPH_HSERSIC_XC',
                     #'SMORPH_HSERSIC_YC',
                     #'SMORPH_HSERSIC_ELLIP',                     
                     #'SMORPH_HSERSIC_THETA',
                     #'SMORPH_HSERSIC_CHISQ',
                     #'SMORPH_HSERSIC_FLAG',
                     'SMORPH_HSKY_MEAN',
                     'SMORPH_HSKY_MED',
                     'SMORPH_HSKY_STD',
                     'SMORPH_HSNR_PIXEL',
                     'SMORPH_SNR_PIXEL',                     
                     'SMORPH_HFLAG',
                     'ELLIP_RA',
                     'ELLIP_DEC']


# replace
for c in csgr_keep_columns:
    if c in hgui.colnames:
        hgui[c]=csgr[c][matchindex]
    else:
        print('adding column ',c)
        col = Column(csgr[c][matchindex],name=c)
        hgui.add_column(col)

# remove unwanted columns
for c in hgui.colnames:
    if c.startswith('GAL_2') | c.startswith('GAL_SERSASYM') | c.startswith('GAL_H'):
        hgui.remove_column(c)

###########################################################
# add column "POINTING" in the csgr-output file as galid in hgui
# this is the directory name and allows for easy id of each unique observation
###########################################################
c = Column(csgr['POINTING'][matchindex],name='galid',dtype='S60')
hgui.add_column(c)




###########################################################
# add column "remove" from gsheet
# this is the directory name and allows for easy id of each unique observation
###########################################################
c = Column(np.zeros(len(hgui), 'bool'),name='badflag')
hgui.add_column(c)

c = Column(np.zeros(len(hgui),'f'),name='remove')
hgui.add_column(c)

#c = Column(np.ones(len(hgui), 'bool'),name='cleanflag')
#hgui.add_column(c)
# search though gtab

# if remove == 1, then
# match gtab galid with hgui galid
# set hgui remove_flag to True

remove = gtab['remove'] > .1


#remove = (gtab['remove'] > 0.1) | (gtab['remove'] == 1)
badgals = gtab['galid'][remove]
remove_val = gtab['remove'][remove]
print("first galaxy to remove in gtab = ",badgals[0])

for i,g in enumerate(badgals):
    for j,hgal in enumerate(hgui['galid']):
        # I saved galid in get_csgr_phot as an S40, but this is not long enough for some
        if badgals[i].startswith(hgal):
            #print("found a match ",badgals[i],hgal)
            hgui['badflag'][j] = True
            hgui['galid'][j] = g
            hgui['remove'][j] = remove_val[i]            
            #print(hgui['badflag'][j],hgui['galid'][j])
            break
# let's just save the remove column from the spreadsheet



#c = Column((hgui['remove'] < 0.1),name='clean')
#hgui.add_column(c)

# # this will include galaxies with remove = 0.5
# remove = (gtab['remove'] == 0.5) | (gtab['remove'] == 1)
# testgal = gtab['galid'] == 'VFID2715-KUG1138+327-BOK-20220426-VFID2762' 
# print("test, this should be 0.5",gtab['remove'][testgal],gtab['remove'][testgal] == 0.5, remove[testgal]) 
# badgals = gtab['galid'][remove]
# #print(badgals)
# for g in badgals: # need to skip first line, which is column total
#     try:
#         n = len(g)
#     except TypeError:
#         continue
    
#     match = g == hgui['galid']
#     print(hgui['galid'][match])
#     if np.sum(matchflag) < 1:
#         print("warning, no match for ",g)
#     #if g == 'VFID2715-KUG1138+327-BOK-20220426-VFID2762':
#     #    print("hey", hgui['galid'][match])
#     #print(match)
#     hgui['cleanflag'][match] = False
    
# print("this should be false: ",hgui['cleanflag'][hgui['VFID']=='VFID2715'])
###########################################################
# sort file by galid
###########################################################

# sort the output file?  or just write it out...
sorted_indices = np.argsort(hgui['VFID'])
hgui = hgui[sorted_indices]


###########################################################
# add a column with index in the full vf_v2 tables
###########################################################

fullvfid = np.zeros(len(hgui),'i')
for i in range(len(hgui)):
    fullvfid[i] = int(hgui['VFID'][i].replace('VFID',''))

c = Column(fullvfid,name='VFINDEX',description='row in vf_v2 catalogs')
hgui.add_column(c)


###########################################################
# write out resulting merged table
###########################################################

today = date.today()
str_date_today = today.strftime('%Y-%b-%d')

outtab_name = 'hgui_csgrphot_combined_{}.fits'.format(str_date_today)

hgui.write(outtab_name,format='fits',overwrite=True)

