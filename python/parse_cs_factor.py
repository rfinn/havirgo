#!/usr/bin/env python

"""
SETTING UP TO RUN:

* download halpha-galaxies-including-duplicate 
https://docs.google.com/spreadsheets/d/12eDw8tZsiyt8atk3yj07o5h5HK7RtzxnODJaXkFzPYk/edit?usp=sharing

* move and rename file to remove spaces

cd /home/rfinn/research/Virgo/halpha-tables

mv ~/Downloads/halpha-galaxies-including-duplicates\ -\ hagalaxies-including-duplicates.csv hagalaxies-including-duplicates.csv




PROCEDURE:
* read in halpha-galaxies-including-duplicates

* read in output from halpha gui

* for galaxies with no CS2
  - get the wrong filter correction for galaxies with no CS2
    - this should get pulled from vf_v2_halpha.fits
  - get the continuum oversubtraction correction (this depends on the telescope/filter combo)
  - then computs CS2 as 

     CS2 = CONSCALE/FILT_TRANS_CORR/CONT_OVERSUBTRACTION_CORR

* write out a file containing : DIRNAME (galid) CS2

* this can then be fed into subtract continuum

"""

from astropy.table import Table
import os
import numpy as np

######################################################################
###  FILTER DEFINITIONS
######################################################################
# Halpha filter width in angstrom
filter_width_AA = {'BOK':80.48,'HDI':80.48,'INT':95,'MOS':80.48,'INT6657':80}

# central wavelength in angstroms
filter_lambda_c_AA = {'BOK':6620.52,'HDI':6620.52,'INT':6568,'MOS':6620.52,'INT6657':6657}


# integral of filter transmission
# calculated in Halpha-paper1.ipynb
filter_Rlambda = {"KPNO_Ha+4nm": 78.58, "WFC_Ha": 84.21, "WFC_Ha6657": 71.96,\
                  "KPNO_R" : 1341.54, "KPNO_r" : 1283.47, "BASS_r": 1042.18, "WFC_r": 1097.07}

# from oversubtraction of continuum, a la Gavazzi+2006
# take telescope as the key
halpha_continuum_oversubtraction = {'BOK':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["BASS_r"]),\
                            'HDI':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["KPNO_r"]),\
                            'INT':(1 +filter_Rlambda["WFC_Ha"]/filter_Rlambda["WFC_r"]),\
                            'MOS':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["KPNO_R"]),\
                            'INT6657':(1 +filter_Rlambda["WFC_Ha6657"]/filter_Rlambda["WFC_r"])}

######################################################################
###  FUNCTIONS
######################################################################

def get_params_from_name(image_name):
    #print(t)
    tels = ['BOK','HDI','INT','MOS']
    for t in tels:
        if t in image_name:
            telescope = t
            break
    t = os.path.basename(image_name).split('-')
    for item in t:
        if item.startswith('20'):
            dateobs = item
            break
    pointing = t[-1]

    return telescope,dateobs,pointing

def get_conscale2():
    pass

if __name__ == '__main__':

    ##################################################################
    # read in continuum scale factors
    ##################################################################        
    homedir = os.getenv("HOME")
    infile = homedir+'/research/Virgo/halpha-tables/hagalaxies-including-duplicates.csv'
    hagals = Table.read(infile)

    # totals for each column are in the first row, so remove that
    hagals = hagals[1:]

    ##################################################################
    # read in vf_v2_halpha.fits to get the filter correction
    # this is wrong for galaxies with duplicate observations
    # but I need this to recontruct CS2 from CONSCALE
    ##################################################################    
    halphagui_table = os.path.join(os.getenv("HOME"),'research/Virgo/tables-north/v2/vf_v2_halpha.fits')
    
    vhalpha = Table.read(halphagui_table)

    nocs2_flag = hagals['CS2'].mask
    nocs2_indices = np.arange(len(hagals))[nocs2_flag]

    alltelescope = []
    allobservation = []
    allcontosub = np.zeros(len(hagals))
    newcs2 = hagals['CS2'].copy()
    print("i   CS   OSUB  FILT  CS2")
    print("---------------------------")    
    for i in range(len(hagals)):
        dirname = hagals['galid'][i]
        
        telescope,dateobs,pointing = get_params_from_name(dirname)
        vfid = dirname.split('-')[0]

        if vfid == 'VFID3992': # outside footprint of INT, should not be in sample
            badgal_index = i
            continue
        alltelescope.append(telescope)
        allobservation.append(telescope+'-'+dateobs+'-'+pointing)
        allcontosub[i] = halpha_continuum_oversubtraction[telescope]
        if nocs2_flag[i]:

            
            ##################################################################
            # get continuum oversubtraction for each galaxy
            ##################################################################        



            vfindex = vhalpha['VFID'] == vfid

            if hagals['CONSCALE'].mask[i]:
                cs = float(hagals['CSAUTO'][i])
            else:
                cs = float(hagals['CONSCALE'][i])

            t = cs/halpha_continuum_oversubtraction[telescope]/vhalpha['FILT_COR'][vfindex][0]

            if t == np.nan:
                print(f"{dirname}: {i:03d} {cs:.2f} {halpha_continuum_oversubtraction[telescope]:.3f} {vhalpha['FILT_COR'][vfindex][0]:.3f} {t:.3f}")
            newcs2[i] = t

    alltelescope = np.array(alltelescope)
    allobservation = np.array(allobservation)

##################################################################
# remove row with wacky galaxy outside INT FOV
##################################################################        
goodflag = np.arange(len(hagals)) != badgal_index

hagals = hagals[goodflag]
newcs2 = newcs2[goodflag]

##################################################################
# write out file to use with parallel and subtract_continuum.py
##################################################################        

outfile = open(homedir+'/research/Virgo/halpha-tables/conscale_factors.txt','w')

for i in range(len(hagals)):
    s = f"{hagals['galid'][i]} {newcs2[i]} \n"
    outfile.write(s)

outfile.close()
