#!/usr/bin/env python

'''
GOAL:
* create index web page that points to html pages for all cutouts 

USAGE:
* run from html-dev/cutouts directory


'''

import os
import numpy as np
import glob

from astropy.io import fits
homedir = os.getenv("HOME")

###########################################################
####  FUNCTIONS
###########################################################

def check_gal_list(galnames,outdir):
    vfindices = np.arange(len(vfmain))
    galindex = 1    
    ids = []
    for i,g in enumerate(galnames):
        #print(g)
        vfid = g.split('-')[0]
        vfindex = vfindices[vfmain['VFID'] == vfid][0]
        ra = vfmain['RA'][vfindex]
        dec = vfmain['DEC'][vfindex]
        ids.append(vfmain['VFID'][vfindex])

        jpg_path = os.path.join(outdir,g)
        search_path = os.path.join(jpg_path,'*legacy*.jpg')
        #print(search_path)
        #legacy_jpg = glob.glob(search_path)[0]            
        try:
            #print()
            #print("looking for legacy image")
            #print(glob.glob(search_path))
            legacy_jpg = glob.glob(search_path)[0]
            legacy_flag = True                
        except:
            legacy_flag = False
            legacy_jpg = None
            print()
            print()            
            print('WARNING: no legacy image for ',g)
            print("trying to download again")
            print()
            print()
            #print('\t Skipping galaxy for now')
            #continue

            # try to rebuild website
            # download cutouts
            cutout_dir = g
            current_dir = os.getcwd()

            # move to cutouts data dir
            os.chdir("/data-pool/Halpha/halphagui-output-20230818/cutouts/")
            s = f"python ~/github/HalphaImaging/python3/generate_all_cutout_plots.py --onegal {cutout_dir}"
            os.system(s)

            # build webpage
            s = f"python ~/github/havirgo/python/build_web_cutouts2.py --oneimage {cutout_dir}"
            os.system(s)

            os.chdir(current_dir)
            try:
                legacy_jpg = glob.glob(search_path)[0]
                legacy_flag = True                
            except:
                legacy_flag = False
                legacy_jpg = None
                print('ERROR: no legacy image for ',g)


if __name__ == '__main__':
    # work through coadd directory
    #global vmain

    VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'    
    vfmain = fits.getdata(VFMAIN_PATH)

    # updating for v2 catalogs
    VFFIL_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_environment.fits'    
    vffil = fits.getdata(VFFIL_PATH)

    VFMAGPHYS_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_magphys_10-Jul-2023.fits'    
    vfmagphys = fits.getdata(VFMAGPHYS_PATH)
    
    outdir = '/data-pool/Halpha/html_dev/cutouts/'    

    # this should contain a list of all the galaxy folders
    flist1 = os.listdir(outdir)
    flist1.sort()
    galnames=[]
    for i,subdir in enumerate(flist1): # loop through list
        

        #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
        if (os.path.isdir(subdir)) & (subdir.startswith('VF')):
            #print('adding ',subdir)
            galnames.append(subdir)
    print('number of subdirectories = ',len(galnames))
    check_gal_list(galnames,outdir)
