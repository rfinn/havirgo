#!/usr/bin/env python

"""
GOAL:
* search cutout directories and find galaxies that don't have legacy cutouts
* run runone_cutout_web.py for each galaxy

RATIONALE:
for whatever reason, getting the cutout images fails in many cases.
Therefore, when I make a new version of the cutouts/data, the most time-consuming
part is having to go back and get images for the cutouts which failed.

This script will streamline the process in one program.

USAGE:
* run from html-dev/cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)

"""

import os
import numpy as np
import glob

from astropy.io import fits
homedir = os.getenv("HOME")

                  

# wrap

if __name__ == '__main__':

    # location of cutout directory on draco
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


    for i,g in enumerate(galnames):
        #print(g)
        vfid = g.split('-')[0]

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
            print('WARNING: no legacy image for ',g)

            # run runone_cutout_web.py but without syncing
            #os.system(f"python ~/github/havirgo/runone_cutout_web.py {g} nosync")
            print(f"will run: python ~/github/havirgo/runone_cutout_web.py {g} nosync")

# sync files to faculty web when done
# sync to web
syncfiles=False
if syncfiles:
    print()
    print("syncing to facultyweb - enter Siena password when prompted...")
    print()
    os.chdir("/data-pool/Halpha/html_dev/")
    os.system("rsync -avz cutouts facultyweb.siena.edu:public_html/virgo/. ")



