#!/usr/bin/env python

"""
GOAL:

reproject legacy images to halpha pixel scale

PROCEDURE:
input: get directory name

Run from the cutouts directory that contains a subdirectory for each galaxy.  
Note galaxies that were observed more than once will have a directory for each observation.


To run with parallel from the cutouts directory:

create a list of all the subdirectories:

ls -d VFID* > virgo-cutouts.txt

Then call parallel:

parallel --eta  python ~/github/havirgo/python/reproject_images_2ha.py :::: virgo-cutouts.txt 

"""
import sys
import os
from astropy.io import fits
import glob
from reproject import reproject_interp

def reproject_image(infile, reffile, outname):
    """reproject infile to reffile image"""
    if os.path.exists(outname):
        print("reprojected image exists - not redoing it")
        return
    
    hinfile = fits.open(infile)
    href = fits.open(reffile)
    # reproject input to referece image
    outim,footprint = reproject_interp(hinfile,href[0].header)

    fits.writeto(outname,outim,href[0].header,overwrite=True)
    hinfile.close()
    href.close()

if __name__ == '__main__':
    dirname = sys.argv[1]

    #get CS image
    reffile = os.path.join(dirname,dirname+'-CS.fits')

    if not os.path.exists(reffile):
        print("can't find CS image - exiting ",reffile)
        sys.exit()
    #get legacy/*.fits
    # get each filter separately to avoid creating -ha-ha-ha.fits when running multiple times...
    legacyr = glob.glob(os.path.join(dirname,'legacy/*r.fits'))
    legacyg = glob.glob(os.path.join(dirname,'legacy/*g.fits'))
    legacyz = glob.glob(os.path.join(dirname,'legacy/*z.fits'))        
    legacy_images = legacyr+legacyg+legacyz #glob.glob(os.path.join(dirname,'legacy/*r.fits'))
    legacy_images.sort()

    for infile in legacy_images:
        outname = infile.replace('.fits','-ha.fits')
        
        # reproject onto halpha wcs
        reproject_image(infile,reffile,outname)
