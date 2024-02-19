#!/usr/bin/env python

"""
GOAL:
get elliptical aperture photometry for:
* stellar mass image
* sfr image
* ssfr image

PROCEDURE:
* use the R-band image as a reference
* then run photutils using same aperture as in R-band image
* photwrapper should be able to handle this
* all we need are the profiles for now - can worry about size measurements later


USAGE:


NOTES:
self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = None,image2_filter=self.hafilter, filter_ratio=self.filter_ratio,psf=self.psf_image_name,psf_ha=self.psf_haimage_name)


REQUIRED MODULES:
halphagui.photwrapper

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

from matplotlib import pyplot as plt

homedir = os.getenv("HOME")
sys.path.append(homedir+"/github/halphagui/")

from halphamain import create_output_table

from photwrapper import ellipse
#from fit_profile import profile, dualprofile, rprofile, haprofile, ratio_error


# define colors - need this for plotting line and fill_between in the same color
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def plot_profiles(subdirname,rmax=None):
    
    mstar1 = subdirname+'-logmstar-vr_phot.fits'
    mstar2 = subdirname+'-logmstar-vcosmic_phot.fits'

    sfr1 = subdirname+'-sfr-vr_phot.fits'
    sfr2 = subdirname+'-sfr-vcosmic_phot.fits'

    ssfr = subdirname+'-ssfr_phot.fits'

    tables = [mstar1,mstar2,sfr1,sfr2,ssfr]
    isubplot = [1,1,2,2,3]
    icolor = [0,1,0,1,0]
    ytext = [0.5,0.4,0.5,0.4]
    labels = ['Mstar-vr','Mstar-Vcosmic','SFR-vr','SFR-Vcosmic','sSFR']


    fig = plt.figure(figsize=(12,3))
    plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95,wspace=.4)


    alphas = [1,.4,.6,.4,.6]
    for i,t in enumerate(tables[:-1]):
        ptab = Table.read(t)
        x = ptab['sma_arcsec']
        y0 = ptab['flux']
        yerr = ptab['flux_err']

        if rmax is not None:
            rflag = x < rmax
            x = x[rflag]
            y0 = y0[rflag]        
            yerr = yerr[rflag]
        else:
            rflag = np.ones(len(x),'bool')
        
        y1 = y0+yerr
        y2 = y0-yerr
        
        plt.subplot(1,3,isubplot[i])

        if i < 2:
            plt.fill_between(x,y1,y2,alpha=alphas[i],color=mycolors[icolor[i]])
        # also plot line because you can't see the result when the error is small
        # this should fix issue #18 in Virgo github

        if i == 0:
            xmin,xmax = plt.xlim()
        else:
            plt.xlim(xmin,xmax)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        #total = np.max(y0[x < rmax])
        total = np.max(y0)
        
        shortlab = labels[i].split('-')[0]
        label=f"{labels[i]} ({np.log10(total):.2f})"
        plt.plot(x,y0,'-',label=label,lw=2,color=mycolors[icolor[i]])        
        plt.ylabel(shortlab,fontsize=16)

        #plt.gca().set_yscale('log')
        #plt.gca().set_xscale('log')
        plt.legend(loc='lower right')
        plt.gca().set_yscale('log')
   
    # add subplot from mstar and sfr profiles
    plt.subplot(1,3,3)
    mstartab = Table.read(tables[0])
    sfrtab = Table.read(tables[2])
    y0 = sfrtab['flux']/mstartab['flux']
    i += 1
    plt.plot(x,np.log10(y0[rflag]),'-',lw=2,color=mycolors[0])
    plt.xlim(xmin,xmax)
    plt.ylabel('log10(sSFR)',fontsize=16)
    plt.xlabel('SMA (arcsec)',fontsize=16)
    #plt.gca().set_yscale('log')
   
    outname = subdirname+'-mstar-sfr-profiles.png'
    plt.savefig(outname)

if __name__ == "__main__":
    #print("got here!")
    
    subdirname = sys.argv[1]
    topdir = os.getcwd()

    ###################################################################
    # get galaxy properties from VFID
    ###################################################################    
    vfid = subdirname.split('-')[0]

    maintab = homedir+"/research/Virgo/tables-north/v2/vf_v2_main.fits"
    ephottab = homedir+"/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits"    
    # read in vf_main
    mtab = Table.read(maintab)
    # get redshift
    galindex = np.arange(len(mtab))[mtab['VFID'] == vfid]
        
    # need RA, DEC, radius, BA, PA, like from halphagui
    ra = mtab['RA'][galindex][0]
    dec = mtab['DEC'][galindex][0]
    rad = mtab['radius'][galindex][0]    

    
    ###################################################################    
    # move to subdirectory
    ###################################################################    
    os.chdir(subdirname)


                      
    plot_profiles(subdirname, rmax=1.5*rad)
                      
    os.chdir(topdir)
