#!/usr/bin/env python

"""
GOAL: overplot radial profiles from photoutils for galaxies with duplicate observations

get galaxies with duplicate observations

for each galaxy with duplicates, create a plot of profiles overlaid, maybe r-band and halpha separate?

which files: 

Halpha:
-CS-gr_phot.fits
-CS_phot.fits

R-band:
-R_phot.fits

run from cutouts directory on 
"""
import os
import numpy as np
import glob
import collections
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from PIL import Image

def make_plots(subdirs,vf):
    # photutils flux
    fig = plt.figure(figsize=(12,12))


    # define colors - need this for plotting line and fill_between in the same color
    mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # plot enclosed flux        

    plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
    #labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
    alphas = [.4,.4,.4,.4]

    
    for i,sd in enumerate(subdirs):
        fileroot = f"{sd}/{sd}"
        cs_gr_phot = fileroot+"-CS-gr_phot.fits"
        csgrphot = fits.getdata(cs_gr_phot)

        cs_phot = fileroot+"-CS_phot.fits"
        csphot = fits.getdata(cs_phot)
        
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_phot =fileroot+'-R_phot.fits'
        rphot = fits.getdata(r_phot)
        
        tabs = [csgrphot, csphot, rphot]
        lss = ['-','--','-']
        labels = ['-gr','','']        
        for j,t in enumerate(tabs):
            y0 = t['flux_erg']            
            y1 = y0 + t['flux_erg_err']
            y2 = y0 - t['flux_erg_err']

            sb0 = t['sb_erg_sqarcsec']
            sb1 = t['sb_erg_sqarcsec']+t['sb_erg_sqarcsec_err']
            sb2 = t['sb_erg_sqarcsec']-t['sb_erg_sqarcsec_err']
            if j < 2:
                plt.subplot(2,2,2)
            else:
                plt.subplot(2,2,1)
            plt.fill_between(rphot['sma_arcsec'],y1,y2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(rphot['sma_arcsec'],y0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])


            if j < 2:
                plt.subplot(2,2,4)
            else:
                plt.subplot(2,2,3)
            plt.fill_between(t['sma_arcsec'],sb1,sb2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],sb0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])

            
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        if i == 0:
            plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        elif i == 2:
            plt.ylabel('SB (erg/s/cm^2/Hz/arcsec^2)',fontsize=16)
            
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend(loc='lower right')
        if i%2 == 0:
            plt.title("Rband",fontsize=16)
        else:
            plt.title("Halpha",fontsize=16)
    plt.savefig(f"duplicates/{vf}_duplicate_profiles.png")
    plt.close(fig)


def make_plots_mags(subdirs,vf):
    # photutils flux
    fig = plt.figure(figsize=(12,12))


    # define colors - need this for plotting line and fill_between in the same color
    mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # plot enclosed flux        

    plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
    #labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
    alphas = [.4,.4,.4,.4]

    
    for i,sd in enumerate(subdirs):
        fileroot = f"{sd}/{sd}"
        cs_gr_phot = fileroot+"-CS-gr_phot.fits"
        csgrphot = fits.getdata(cs_gr_phot)

        cs_phot = fileroot+"-CS_phot.fits"
        csphot = fits.getdata(cs_phot)
        
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_phot =fileroot+'-R_phot.fits'
        rphot = fits.getdata(r_phot)
        
        tabs = [csgrphot, csphot, rphot]
        lss = ['-','--','-']
        labels = ['-gr','','']        
        for j,t in enumerate(tabs):
            y0 = t['mag']
            y1 = t['mag']+t['mag_err']
            y2 = t['mag']-t['mag_err']
 
            sb0 = t['sb_mag_sqarcsec']
            sb1 = t['sb_mag_sqarcsec']+t['sb_mag_sqarcsec_err']
            sb2 = t['sb_mag_sqarcsec']-t['sb_mag_sqarcsec_err']

            if j < 2:
                plt.subplot(2,2,2)
            else:
                plt.subplot(2,2,1)
            plt.fill_between(rphot['sma_arcsec'],y1,y2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(rphot['sma_arcsec'],y0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])


            if j < 2:
                plt.subplot(2,2,4)
            else:
                plt.subplot(2,2,3)
            plt.fill_between(t['sma_arcsec'],sb1,sb2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],sb0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])

            
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        if i == 0:
            plt.ylabel('magnitude (AB)',fontsize=16)
        elif i == 2:
            plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
            
        #plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()  
        plt.legend(loc='lower right')
        if i%2 == 0:
            plt.title("Rband",fontsize=16)
        else:
            plt.title("Halpha",fontsize=16)
    plt.savefig(f"duplicates/{vf}_duplicate_profiles_mag.png")
    plt.close(fig)


def make_plots_mags_cutouts(subdirs,vf):
    # photutils flux
    fig = plt.figure(figsize=(16,16))


    # define colors - need this for plotting line and fill_between in the same color
    mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # plot enclosed flux        

    plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
    #labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
    alphas = [.4,.4,.4,.4]

    
    for i,sd in enumerate(subdirs):
        fileroot = f"{sd}/{sd}"
        cs_gr_phot = fileroot+"-CS-gr_phot.fits"
        csgrphot = fits.getdata(cs_gr_phot)

        cs_phot = fileroot+"-CS_phot.fits"
        csphot = fits.getdata(cs_phot)
        
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_phot =fileroot+'-R_phot.fits'
        rphot = fits.getdata(r_phot)
        
        tabs = [csgrphot, csphot, rphot]

        lss = ['-','--','-']
        labels = ['-gr','','']        
        for j,t in enumerate(tabs):
            y0 = t['mag']
            y1 = t['mag']+t['mag_err']
            y2 = t['mag']-t['mag_err']
 
            sb0 = t['sb_mag_sqarcsec']
            sb1 = t['sb_mag_sqarcsec']+t['sb_mag_sqarcsec_err']
            sb2 = t['sb_mag_sqarcsec']-t['sb_mag_sqarcsec_err']

            if j < 2:
                plt.subplot(3,3,3)
            else:
                plt.subplot(3,3,2)
            plt.fill_between(rphot['sma_arcsec'],y1,y2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(rphot['sma_arcsec'],y0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])


            if j < 2:
                plt.subplot(3,3,6)
            else:
                plt.subplot(3,3,5)
            plt.fill_between(t['sma_arcsec'],sb1,sb2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],sb0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])

            
    for i in [2,3,5,6]:
        plt.subplot(3,3,i)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        if i == 2:
            plt.ylabel('magnitude (AB)',fontsize=16)
        elif i == 5:
            plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
            
        #plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()  
        plt.legend(loc='lower right')
        if i%2 == 0:
            plt.title("Rband",fontsize=16)
        else:
            plt.title("Halpha",fontsize=16)

    # plot jpg in subplot 1

    #jpgfile = glob.glob(fileroot+"/legacy/*.jpg")
    legdir = subdirs[0] + "/legacy/"
    legacy_jpg = glob.glob(legdir+"*.jpg")[0]
    legacy_g = glob.glob(legdir+"*-g.fits")[0]
    jpeg_data = Image.open(legacy_jpg)

    header = fits.getheader(legacy_g)
    imwcs = wcs.WCS(header)
    plt.subplot(3,3,1,projection=imwcs)
    plt.imshow(jpeg_data, origin='lower')
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('Dec (deg)',fontsize=16)
    
    
    plt.savefig(f"duplicates/{vf}_duplicate_profiles_mag_cutouts.png")
    plt.close(fig)
    
if __name__ == '__main__':
    # get list of directories

    if not os.path.exists('duplicates'):
        os.mkdir('duplicates')
    infile = open('conscale_factors.txt','r')

    dirlist = []
    vfid = []
    for line in infile:
        t = line.split()
        dirlist.append(t[0])
        r = t[0].split('-')
        vfid.append(r[0])
    infile.close()

    # get duplicates
    duplist = ([item for item, count in collections.Counter(vfid).items() if (count > 1)])

    # print number of duplicates
    print(f"got {len(duplist)} galaxies with duplicate observations")

    vfid = np.array(vfid)
    dirlist = np.array(dirlist)
    for vf in duplist:
        # get list of subdirectories starting with vf
        flag = vfid == vf

        subdirs = dirlist[flag]
        print("testing subdirs: ",vf, subdirs)
        #make_plots(subdirs, vf)
        #make_plots_mags(subdirs, vf)
        make_plots_mags_cutouts(subdirs, vf)
        break
        
