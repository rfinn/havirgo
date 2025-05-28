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
import collections


def make_plots(subdirs,vf):
    # photutils flux
    fig = plt.figure(figsize=(12,12))


    # define colors - need this for plotting line and fill_between in the same color
    mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # plot enclosed flux        

    plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
    #labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
    alphas = [1,.4,.6,.4]

    
    for i,sd in enumerate(subdirs):
        fileroot = f"{sd}/{sd}"
        cs_gr_phot = fileroot+"-CS-gr_phot.fits"
        csgrphot = fits.getdata(cs_gr_phot)

        cs_phot = fileroot+"-CS_phot.fits"
        csphot = fits.getdata(cs_phot)
        
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_phot = self.rimage.replace('.fits','_phot.fits')
        rphot = fits.getdata(r_phot)
        
        tabs = [csgrphot, csphot, rphot]
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
            plt.fill_between(rphot['sma_arcsec'],y1,y2,label=labels[j],alpha=alphas[i],color=sd)
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(rphot['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])


            if j < 2:
                plt.subplot(2,2,4)
            else:
                plt.subplot(2,2,3)
            plt.fill_between(t['sma_arcsec'],sb1,sb2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],sb0,'-',lw=2,color=mycolors[i])

            
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        if i < 2:
            plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        else:
            plt.ylabel('SB (erg/s/cm^2/Hz/arcsec^2)',fontsize=16)
            
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend(loc='lower right')
        if i%2 == 1:
            plt.title("Rband")
        else:
            plt.title("Halpha")
    plt.savefig(f"duplicates/{vf}_duplicate_profiles.png")
    plt.close(fig)

    
if __name__ == '__main__':
    # get list of directories

    if not os.path.exists('duplicates'):
        os.mkdir('duplicates')
    infile = ('conscale_factors.txt','r')

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

    for vf in duplist:
        # get list of subdirectories starting with vf
        flag = vfid == vf

        subdirs = dirlist[flag]
        make_plots(subdirs, vf)
        
