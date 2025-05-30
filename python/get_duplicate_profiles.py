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
from astropy.visualization import simple_norm
from astropy.table import Table

from PIL import Image

from build_web_cutouts2 import display_image


# define colors - need this for plotting line and fill_between in the same color
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


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


def plot_mstar_sfr_profiles(subdirlist, ncol, nrow, isubplot=[5,5,10,10,15]):
    """
    retrofitting this from the build_web_cutouts

    goal is to plot phot profiles of the stellar mass and sfr images

    PARAMS
    @param subdirlist list containing name of the subdirectory containing the mstar and sfr phot files; will contain more than one dir if there are duplicate observations
    @param ncol number of columns in larger figure
    @param ncow number of rows in larger figure
    @param isubplot list containing subplots for the 5 figures

    RETURN
    nothing, just adds to the figure
    """
    
    for subdirname in subdirlist:
        mstar1 = f"{subdirname}/{subdirname}-mstar-vr_phot.fits"
        mstar2 = f"{subdirname}/{subdirname}-mstar-vcosmic_phot.fits"

        sfr1 = f"{subdirname}/{subdirname}-sfr-vr_phot.fits"
        sfr2 = f"{subdirname}/{subdirname}-sfr-vcosmic_phot.fits"

        ssfr = f"{subdirname}/{subdirname}-ssfr_phot.fits"

        tables = [mstar1,mstar2,sfr1,sfr2,ssfr]
        #isubplot = [1,1,2,2,3]
        icolor = [0,1,0,1,0]
        ytext = [0.5,0.4,0.5,0.4]
        labels = ['Mstar-vr','Mstar-Vcosmic','SFR-vr','SFR-Vcosmic','sSFR']


        #fig = plt.figure(figsize=(12,3))
        #plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95,wspace=.4)


        alphas = [1,.4,.6,.4,.6]
        for i,t in enumerate(tables[:-1]):
            ptab = Table.read(t)
            x = ptab['sma_arcsec']
            #y0 = ptab['flux']
            #yerr = ptab['flux_err']

            y0 = ptab['sb']
            yerr = ptab['sb_err']
            
            # cut the profiles at SNR > 3
            if i == 0:
                snrflag = np.abs(y0/yerr) > 3
                x = x[snrflag]
                y0 = y0[snrflag]        
                yerr = yerr[snrflag]

            y1 = y0+yerr
            y2 = y0-yerr

            plt.subplot(nrow,ncol,isubplot[i])

            #if i < 2:
            plt.fill_between(x,y1,y2,alpha=alphas[i])#,color=mycolors[icolor[i]])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(x,y0)
            if i == 0:
                xmin,xmax = plt.xlim()
            else:
                plt.xlim(xmin,xmax)
            plt.xlabel('SMA (arcsec)',fontsize=16)
            #total = np.max(y0[x < xmax])
            shortlab = labels[i].split('-')[0]
            #label=f"{labels[i]} ({np.log10(total):.2f})"
            #plt.plot(x,y0,'-',label=label,lw=2,color=mycolors[icolor[i]])        
            plt.ylabel(shortlab,fontsize=16)
            if i < 4:
                plt.gca().set_yscale('log')
            plt.gca().set_xscale('log')
            plt.legend(loc='lower right')
            #plt.gca().set_yscale('log')
   
        # add subplot from mstar and sfr profiles
        plt.subplot(nrow,ncol,isubplot[-1])
        mstartab = Table.read(tables[0])
        sfrtab = Table.read(tables[2])
        y0 = sfrtab['flux']/mstartab['flux']
        y0 = (sfrtab['sb']/1.e3)/(mstartab['sb']*1.e7)        
        i += 1
        plt.plot(x,np.log10(y0),'-',lw=2)#,color=mycolors[0])
        #plt.xlim(xmin,xmax)
        plt.gca().set_xscale('log')        
        plt.ylabel('log10(sSFR)',fontsize=16)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        #plt.gca().set_yscale('log')
   

    
def make_plots_mags_cutouts(subdirs,vf, singleflag=False):
    # photutils flux




    # plot enclosed flux        

    plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
    #labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
    alphas = [.4,.4,.4,.4]

    ncol = 5
    nrow = 3
    figs = (20,12)
    if len(subdirs) == 3:
        nrow = 4
        figs = (20,16)
    fig = plt.figure(figsize=figs)
    plt.subplots_adjust(bottom=.05, top=.95,left=.1, right=.95)
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
                plt.subplot(nrow,ncol,7)
            else:
                plt.subplot(nrow,ncol,2)
            plt.fill_between(rphot['sma_arcsec'],y1,y2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(rphot['sma_arcsec'],y0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])


            if j < 2:
                plt.subplot(nrow,ncol,8)
            else:
                plt.subplot(nrow,ncol,3)
            plt.fill_between(t['sma_arcsec'],sb1,sb2,alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],sb0,lss[j],lw=2,label=sd+labels[j],color=mycolors[i])

            
    for i in [2,7,3,8]:
        plt.subplot(nrow,ncol,i)
        if i > 5:
            plt.xlabel('SMA (arcsec)',fontsize=16)
        if i == 2:
            plt.ylabel('magnitude (AB)',fontsize=16)
        elif i == 5:
            plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
            
        #plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()  
        plt.legend(loc='lower right',fontsize=8)
        if (i == 2) | (i == 3):
            plt.title("Rband",fontsize=16)
        elif (i == 7) | (i == 8):
            plt.title("Halpha",fontsize=16)

    # plot jpg in subplot 1
    #jpgfile = glob.glob(fileroot+"/legacy/*.jpg")
    legdir = subdirs[0] + "/legacy/"
    try:
        legacy_jpg = glob.glob(legdir+"*.jpg")[0]
        jpeg_data = Image.open(legacy_jpg)    
        try:
            legacy_g = glob.glob(legdir+"*-g.fits")[0]
            header = fits.getheader(legacy_g)
            imwcs = wcs.WCS(header)
            plt.subplot(nrow,ncol,1,projection=imwcs)
        
        except IndexError:
            print(f"WARNING: no legacy g-band image for {subdirs[0]}")
            plt.subplot(nrow,ncol,1)
        plt.imshow(jpeg_data, origin='lower')
        plt.xlabel('RA (deg)',fontsize=16)
        plt.ylabel('Dec (deg)',fontsize=16)
    except IndexError:
        print(f"WARNING: no legacy jpg image for {subdirs[0]}")
    
    # plot cs subtracted images
    nsubplots = [4,7,8,9,11,12]
    # updating subplots after adding mstar images and profiles
    nsubplots = [6, 11, 12, 13, 17, 18] 
    np = 0
    tels = ['BOK', 'INT', 'HDI', "MOS"]
    for i,sd in enumerate(subdirs):
        for t in tels:
            if t in sd:
                thistel = t
                break
        fileroot = f"{sd}/{sd}"
        cs_gr_phot = fileroot+"-CS-gr.fits"
        csgrdata, csgrheader = fits.getdata(cs_gr_phot, header=True)
        csgrwcs = wcs.WCS(csgrheader)
        
        cs_phot = fileroot+"-CS.fits"
        csdata, csheader = fits.getdata(cs_phot, header=True)
        cswcs = wcs.WCS(csheader)
        maskfile = fileroot+"-R-mask.fits"
        mask = fits.getdata(maskfile)
        mask = mask > 0
        #norm = simple_norm(clipped_data, stretch=stretch,max_percent=percentile2,min_percent=percentile1)

        plt.subplot(nrow,ncol,nsubplots[np],projection=cswcs)
        display_image(csdata,stretch='asinh',percentile1=.5,percentile2=99.5,mask=mask)
        plt.xlabel(sd + "-CS",fontsize=8)
        plt.text(0.95, 0.92, thistel+" CS Halpha", transform=plt.gca().transAxes, color='white',fontsize=14, horizontalalignment='right')
        np += 1
        
        plt.subplot(nrow,ncol,nsubplots[np],projection=csgrwcs)
        display_image(csgrdata,stretch='asinh',percentile1=.5,percentile2=99.5,mask=mask)
        plt.xlabel(sd + "-CSgr",fontsize=8)
        plt.text(0.95, 0.92, thistel+" CS-gr Halpha", transform=plt.gca().transAxes, color='white',fontsize=14, horizontalalignment='right')        
        np += 1

    # add mstar and sfr images
    # plot cs subtracted images
    nsubplots = [4,9,14]
    np = 0

    for i,sd in enumerate(subdirs): # just need to plot one image - hopefully zero is a good choice
        fileroot = f"{sd}/{sd}"
        
        mstarfilename = fileroot+"-mstar-vr.fits"
        mstardata, mstarheader = fits.getdata(mstarfilename, header=True)
        mstarwcs = wcs.WCS(mstarheader)
        
        sfrfilename = fileroot+"-sfr-vr.fits"
        sfrdata, sfrheader = fits.getdata(sfrfilename, header=True)
        sfrwcs = wcs.WCS(sfrheader)

        ssfrfilename = fileroot+"-ssfr.fits"
        ssfrdata, ssfrheader = fits.getdata(ssfrfilename, header=True)
        ssfrwcs = wcs.WCS(ssfrheader)
        
        maskfile = fileroot+"-R-mask.fits"
        mask = fits.getdata(maskfile)
        mask = mask > 0
        #norm = simple_norm(clipped_data, stretch=stretch,max_percent=percentile2,min_percent=percentile1)

        plt.subplot(nrow,ncol,nsubplots[np],projection=mstarwcs)
        display_image(mstardata,stretch='asinh',percentile1=.5,percentile2=99.5,mask=mask,zoom=2)
        plt.xlabel("Mstar * 1e7",fontsize=8)
        #plt.text(0.95, 0.92, "Mstar", transform=plt.gca().transAxes, color='white',fontsize=14, horizontalalignment='right')
        np += 1
        
        plt.subplot(nrow,ncol,nsubplots[np],projection=sfrwcs)
        display_image(sfrdata,stretch='asinh',percentile1=.5,percentile2=99.5,mask=mask,zoom=2)
        plt.xlabel("SFR",fontsize=8)
        #plt.text(0.95, 0.92, thistel+" CS-gr Halpha", transform=plt.gca().transAxes, color='white',fontsize=14, horizontalalignment='right') 
        np += 1

        plt.subplot(nrow,ncol,nsubplots[np],projection=ssfrwcs)
        display_image(ssfrdata,stretch='asinh',percentile1=.5,percentile2=99.5,mask=mask,zoom=2)
        plt.xlabel("sSFR/1e10",fontsize=8)
        #plt.text(0.95, 0.92, thistel+" CS-gr Halpha", transform=plt.gca().transAxes, color='white',fontsize=14, horizontalalignment='right') 
        np += 1
        break
        
    # plot mstar, sfr profiles
    plot_mstar_sfr_profiles(subdirs, ncol, nrow, isubplot=[5,5,10,10,15])        
    if singleflag:
        outfile = f"duplicates/{vf}_profiles_mag_cutouts.png"
    else:
        outfile = f"duplicates/{vf}_profiles_duplicate_mag_cutouts.png"
    plt.savefig(outfile)
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

    singlelist = ([item for item, count in collections.Counter(vfid).items() if (count <  2)])

    # print number of duplicates
    print(f"got {len(duplist)} galaxies with duplicate observations")

    vfid = np.array(vfid)
    dirlist = np.array(dirlist)
    #for vf in duplist:
    #    # get list of subdirectories starting with vf
    #    flag = vfid == vf
    #    subdirs = dirlist[flag]
    #    print("testing subdirs: ",vf, subdirs)
    #    #make_plots(subdirs, vf)
    #    #make_plots_mags(subdirs, vf)
    #    make_plots_mags_cutouts(subdirs, vf, singleflag = True)
    #    break
    
    for vf in singlelist:
        # get list of subdirectories starting with vf
        flag = vfid == vf

        subdirs = dirlist[flag]
        print("testing subdirs: ",vf, subdirs)
        #make_plots(subdirs, vf)
        #make_plots_mags(subdirs, vf)
        make_plots_mags_cutouts(subdirs, vf, singleflag = True)
        break
        
