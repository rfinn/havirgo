#!/usr/bin/env

import numpy as np

import os
homedir = os.getenv("HOME")

import sys
sys.path.append(homedir+'/github/HalphaImaging/python3/')
import plot_cutouts_ha as pch
sys.path.append(homedir+'/github/Virgo/programs/')
from readtablesv2 import vtables

import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
from astropy.table import Table

from PIL import Image

import glob

homedir = os.getenv("HOME")

import warnings
warnings.filterwarnings('ignore')


#%matplotlib inline
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plotdir = homedir+'/research/Virgo/plots/halpha/'

zoom_coords_INT = {'VFID5889':[100,1300,850,2100],\
                   'VFID5851':[255,1500,740,2120],\
                   'VFID5855':[820,1200,700,1400],\
                   'VFID5842':[900,1300,700,1450],\
                   'VFID5859':[50,150,30,145],\
                   'VFID5892':[240,780,250,690]
                   }


# new coords to fit HI contours
zoom_coords_INT = {'VFID5889':[100,1300,850,2100],\
                   'VFID5851':[255,1500,740,2120],\
                   'VFID5855':[730,1290,640,1460],\
                   'VFID5842':[800,1400,650,1500],\
                   'VFID5859':[1,175,1,175],\
                   'VFID5892':[240,780,250,690]
                   }


HIdir = homedir+'/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/'
HI_file = {'VFID5889':None,\
           'VFID5851':HIdir+'J1355_fin_lw05_bpcorr_2_mom0.fits',\
           'VFID5855':HIdir+'J1355_fin_lw05_bpcorr_6_mom0.fits',\
           'VFID5842':HIdir+'J1355_fin_lw05_bpcorr_5_mom0.fits',\
           'VFID5859':None,\
           'VFID5892':None
}


zoom_coords_HDI = {'VFID5889':[250,1200,750,1700],\
                   'VFID5842':[400,1300,550,1130],\
                   'VFID5892':[200,600,190,550],\
                   'VFID5859':[35,115,25,115]
                   }


acbfrac = {'VFID5889':.045,\
           'VFID5851':.045,\
          'VFID5855':0.08,\
          'VFID5842':0.08,\
          'VFID5859':0.05,\
          'VFID5892':.075
          }


afigsize = {'VFID5889':[16,3.5],\
            'VFID5851':[16,3.5],\
           'VFID5855':[16,4.5],\
           'VFID5842':[16,4.5],\
           'VFID5859':[16,3.5],\
           'VFID5892':[16,3.]
           }

alevels = {'VFID5889':[4.5],\
           'VFID5851':[4.5],\
           'VFID5855':[4],\
           'VFID5842':[4],\
           'VFID5859':[3.55],\
           'VFID5892':[4.3]
           }


HIfiles = {'VFID5859':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_5_mom0.fits',\
           'VFID5855':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_6_mom0.fits',\
           'VFID5851':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_2_mom0.fits',\
           }

###########################################
## FUNCTIONS
###########################################
def plot_HI_contours(ax,HIfilename,xmin=None,xmax=None,ymin=None,ymax=None,levels=None):
    """ plot HI contours, given current axis + reference image header """
    # borrowing from alma 2023 proposal
    ncontour=5

    hdu = fits.open(HIfilename)[0]
    HI_WCS = WCS(hdu.header)

    if levels is None:
        levels = 3**np.arange(ncontour)+1
    ax.contour(hdu.data,transform=ax.get_transform(HI_WCS),levels=levels,colors='white',alpha=.5,lw=1)    

def update_ngc5348_mask():
    imname = ''
    hdu = fits.open()

    # mask out all columns with with x > 1290
    flag 

    # mask out rows with y > 1680

        


def get_rad_fluxfrac(pfit,frac=0.9):
    from scipy.interpolate import interp1d
    N = 10
    x = np.linspace(0,N*pfit.r_eff.value,100*N)
    flux = pfit(x)*2*np.pi*x*(x[1]-x[0])
    integral2 = np.cumsum(flux)
    
    interp_profile = interp1d(integral2,x)
    return interp_profile(frac*integral2[-1])

def get_legacy_images():

    # getting image center from topcat, after matching with Tempel+2017 Table 1
    # and creating subset for GroupID ==
    ramin,ramax = 208.02354, 209.30849
    decmin,decmax = 4.30731, 5.77169
    # not using the full range of coordinates b/c dwarfs are far flung
    #ra = 0.5*(ramin+ramax)

    # selected RA and DEC from legacy viewer
    ra = 208.744
    dec = 5.1388
    

    size = np.max([decmax-decmin,ramax-ramin])
    # using RA of source E of two biggies
    #size = 2*(209.2316 - ra)
    print(f"imsize = {size} (deg)")
    # convert from deg to arcsec
    imsize = size*3600*1.2
    

    fits_name, jpeg_name = pch.get_legacy_images(ra,dec,galid='NGC5364_group',imsize=imsize,band='r',subfolder='legacy',makeplots=True)
    return  fits_name, jpeg_name



# trying to display stellar mass and sSFR images better
def display_image(image,percent=99.9,lowrange=False,mask=None,sigclip=True,cmap='gray_r'):
    #from astropy.visualization import SqrtStretch
    #from astropy.visualization.mpl_normalize import ImageNormalize

    from astropy.stats import sigma_clip
    from astropy.visualization import simple_norm

    lowrange=False
    if sigclip:
        clipped_data = sigma_clip(image,sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap=cmap,origin='lower')
    #v1,v2=scoreatpercentile(image,[.5,99.5])            
    #plt.imshow(image, cmap='gray_r',vmin=v1,vmax=v2,origin='lower')    

def plot_mstar_sfr(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],cbfrac=.08,cbaspect=20,\
                   clevels=[4],contourFlag=True):
    #%matplotlib inline
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name


    cwd = os.getcwd()
    os.chdir(dirname)
    massim = dirname+"-logmstar-vr.fits"
    sfrim = dirname+"-sfr-vr.fits"
    ssfrim = dirname+"-ssfr.fits"
    mask = dirname+'-R-mask.fits'
    titles = ['log Mstar','SFR','log sSFR']
    vmin = [2,0,-11.5]
    vmax = [6,1.e-6,-9]
    allim = [massim,sfrim,ssfrim]

    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)
    if 'INT' in dirname:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
        except KeyError:
            print("no xmin,xmax in dictionary")
    else:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
        except KeyError:
            print("no xmin,xmax in dictionary")
    


    myfigsize=afigsize[vfid]
    mycbfrac=acbfrac[vfid]
    myclevels=alevels[vfid]

    if 'VFID5892' in dirname:
        cbaspect = 10
        
    plt.figure(figsize=(myfigsize[0],myfigsize[1]))

    plt.subplots_adjust(wspace=0.1)
    maskdat = fits.getdata(mask)

    allax = []
    for i, im in enumerate(allim):
        plt.subplot(1,4,i+2)
        dat = fits.getdata(im)
        mdat = np.ma.array(dat,mask=maskdat)
        if xmin == None:
            mdat = mdat
        else:
            mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],origin='lower',interpolation='nearest')
        else:
            display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
        
        if i in [1,2]:
            plt.xticks([],[])
            plt.yticks([],[])

        elif not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        plt.title(titles[i],fontsize=20)
        plt.colorbar(fraction=mycbfrac,aspect=cbaspect)
        if i == 1:
            t = dirname.split('-')
            plt.xlabel(t[0]+' '+t[1],fontsize=20)
        # plot contours from mass
        allax.append(plt.gca())
    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band
    plt.subplot(1,4,1,projection=legwcs)
    plt.imshow(jpeg_data)
    # set limits in ra,dec
    x,y = legwcs.world_to_pixel(sky)
    # convert ramin,ramax and decmin,decmax to (x,y)
    print(sky)
    plt.axis([x[0],x[1],y[0],y[1]])
    #print(x[0],x[1],y[0],y[1])
    if contourFlag:
        # get contours from logmstar image
        hdu = fits.open(massim)
        contour_data = hdu[0].data
        contour_header = hdu[0].header
        hdu.close()
        mcontour_data = np.ma.array(contour_data,mask=maskdat)
        for ax in allax:
            ax.contour((mcontour_data[ymin:ymax,xmin:xmax]),levels=myclevels, colors='k',linestyles='-',linewidths=1)#,transform=ax.get_transform(WCS(contour_header))
    #plt.contour(contour_data[y[0]:y[1],x[0]:x[1]],levels = [4,5,6],colors='c',alpha=0.5)
    plt.title("Legacy grz",fontsize=20)
   
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-ssfr.png')


    os.chdir(cwd)


def plot_mstar_sfr_profiles(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,\
                            Re_mstar=None,Re_sfr=None,R90_mstar=None,R90_sfr=None):
    """
    same plot as mstar_sfr, but swap out ssfr for radial profiles in the 4th panel

    """
    #%matplotlib inline
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name


    cwd = os.getcwd()
    os.chdir(dirname)
    massim = dirname+"-logmstar-vr.fits"
    sfrim = dirname+"-sfr-vr.fits"

    
    ssfrim = dirname+"-ssfr.fits"
    mask = dirname+'-R-mask.fits'
    titles = ['log Mstar','SFR','log sSFR']
    vmin = [2,-.2e-5,-11.5]
    vmax = [6,.6e-5,-9]
    allim = [massim,sfrim,ssfrim]
    allim = [massim,sfrim]    

    
    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)
    if 'INT' in dirname:
        xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
    else:
        xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
    


    myfigsize=afigsize[vfid]
    mycbfrac=acbfrac[vfid]
    myclevels=alevels[vfid]

    if 'VFID5892' in dirname:
        cbaspect = 10
        
    plt.figure(figsize=(myfigsize[0],myfigsize[1]))

    plt.subplots_adjust(wspace=0.2,bottom=.15)
    maskdat = fits.getdata(mask)

    if contourFlag:
        # get contours from logmstar image
        hdu = fits.open(massim)
        contour_data = hdu[0].data
        contour_header = hdu[0].header
        contour_WCS = WCS(contour_header)
        hdu.close()
        mcontour_data = np.ma.array(contour_data,mask=maskdat)

    
    allax = []
    for i, im in enumerate(allim):
        hdu = fits.open(im)[0]
        dat = fits.getdata(im)
        imwcs = WCS(fits.getheader(im))

        dat = hdu.data
        imwcs = WCS(hdu.header)
        plt.subplot(1,4,i+2,projection=imwcs)


        mdat = np.ma.array(dat,mask=maskdat)
        #if xmin is None:
        #    mdat = mdat
        #else:
        #    mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],origin='lower',interpolation='nearest')

        else:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            plt.imshow(mdat,cmap='viridis',vmin=vmin[i],vmax=vmax[i])
        

        #############################################################
        # add HI contour to halpha image
        #############################################################
        if i == 1:
            # check if HI moment zero map is available
            HIfilename = HI_file[vfid]
            if HIfilename is not None:
                print("HIfilename = ",HIfilename)
                plot_HI_contours(plt.gca(),HIfilename)


        #############################################################
        # add contours from stellar mass image
        #############################################################    
        if contourFlag:
            ax = plt.gca()
            ax.contour(mcontour_data,levels=myclevels, colors='k',linestyles='-',linewidths=1,transform=ax.get_transform(contour_WCS))

                
        plt.title(titles[i],fontsize=18)
        plt.colorbar(fraction=mycbfrac,aspect=cbaspect)
        # plot contours from mass

        if xmin is not None:
            plt.axis([xmin,xmax,ymin,ymax])
        if i in [0,1]:
            plt.xticks([],[])
            plt.yticks([],[])

        elif not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        
        allax.append(plt.gca())


    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band

    
    plt.subplot(1,4,1,projection=legwcs)
    plt.imshow(jpeg_data)

    if xmin is not None:
        # plot the legacy image in panel 1
        xcoords = np.array([xmin,xmax])
        ycoords = np.array([ymin,ymax])
    
        # get ramin,ramax and decmin,decmax from SFR image
        sfrim = dirname+"-sfr-vr.fits"
        header = fits.getheader(sfrim)
        #print(header)
        wcs = WCS(header)
        sky = wcs.pixel_to_world(xcoords,ycoords)
    
        # set limits in ra,dec
        x,y = legwcs.world_to_pixel(sky)
        # convert ramin,ramax and decmin,decmax to (x,y)
        #print(sky)
        plt.axis([x[0],x[1],y[0],y[1]])
    t = dirname.split('-')
    #plt.text(.05,.02,t[0],fontsize=20,transform=plt.gca().transAxes,horizontalalignment='left',color='white')
    
    plt.title(f"{t[0]} Legacy grz",fontsize=18)
    
    #print(x[0],x[1],y[0],y[1])
    #############################################################
    # add HI contour to legacy image
    #############################################################    
    # check if HI moment zero map is available
    HIfilename = HI_file[vfid]
    if HIfilename is not None:
        print("HIfilename = ",HIfilename)
        plot_HI_contours(plt.gca(),HIfilename)

        

    

    #################################################################
    # plot profiles in the 4th panel
    #################################################################    
    plt.subplot(1,4,4)
    #photometry files
    massphot = massim.replace('.fits','_phot.fits')
    sfrphot = sfrim.replace('.fits','_phot.fits')
    mtab = Table.read(massphot)
    stab = Table.read(sfrphot)    
    tables = [mtab,stab]
    labels = ['$M_\star$','$SFR$']    
    
    for i,t in enumerate(tables):
        x = t['sma_arcsec']
        y = t['sb']
        yerr = t['sb_err']
        if ('VFID5859' in dirname) | ('VFID5892' in dirname):
            norm_factor = np.median(y[0:5])
        else:
            norm_factor = np.median(y[0:15])
        plt.plot(x,y/norm_factor,'bo',c=mycolors[i],label=labels[i])
        #plt.fill_between(x,(y+yerr)/norm_factor,y2=(y-yerr)/norm_factor,color=mycolors[i],alpha=.5)
    # add R50 if it's provided
    if Re_mstar is not None:
        plt.axvline(x=Re_mstar,ls='--',lw=2,color=mycolors[0],label='$R_e(M_\star)$')
    if Re_sfr is not None:
        plt.axvline(x=Re_sfr,ls='--',color=mycolors[1],label='$R_e(SFR)$')
    if R90_mstar is not None:
        plt.axvline(x=R90_mstar,ls='-',lw=2,color=mycolors[0],label='$R_{90}(M_\star)$')
    if R90_sfr is not None:
        plt.axvline(x=R90_sfr,ls='-',color=mycolors[1],label='$R_{90}(SFR)$')
    plt.legend(bbox_to_anchor=(1.02,0.95))
    plt.gca().set_yscale('log')
    plt.ylim(.003,3.5)
    #print("rmax = ",rmax)    
    if rmax is not None:
        
        plt.xlim(-1,rmax)

    #plt.ylim(-.1,1.5)
    plt.xlabel("SMA (arcsec)",fontsize=16)
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-profiles-4panel.png')


    os.chdir(cwd)



    

def fit_profiles(rp,hp,weights=None,rmax=None,fixN=False,labels = ['logMstar','logSFR'],log1Flag=False):
    """
    PARAMS: 
    * rp = r-band photometry file
    * hp = halpha photometry file
    
    RETURN:
    * rfit = r-band model
    * hfit = halpha model
    """
    # this is using R band and CS Halpha
    import scipy.optimize
    from astropy.modeling import models, fitting, functional_models
    
    ###################################################
    # fit r-band
    ###################################################
    #flag = rp['sma_arcsec'] < 140
    flag = rp['sb'] > 0 & (rp['sb'] != np.nan)
    flag = np.ones(len(rp),'bool')
    x = rp['sma_arcsec'][flag]
    y = rp['sb'][flag]
    if log1Flag:
        y = 10.**y

    t_init = models.Sersic1D(amplitude=1,r_eff=50,n=2,fixed={'n':fixN})
    fit_t = fitting.LevMarLSQFitter()

    # add weight that scales with the area of the ellipse
    # should be 2 pi r dr, but just using x for now.
    rfit = fit_t(t_init, x, y, weights=x*2*np.pi,maxiter=400)

    ###################################################
    # fit halpha
    ###################################################
    #flag = (hp['sb'] > 0) #& (hp['sma_arcsec'] > 5)
    # force it to fit zeros as well
    flag = np.ones(len(hp),'bool')
    x = hp['sma_arcsec'][flag]
    y = hp['sb'][flag]

    h_init = models.Sersic1D(r_eff=50,n=2,fixed={'n':fixN})

    fit_t = fitting.LevMarLSQFitter()
    # add weight that scales with the area of the ellipse
    hfit = fit_t(h_init, x, y, weights=x,maxiter=200)

    ###################################################
    # plot results
    ###################################################
    fits = [hfit,rfit]
    phots = [hp,rp]
    names = ['Halpha','r-band']
    names = labels
    names.reverse()
    
    plt.figure()

    for i,f in enumerate(fits):
        plt.subplot(2,1,i+1)
        #print("rmax = ",rmax)
        if rmax is not None:
            plt.xlim(0,rmax)
        #plt.gca().set_yscale('log')
        if (i == 1) & (log1Flag):
            plt.plot(phots[i]['sma_arcsec'],(phots[i]['sb'])/fits[i].amplitude.value,'bo',c=mycolors[0],label=labels[i])
            plt.plot(phots[i]['sma_arcsec'],fits[i]((phots[i]['sma_arcsec']))/fits[i].amplitude.value,label='SERSIC',c=mycolors[1])

        else:
            plt.plot(phots[i]['sma_arcsec'],phots[i]['sb']/fits[i].amplitude.value,'bo',c=mycolors[0],label=labels[i])
            plt.plot(phots[i]['sma_arcsec'],fits[i](phots[i]['sma_arcsec'])/fits[i].amplitude.value,label='SERSIC',c=mycolors[1])
        plt.axhline(ls=':',color='k')
        #plt.axhline(y=1,ls='--',color='k',alpha=.5)
        plt.text(0.5,0.8,f"{names[i]}: Re = {fits[i].r_eff.value:.1f}, n={f.n.value:.1f}",transform=plt.gca().transAxes)
        if i == 1:
            plt.xlabel("SMA (arcsec)",fontsize=16)
            plt.text(-.1,0.5,"Surface Brightness ",fontsize=16,transform=plt.gca().transAxes,rotation=90)
    size_ratio = hfit.r_eff.value/rfit.r_eff.value
    print(f"Ratio of Re({names[0]})/Re({names[1]}) = {size_ratio:.3}")
    return rfit, hfit

def plot_cog(rp,hp,rfit,hfit,rmax=None,labels = ['logMstar','logSFR']):
    # plot enclosed flux to see if Re is approx half light radius

    # this is using R band and CS Halpha
    import scipy.optimize
    from astropy.modeling import models, fitting, functional_models
    #Sersic1D
    plt.figure()



    # restrict to positive fluxes
    flag = rp['sb'] > 0
    rp = rp[flag]
    
   
    flag = hp['sb'] > 0
    hp = hp[flag]
    
    
    ###################################################
    # plot results
    ###################################################
    fits = [hfit,rfit]
    phots = [hp,rp]
    names = ['Halpha','r-band']
    names = labels
    names.reverse()
    
    plt.figure()
    radii50=[]
    for i,f in enumerate(fits):
        plt.subplot(2,1,i+1)
        if rmax is not None:
            plt.xlim(0,rmax)
        #plt.gca().set_yscale('log')

        plt.plot(phots[i]['sma_arcsec'],phots[i]['flux'],'bo',c=mycolors[0],label=labels[i])
        
        plt.axhline(y=0.5*np.max(phots[i]['flux']),ls=':',color='k')
        plt.axvline(x=fits[i].r_eff.value,ls='--',color='k')

    
        # plot R50, where the enclosed flux reaches the 0.5*total
        r50flag = phots[i]['flux'] > 0.5*np.max(phots[i]['flux'])
        r50 = phots[i]['sma_arcsec'][r50flag][0]
        radii50.append(r50)
        plt.axvline(x=r50,ls=':',color='k')                                     
        plt.text(0.05,0.8,f"{names[i]}: Re = {fits[i].r_eff.value:.1f}, n={f.n.value:.1f}, R50={r50:.1f}",transform=plt.gca().transAxes)
                                 
        if i == 1:

            plt.xlabel("SMA (arcsec)",fontsize=16)
            plt.text(-.15,1,"Enclosed Flux",fontsize=16,transform=plt.gca().transAxes,rotation=90,horizontalalignment='center')
    return radii50
    
def fit1profile(dirname='VFID5842-NGC5356-INT-20190206-p120',rmax=None):
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    #print("in fit1profile, rmax = ",rmax)

    os.chdir(dirname)

    rimage = dirname+'-logmstar-vr.fits'
    haimage = dirname+'-sfr-vr.fits'
    rphot = dirname+'-logmstar-vr_phot.fits'
    haphot = dirname+'-sfr-vr_phot.fits'


    rp = Table.read(rphot)
    hp = Table.read(haphot)
    # mstar image is logged, so correct for that
    print("in fit1profile, rmax = ",rmax)
    mfit,sfit = fit_profiles(rp,hp,rmax=rmax,labels=['logMstar','logSFR'],log1Flag=False)
    vfid = dirname.split('-')[0]
    plt.savefig(vfid+'-mstar-sfr-profiles.png')
        
    radii50 = plot_cog(rp,hp,mfit,sfit,rmax=rmax,labels=['logMstar','logSFR'])
    plt.savefig(vfid+'-mstar-sfr-cog.png')    
    # use this to run on R and CS Halpha
    rphot = dirname+'-R_phot.fits'
    haphot = dirname+'-CS-gr_phot.fits'

    rp = Table.read(rphot)
    hp = Table.read(haphot)
    rfit,hfit = fit_profiles(rp,hp,rmax=rmax,labels=['r','halpha'])
    plt.savefig(vfid+'-r-halpha-profiles.png')    
    return mfit,sfit,rfit,hfit


def plot_sfr_indicators(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],cbfrac=.08,cbaspect=20,clevels=[4,5],contourFlag=False):
    #%matplotlib inline
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name


    cwd = os.getcwd()
    os.chdir(dirname)
    massim = dirname+"-logmstar-vr.fits"
    sfrim = dirname+"-sfr-vr.fits"
    ssfrim = dirname+"-ssfr.fits"
    mask = dirname+'-R-mask.fits'
    vfid = dirname.split('-')[0]
    print(dirname+'/galex/*nuv*.fits')
    print("sfrim = ",sfrim)
    try:
        nuvim = glob.glob('galex/*nuv*.fits')[0]
    except IndexError:
        nuvim = None
        print("problem getting nuv image ",dirname)
    # look for coadded image first
    t = glob.glob('unwise/*w3-coadd.fits')
    if len(t) > 0:
        w3im = t[0]
    else:
        t = glob.glob(dirname+'/unwise/*w3-img-m.fits')
        w3im = t[0]
        
        
    titles = ['NUV',r'$H\alpha$',r'$WISE \ 12\mu m$']
    #vmin = [2,0,-11.5]
    #vmax = [6,1.e-6,-9]
    allim = [nuvim,sfrim,w3im]
    
    print(nuvim,w3im)
    
    if xmin is not None:
        # get ra  and dec corresponding to xmin,ymin xmax,ymax from sfr image
        # plot the legacy image in panel 1
        xcoords = np.array([xmin,xmax])
        ycoords = np.array([ymin,ymax])

        header = fits.getheader(sfrim)
        wcs = WCS(header)
        sky = wcs.pixel_to_world(xcoords,ycoords)
        print("skycoords = ",sky)
        
    plt.figure(figsize=(figsize[0],figsize[1]))

    plt.subplots_adjust(wspace=0.1)
    maskdat = fits.getdata(mask)

    allax = []
    for i, im in enumerate(allim):
        plt.subplot(1,4,i+2)
        dat = fits.getdata(im)
        
        if dat.shape == maskdat.shape:
            mdat = np.ma.array(dat,mask=maskdat)
        else:
            mdat = dat
        if xmin is None:
            mdat = mdat
        #elif i == 0: # galex images don't have WCS?
        #    mdat = mdat
        elif i in [0,2]: # convert coords for galex and unwise
            print(i,im)
            imwcs = WCS(fits.getheader(im))
            # set limits in ra,dec
            x,y = imwcs.world_to_pixel(sky)
            print(x,y)
            print(int(np.rint(y[0])),int(np.rint(y[1])),int(np.rint(x[0])),int(np.rint(x[1])))
            print(mdat.shape)
            mdat = mdat[int(np.rint(y[0])):int(np.rint(y[1])),int(np.rint(x[0])):int(np.rint(x[1]))] 
        else:
            mdat = mdat[ymin:ymax,xmin:xmax]
            
            
        #plt.imshow(mdat,origin='lower',interpolation='nearest')
        #if i == 2:
        #    plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],origin='lower',interpolation='nearest')
        #else:
        #    display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
        display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
        
        if not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        plt.title(titles[i],fontsize=20)
        plt.colorbar(fraction=cbfrac,aspect=cbaspect)
        if i == 1:
            t = dirname.split('-')
            plt.xlabel(t[0]+' '+t[1],fontsize=20)
        # plot contours from mass
        allax.append(plt.gca())
    # plot the legacy image in panel 1
    xcoords = np.array([xmin,xmax])
    ycoords = np.array([ymin,ymax])
    
    # get ramin,ramax and decmin,decmax from SFR image
    sfrim = dirname+"-sfr-vr.fits"
    header = fits.getheader(sfrim)
    wcs = WCS(header)
    sky = wcs.pixel_to_world(xcoords,ycoords)
    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band
    plt.subplot(1,4,1,projection=legwcs)
    plt.imshow(jpeg_data)
    # set limits in ra,dec
    x,y = legwcs.world_to_pixel(sky)
    # convert ramin,ramax and decmin,decmax to (x,y)
    plt.axis([x[0],x[1],y[0],y[1]])
    #print(x[0],x[1],y[0],y[1])
    if contourFlag:
        # get contours from logmstar image
        hdu = fits.open(massim)
        contour_data = hdu[0].data
        contour_header = hdu[0].header
        hdu.close()
        mcontour_data = np.ma.array(contour_data,mask=maskdat)
        for ax in allax:
            ax.contour((mcontour_data[ymin:ymax,xmin:xmax]),levels=clevels, colors='k',linestyles='-',linewidths=1)#,transform=ax.get_transform(WCS(contour_header))
    #plt.contour(contour_data[y[0]:y[1],x[0]:x[1]],levels = [4,5,6],colors='c',alpha=0.5)
    plt.title("Legacy grz",fontsize=20)
    # read in the
    plt.savefig(os.path.join(plotdir,dirname)+'-sfr-indicators.png')
    plt.show()

    os.chdir(cwd)


class grouptables(vtables):
    def get_group_members(self):
        self.groupMembs = v.kt['PGC1'] == 49547
        print(f"number of Kourkchi group members = {np.sum(self.groupMembs)}")
    def write_latex_table(self):

        flag = self.groupMembs & v.main['HAobsflag']
        col_names = ['VFID','NED Name','vr','logMstar','logSFR','logsSFR','H2 def','HI def', 'HI def Bos']
        col_formats={'logMstar': '%5.2f',\
                     'logSFR': '%5.2f',\
                     'logsSFR': '%5.2f',\
                     'H2 def': '%5.2f',\
                     'HI def': '%5.2f',\
                     'HI def Bos': '%5.2f'}
        latexdict={'preamble': r'\begin{center}',\
                   'tablefoot': r'\end{center}',\
                   'tabletype': 'table*',\
                   'header_start': '\\hline \\hline',\
                   'header_end': '\\hline',\
                   'data_end': '\\hline',\
                   'caption': 'NGC~5364 Group Members within the \\ha \\ footprint \\label{tab:sample}'}
        paperTab = Table([self.main['VFID'],self.main['NEDname'],self.main['vr'],self.magphys['logMstar_med'],self.magphys['logSFR_med'],self.magphys['logsSFR_med'],self.paper1['H2def'],self.paper1['HIdef'],self.a100['HIdef_bos']])[flag]
        paperTab.write(plotdir+'/NGC5364_tab1.tex',format='latex',names=col_names,formats=col_formats,latexdict=latexdict,\
                       fill_values=[(ascii.masked,'\\nodata')])#,(np.ma.masked,'no data')])
        #paperTab.write(plotdir+'/NGC5364_tab1.tex',format='latex',names=col_names,formats=col_formats,latexdict=ascii.latex.latexdicts['ApJ'])
        self.paperTab = paperTab
        pass


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v2/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_v2_', help = 'prefix for tables; default is vf_v2')                               
    args = parser.parse_args()

    if args.tabledir.startswith('/home/rfinn/'):
        homedir = os.getenv("HOME")
        args.tabledir = args.tabledir.replace('/home/rfinn',homedir)
    v = grouptables(args.tabledir,args.tableprefix) 
    v.read_all()

