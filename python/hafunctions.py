#!/usr/bin/env

import numpy as np

import os
homedir = os.getenv("HOME")

import sys
sys.path.append(homedir+'/github/HalphaImaging/python3/')
import plot_cutouts_ha as pch
sys.path.append(homedir+'/github/Virgo/programs/')
from readtablesv2 import vtables


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
import matplotlib.ticker as ticker

from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
from astropy import wcs
from astropy.table import Table
from astropy import convolution
from astropy.nddata import Cutout2D
from astropy import stats as apstats
from astropy import units as u

from scipy.stats import spearmanr

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from PIL import Image

import glob


import warnings
warnings.filterwarnings('ignore')

UNWISE_PIXSCALE = 2.75

# why am I using 1 arcsec for the legacy images???
LEGACY_PIXSCALE = 0.262


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plotdir = homedir+'/research/Virgo/plots/halpha/'






#zoom_coords_INT = {'VFID5889':[100,1300,850,2100],\
#                   'VFID5851':[255,1500,740,2120],\
#                  'VFID5855':[820,1200,700,1400],\
#                   'VFID5842':[900,1300,700,1450],\
#                   'VFID5859':[50,150,30,145],\
#                   'VFID5892':[240,780,250,690],\
#                   #'VFID6018':[80,580,70,320],\
#                   #'VFID6033':COdir+'n5566-co10-mean.fits',\
#                   #'VFID6091':COdir+'n5577-co10-mean.fits',\
#                   #'VFID6362':COdir+'u9661-co10-mean.fits',\                   
#                   }


# new coords to fit HI contours
zoom_coords_INT = {'VFID5889':[100,1300,850,2100],\
                   'VFID5709':[550,1050,650,950],\
                   'VFID5851':[255,1500,740,2120],\
                   'VFID5855':[730,1290,640,1460],\
                   'VFID5842':[800,1400,700,1450],\
                   'VFID5859':[1,175,1,175],\
                   'VFID5892':[240,780,250,690],\
                   'VFID5879':[1,207,1,207],\
                   'VFID5844':[1,194,1,194],\
                   'VFID6018':[500,1300,550,1200],\
                   'VFID6033':[1600,2750,1500,2800],\
                   'VFID6091':[900,1600,1020,1520],\
                   'VFID6362':[150,240,150,240],\
                       }

CO_zoom_coords_INT = {'VFID5889':[100,1300,850,2100],\
                   'VFID5709':[550,1050,650,950],\
                   'VFID5851':[255,1500,740,2120],\
                   'VFID5855':[830,1190,690,1410],\
                   'VFID5842':[920,1250,730,1430],\
                   'VFID5859':[1,175,1,175],\
                   'VFID5892':[240,780,250,690],\
                   'VFID5879':[1,207,1,207],\
                   'VFID5844':[1,194,1,194],\
                   'VFID6018':[520,1230,620,1100],\
                   'VFID6033':[1800,2500,1700,2600],\
                   'VFID6091':[980,1520,1080,1470],\
                   'VFID6362':[70,320,70,320],\
                       }

cofigsize = {
            'VFID5855':[10,10],\
            'VFID5842':[9,10],\
            'VFID5709':[10,4.75],\
            'VFID6018':[10,5.2],\
            'VFID6033':[10,8],\
            'VFID6091':[10,5.25],\
            'VFID6362':[10,6.5],\
           }

allratio90 = {'VFID5859':.6,\
              'VFID5892':.35,\
              'VFID5855':1.3,\
              'VFID5842':.6,\
              'VFID5889':1.11,\
              'VFID5851':.11,\
              'VFID5844':0,\
              'VFID5879':0}


# HI masses from Gautam, sent over slack on 6/3/2024
meerkat_HI_mass = {'VFID5859':7.84,\
                   'VFID5892':7.65,\
                   'VFID5855':8.98,\
                   'VFID5842':8.79,\
                   'VFID5889':9.43,\
                   'VFID5851':np.nan,\
                   'VFID5844':np.nan,\
                   'VFID5879':np.nan}

#'VFID5851': 
HIdir = homedir+'/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/'
HIdir2 = homedir+'/research/Virgo/MeerKAT/final_24Jan2025/'
HI_file = {'VFID5709': HIdir2+'J1406_p_0601_lw00_final_1_mom0.fits',\
               'VFID6033': HIdir2+'J1420_p_0336_lw00_final_1_mom0.fits',\
               'VFID6091': HIdir2+'J1420_p_0336_lw00_final_3_mom0.fits',\
               'VFID6018': HIdir2+'J1420_p_0336_lw00_final_6_mom0.fits',\
               'VFID6362': HIdir2+'J1502_p_0150_lw00_final_2_mom0.fits',\
               'VFID5889':HIdir+'ngc5363_ngc5364_ngc5360.fits',\
               'VFID5851':HIdir+'J1355_fin_lw05_bpcorr_2_mom0.fits',\
               #'VFID5851':None, \
               'VFID5855':HIdir+'J1355_fin_lw05_bpcorr_6_mom0.fits',\
               'VFID5842':HIdir+'J1355_fin_lw05_bpcorr_5_mom0.fits',\
               'VFID5859':HIdir+'J1355_p_0512_lw05_final_15_mom0.fits',\
               'VFID5892':HIdir+'ngc5363_ngc5364_ngc5360.fits',\
               'VFID5879':None,\
               'VFID5844':None 
}


COdir = homedir+'/research/Virgo/alma/combes_data/'
COveldir = homedir+'/research/Virgo/alma/combes_data/masked_velocity_maps/'
COwidthdir = homedir+'/research/Virgo/alma/combes_data/masked_velocity_width_maps/'
COmaskdir = homedir+'/research/Virgo/alma/combes_data/masked/'
CO_file = {'VFID5889':None,\
           'VFID5851':None,\
           #'VFID5851':None, \
           'VFID5709':COdir+'n5470-co10-mean.fits',\
           'VFID5855':COdir+'n5348-co10-mean.fits',\
           #'VFID5855':COmaskdir+'n5348_intensity_COmap_masked.fits',\
           #'VFID5855':COmaskdir+'n5348_intensity_COmap_masked.fits',\
           'VFID5842':COdir+'n5356-co10-mean.fits',\
           #'VFID5842':COmaskdir+'n5356_intensity_COmap_masked.fits',\
           'VFID5859':None,\
           'VFID5892':None,\
           'VFID5879':None,\
           'VFID5844':None,\
           'VFID6018':COdir+'n5560-co10-mean.fits',\
           'VFID6033':COdir+'n5566-co10-mean.fits',\
           'VFID6091':COdir+'n5577-co10-mean.fits',\
           'VFID6362':COdir+'u9661-co10-mean.fits',\
           
}


alma_galaxies = ['VFID5709', 'VFID5855', 'VFID5842','VFID6018','VFID6033','VFID6091','VFID6362']


CO_mask = {'VFID5889':None,\
           'VFID5851':None,\
           #'VFID5851':None, \
           #'VFID5855':COmaskdir+'n5348_mask.fits',\ 
           'VFID5709':COmaskdir+'n5470_intensity_COmap_masked.fits',\
           'VFID5855':None,\
           'VFID5855':COmaskdir+'n5348_intensity_COmap_masked.fits',\
           #'VFID5842':None,\
           #'VFID5842':COmaskdir+'n5356-co10-mean.fits',\
           'VFID5842':COmaskdir+'n5356_intensity_COmap_masked.fits',\
           'VFID5859':None,\
           'VFID5892':None,\
           'VFID5879':None,\
           'VFID5844':None,\
           'VFID6018':COmaskdir+'n5560_intensity_COmap_masked.fits',\
           'VFID6033':COmaskdir+'n5566_intensity_COmap_masked.fits',\
           'VFID6091':COmaskdir+'n5577_intensity_COmap_masked.fits',\
           'VFID6362':COmaskdir+'u9661_intensity_COmap_masked.fits',\
}

CO_vel = { 'VFID5709':COveldir+'n5470_velocity_COmap_masked.fits',\
           'VFID5855':COveldir+'n5348_velocity_COmap_masked.fits',\
           'VFID5842':COveldir+'n5356_velocity_COmap_masked.fits',\
           'VFID6018':COveldir+'n5560_velocity_COmap_masked.fits',\
           'VFID6033':COveldir+'n5566_velocity_COmap_masked.fits',\
           'VFID6091':COveldir+'n5577_velocity_COmap_masked.fits',\
           'VFID6362':COveldir+'u9661_velocity_COmap_masked.fits'
}
CO_width = {'VFID5709':COwidthdir+'n5470_width_COmap_masked.fits',\
           'VFID5855':COwidthdir+'n5348_width_COmap_masked.fits',\
           'VFID5842':COwidthdir+'n5356_width_COmap_masked.fits',\
           'VFID6018':COwidthdir+'n5560_width_COmap_masked.fits',\
           'VFID6033':COwidthdir+'n5566_width_COmap_masked.fits',\
           'VFID6091':COwidthdir+'n5577_width_COmap_masked.fits',\
           'VFID6362':COwidthdir+'u9661_width_COmap_masked.fits'
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
           'VFID5892':.075,\
            'VFID5879':.05,\
           'VFID5844':.05            
           
          }



afigsize = {'VFID5889':[14.5,4],\
            'VFID5851':[14,4],\
            'VFID5855':[14,5.2],\
            'VFID5842':[14,5.2],\
            'VFID5859':[14,3.5],\
            'VFID5892':[14,2.8],\
            'VFID5879':[14,3.5],\
            'VFID5844':[14,3.5],\
            'VFID6018':[14,3.5],\
            'VFID6033':[14,3.5],\
            'VFID6091':[14,3.5],\
           }
           
alevels = {'VFID5889':[3.5],\
           'VFID5851':[4.],\
           'VFID5842':[3.9],\
           'VFID5855':[3.9],\
           'VFID5842':[3.9],\
           'VFID5859':[3.5],\
           'VFID5892':[3.5],\
           'VFID5879':[3.75],\
           'VFID5844':[3.75],\
           'VFID6018':[4.2],\
           'VFID6033':[4.5],\
           'VFID6091':[4.2],\
           'VFID6362':[3.9]
           }



HIfiles = {'VFID5859':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_5_mom0.fits',\
           'VFID5855':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_6_mom0.fits',\
           'VFID5851':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_2_mom0.fits',\
           }

###########################################
## FUNCTIONS
###########################################

def add_scale(ax,vr=1100,pscale=.331,barsize=5.,color='k',fontsize=16,xscale=.09,yscale=.9,dytext=.08):
    '''
    ax = axis for drawing on
    vr = cosmic recession velocity
    pscale = pixel scale, in arcsec/pix
    barsize = size of marker in kpc
    fontisize = fontsize for the text
    '''

    from astropy.cosmology import WMAP9 as cosmo
        
    H0 = cosmo.H0.value

    z = vr/3.e5
    # get Mpc/radians
    add = cosmo.angular_diameter_distance(z)

    # convert to kpc/arcsec
    
    add_kpc_arcsec = add.value*1000*np.pi/(180*3600)
    
    # convert barsize from kpc to arcsec
    barsize_arcsec = barsize/add_kpc_arcsec
    
    # convert barsize from arcsec to pixels
    barsize_pixels = barsize_arcsec/pscale
    
    # get size of image
    x1,x2 = plt.gca().get_xlim()
    xline1 = x1 + xscale*(x2-x1)
    y1,y2 = plt.gca().get_ylim()
    yline1 = y1 + yscale*(y2-y1)
    # for size label
    xtext = xline1 + 0.5*barsize_pixels
    ytext = yline1 - dytext*(y2-y1)
    
    # set up arrays for reference line
    xbar = np.array([xline1,xline1+barsize_pixels])
    ybar = np.array([yline1,yline1])  
    # plot line
    plt.plot(xbar,ybar,'r-',lw=2,color=color)
    
    # label line
    plt.text(xtext,ytext,'{:.0f} kpc'.format(barsize),horizontalalignment='center',fontsize=fontsize,c=color)
    #plt.arrow(xline1,yline1,barsize_pixels,0)
    # draw line to show size of barsize


def convert_alma_header(inheader,outheader):
    """ remove all the header keywords that relate to 3rd axis """
    
    delcards = ['CTYPE','CRVAL','CDELT','CRPIX','CROTA','CUNIT']

    #for c in outheader.cards:
    for c in delcards:
        for i in range(2):
            ci = c+str(i+1)
            outheader.set(ci,inheader[ci])
    
    #outheader['NAXIS'] = 2
    #outheader['Number of WCS axes'] = 2
    return outheader


def convert_alma3d_alma2d(alma3d_image,nhdu=0):
    hdu = fits.open(alma3d_image)

   

    outdata = hdu[0].data[nhdu]

    newhdu = fits.PrimaryHDU(outdata)
    
    outheader = convert_alma_header(hdu[0].header,newhdu.header)
    newhdu.header = outheader
    return outdata, outheader


def get_cutouts(ra,dec,vfid,imsize):
    """
    imsize = size of cutout in arcsec

    """
    imsize_pixels = np.int(np.round(imsize/LEGACY_PIXSCALE))
    get_legacy_images(ra,dec,galid='VFID0',pixscale=0.262,imsize=imsize_pixels,band='r',makeplots=False,subfolder=None)


    
    get_wise_images(ra,dec,size=imsize)

    #get_unwise_image(ra,dec,galid=vfid,pixscale=2.75,imsize='60',bands='34',makeplots=False,subfolder=None):
    get_galex_images(ra,dec,size=imsize)    
    
def get_kpc_per_arcsec(vr):
    """
    PARAMS:
    * vr = recession velocity to use in distance calculation

    RETURN:
    * kpc_arcsec : physical scale corresponding to 1 arcsec

    """
    z = vr/3.e5
    
    from astropy.cosmology import WMAP9 as cosmo
    # get Mpc/radians
    add = cosmo.angular_diameter_distance(z)

    # convert to kpc/arcsec
    add_kpc_arcsec = add.value*1000*np.pi/(180*3600)
    
    return add_kpc_arcsec

def get_kpc_per_pixel(imwcs,vr):
    """
    PARAMS:
    * imwcs :  WCS for image
    * vr : recession velocity to use in angular diameter distance

    RETURN:
    * kpc_per_pixel : kpc per pixel at distance of H0*vr

    """
    pscale = wcs.utils.proj_plane_pixel_scales(imwcs) # in deg per pixel
    arcsec_per_pixel = pscale[0]*3600 #arcsec per pixel

    kpc_per_arcsec = get_kpc_per_arcsec(vr) # kpc/arcsec

    return arcsec_per_pixel*kpc_per_arcsec
    


def plot_HI_contours(ax,HIfilename,xmin=None,xmax=None,ymin=None,ymax=None,levels=None,color='white',addbeam=False,ncontour=5):
    """ plot HI contours, given current axis + reference image header """
    # borrowing from alma 2023 proposal
    #ncontour=3

    hdu = fits.open(HIfilename)[0]
    HI_WCS = WCS(hdu.header)
    #ncontour = np.array([1,4,15,29])
    #ncontour
    #ncontour = np.array([0,1,2.25,3.5,4])
    print("levels = ",levels)
    if levels is None:
        #levels = 3**np.arange(ncontour+1)+1
        levels = 3**np.arange(1,ncontour+1)+1
    #levels = np.append(levels,43)
    #levels = np.append(levels,50)
    #levels = np.append(levels,80)        
    ax.contour(hdu.data,transform=ax.get_transform(HI_WCS),levels=levels,colors=color,alpha=.4,lw=1)

def plot_HI_beam(ax,HIfilename,hostim_header,color='white',expandBox=False):
    ###
    # get size of the host image
    ###
    naxis1 = hostim_header['NAXIS1']
    naxis2 = hostim_header['NAXIS2']
    pscale = wcs.utils.proj_plane_pixel_scales(WCS(hostim_header))*u.deg # in deg
    image_x = naxis1*pscale[0]
    image_y = naxis2*pscale[1] 

    ###
    # get dimensions of host image
    ###
    hdu = fits.open(HIfilename)
    HI_WCS = WCS(hdu[0].header)
    PA = hdu[0].header['BPA']#*u.deg
    a = hdu[0].header['BMAJ']*u.deg
    b = hdu[0].header['BMIN']*u.deg
    hdu.close()
    
    #patch_height = (cube_params['bmaj'] / opt_view).decompose()
    #patch_width = (cube_params['bmin'] / opt_view).decompose()
    patch_height = a/image_x
    patch_width = b/image_y    
    patch = {'width': patch_width, 'height': patch_height}
    print(a,image_x,patch['height'], patch['width'], PA)
    ax.add_patch(Ellipse((0.9, 0.9), height=patch['height'], width=patch['width'], angle=PA,\
                         transform=ax.transAxes, edgecolor=color,facecolor=color, linewidth=1,alpha=.8))    

    # add white box to highlight position of beam
    dx = 0.1
    dy = 0.1
    if expandBox:
        dx=.25
        dy=.25
    xc = .9
    yc = 0.9
    rect = Rectangle((xc-0.5*dx,yc-0.5*dy),dx,dy,transform=ax.transAxes,edgecolor=color,facecolor='None',lw=1,alpha=.8)
    ax.add_patch(rect)

def plot_CO_beam(ax,header,frac_loc_mean=.1,boxsize_x=10,boxsize_y=10,maptype='continuum'):
    """
    ax = ax to make plot
    header = CO header file
    #maptype = continuum, intensity
    """
    PA = header['BPA'] #(deg) position angle, positive PA is from the N to the E (anticlock wise)
    if PA > 180 and PA <360:
        PA = PA-360.
    if PA < -180 and PA > -360:
        PA = 360.+PA #I redefine PA to be <180 in case PA in (-360,-180)
    if PA>=-180 and PA <=0 :
        angle = -90.-PA
    if PA>0 and PA<= 180. :
        angle = 90.-PA
    a = header['BMAJ']*60.*60. #major axis (arcsec) 
    b = header['BMIN']*60.*60. #minor axis (arcsec) 
    x0 = frac_loc_beam*boxsize_x
    y0= -frac_loc_beam*boxsize_y 
    #if 'continuum' in maptype or 'intensity' in maptype:
    #    _edgecolor = 'white'
    #else:
    _edgecolor = 'black'
    ell = Ellipse((x0,y0),a,b,angle,fill=False,hatch='////',edgecolor=_edgecolor)#,color='black')
    ell.set_clip_box(ax.bbox) 
    ell.set_alpha = 0.1
    
def plot_CO_contours(ax,COfilename,mask=None,xmin=None,xmax=None,ymin=None,ymax=None,levels=None,color='white',addbeam=False,ncontour=3):
    """ plot CO contours, given current axis + reference image header """
    # borrowing from alma 2023 proposal
    #ncontour=3
    
    #hdu = fits.open(COfilename)
    try:
        codata, coheader = convert_alma3d_alma2d(COfilename)
    except KeyError:
        codata, coheader = fits.getdata(COfilename, header=True)
    #print(coheader)
    CO_WCS = WCS(coheader)

    # checking mask - assuming these are now the most recent mask file from gianluca, not an actual mask
    if mask is not None:
        print()
        print("FYI: using a mask for CO data")
        # read in the mask data
        dat = fits.getdata(mask)

        # gianluca

        #has an image with good values are 1
        # and bad values are nan
        
        # create a masked array
        if mask.find('intensity_COmap') > -1:
            codata = dat
        else:
            codata = codata * dat
    
    #ncontour = np.array([1,4,15,29])
    #ncontour
    #print(CO_WCS)
    ncontour = np.array([0,1,2.25,3.5,4])
    ncontour = np.array([0,1,2.25,3.5,4])


    if levels is None:
        #levels = 3**np.arange(ncontour+1)+1
        levels = np.linspace(.1,1.4,3)
    #levels = np.append(levels,43)
    #levels = np.append(levels,50)
    #levels = np.append(levels,80)
    ax.contour(codata,transform=ax.get_transform(CO_WCS),slices=('x', 'y', 1),levels=levels,colors=color,alpha=.7,lw=1)


def plot_CO_image(ax,COfilename,mask=None,xmin=None,xmax=None,ymin=None,ymax=None,levels=None,color='white',addbeam=False,ncontour=3):
    """ plot CO contours, given current axis + reference image header """
    # borrowing from alma 2023 proposal
    #ncontour=3
    
    #hdu = fits.open(COfilename)
    codata, coheader = convert_alma3d_alma2d(COfilename)
    #print(coheader)
    CO_WCS = WCS(coheader)
    #ncontour = np.array([1,4,15,29])
    #ncontour
    #print(CO_WCS)

    if mask is not None:
        print()
        print("FYI: using a mask for CO data")
        # read in the mask data
        dat = fits.getdata(mask)

        # gianluca has an image with good values are 1
        # and bad values are nan
        
        # create a masked array
        codata = codata * dat
    ax.imshow(codata,transform=ax.get_transform(CO_WCS))
    

def plot_INT_footprint(center_ra,center_dec,plotguide=False):
    #using full detector sizes for now because 
    detector_dra = 4100.*0.33/3600.*1.05 # 2154 pixels * 0.33"/pix, /3600 to get deg
    detector_ddec = 2048.*0.33/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    # draw footprint of chip 4
    rect= plt.Rectangle((center_ra-detector_dra/2.,center_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)
    # draw footprint of chip 3
    # assuming chip 3 is NORTH and a smidge WEST of chip 4
    offset_dec = detector_ddec+17./3600. # 17 arcsec gap in RA between 
    offset_ra = -9.5/3600. # 9.5 arcsec offset toward N
    rect= plt.Rectangle((center_ra+offset_ra-detector_dra/2.,center_dec+offset_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)

    # draw footprint of chip 1
    # assuming chip 1 is SOUTH and a smidge EAST of chip 4
    offset_dec = -1*detector_ddec-22.7/3600. # 17 arcsec gap in RA between 
    offset_ra = +3.18/3600. # 9.5 arcsec offset toward N
    rect= plt.Rectangle((center_ra+offset_ra-detector_dra/2.,center_dec+offset_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)

    # draw footprint of chip 2
    # assuming chip 2 is WEST of chip 4
    offset_dec = detector_ddec/2.-detector_dra-19.2/3600. # hard to explain
    offset_ra =  -.5*detector_dra-23.8/3600.# hard to explain
    # this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
    rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -1.*detector_ddec, detector_dra,fill=False, color='k',label='WFC Footprint')
    plt.gca().add_artist(rect)

    if plotguide:
        # adding guide camera
        offset_dec = -2*detector_ddec-(22.7+98.1)/3600. # hard to explain
        offset_ra =  detector_dra/2-(3.18+649.9)/3600.# hard to explain
        # this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
        rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -7./60., 7./60,fill=False, color='k')
        plt.gca().add_artist(rect)
    
def plot_sfr_contours(dirname,ax,legwcs=None,xmin=None,xmax=None,ymin=None,ymax=None,levels=None):
    from astropy import convolution
    topdir = homedir+'/research/Virgo-dev/cont-sub-gr/'+dirname+'/'
    sfrim = topdir+dirname+"-sfr-vr.fits"
    hdu = fits.open(sfrim)[0]
    
    mask = topdir+dirname+'-R-mask.fits'
    if mask is not None:
        maskdat = fits.getdata(mask)
    else:
        maskdat = np.zeros_like(hdu.data,'i')
    if dirname.find('5889') > -1:
        maskdat[0:600,:]= 1
    mdat = np.ma.array(hdu.data,mask=maskdat)
    
    hdu.data = np.ma.array(hdu.data,mask=maskdat)
    if levels is not None:
        mylevels = levels
    else:
        mylevels = [1e-5,3e-5]
    if xmin is None:
        contour_data = hdu.data
        contour_header = hdu.header
        contour_WCS = WCS(contour_header)
    else:
        xcenter = (xmin+xmax)//2
        ycenter = (ymin+ymax)//2
        xsize = xmax - xmin
        ysize = ymax - ymin
        position = (ycenter,xcenter)
        size = (ysize,xsize)
        print(position)
        cutout = Cutout2D(hdu.data,position,size,wcs=WCS(hdu.header))
        contour_data = cutout.data
        #print(contour_data.wcs.to_header())
        contour_WCS = cutout.wcs
        
        
        # alternate approach - transform extent 
        #contour_data = hdu.data
        #contour_header = hdu.header
        #contour_WCS = WCS(contour_header)
        #xcoords = np.array([xmin,xmax])
        #ycoords = np.array([ymin,ymax])


        #header = hdu.header 
        #wcs = WCS(header)
        #sky = wcs.pixel_to_world(xcoords,ycoords)

        #x,y = legwcs.world_to_pixel(sky)
        
        #xmin,xmax = x
        #ymin,ymax = y
    contour_data = convolution.convolve_fft(contour_data, convolution.Box2DKernel(5), allow_huge=True, nan_treatment='interpolate')
    ax.contour(contour_data,transform=ax.get_transform(contour_WCS),levels=mylevels,colors='fuchsia',lw=5)


    
def add_spine(ax,color='c'):
    from astropy.table import Table
    spinedir = homedir+'/research/Virgo/tables-north/v2/spines/'
    spine  = Table.read(spinedir+'filament_spine_VirgoIII.fits')
    ax.plot(spine['ra'],spine['dec'],'c--',color=color,transform=ax.get_transform('world'))#, ,label='Filament Spine')

def plot_spines(multicolor=True,colorone=None,color=None,legend=True):
    import glob
    sfiles = glob.glob(homedir+'/research/Virgo/tables-north/spines/filament*.fits')
    ncolor = 0
    for i,f in enumerate(sfiles):
        spine  = Table.read(f)
        if legend:
            mylabel = os.path.basename(f).replace('filament_spine_','').replace('.fits','').replace('_Filament','')
        else:
            mylabel = None
        
        if multicolor:
            plt.plot(spine['ra'],spine['dec'],c=mycolors[ncolor],label=mylabel,lw=3)
        
            ncolor += 1
            if ncolor > len(mycolors)-1:
                ncolor = 0
        else:
            if colorone is not None:
                if colorone in f:
                    plt.plot(spine['ra'],spine['dec'],c=color,label=mylabel,lw=3,alpha=.6)
                else:
                    plt.plot(spine['ra'],spine['dec'],c='0.5',label=mylabel,lw=3)
            else:
                plt.plot(spine['ra'],spine['dec'],c='0.5',label=mylabel,lw=3)
    
def update_ngc5348_mask():
    imname = ''
    hdu = fits.open()

    # mask out all columns with with x > 1290
    

    # mask out rows with y > 1680


    
def plot_sfr_mstar():

    plt.figure(figsize=(8,6))
    plt.subplots_adjust(bottom=.15)

    x = v.magphys['logMstar_med']
    y = v.magphys['logSFR_med']

    HIdef = v.a100['HIdef_bos']
    magflag = x > 4
    flag = magflag    
    HIflag = magflag & (v.a100['logMH'] > 0)
    
    ##################################
    # plot the main sequence
    ##################################    
    plt.plot(x[flag],y[flag],'k.',alpha=.05,label='All VFS Galaxies')

    ##################################
    # plot the main sequence
    ##################################    
    xmin = 7
    xmax = 11.1
    xline = np.linspace(xmin,xmax,100)
    yline = 0.8*xline-8.32
    plt.plot(xline,yline,'k-',label='MS (Conger+2024)',alpha=.5)
    plt.plot(xline,yline+0.3,'k--',alpha=.5)
    plt.plot(xline,yline-0.3,'k--',alpha=.5)
    
    
    # plot the group members
    nha = 0
    nnha = 0

    sizeratio90 = np.zeros(len(v.a100))
    vfids = [5859,5892,5855,5842,5889,5851,5844,5879]
    ratio90 = [0.6,0.35,1.3,0.6,1.11,0.11,0,0]
                  

    for i,vf in enumerate(vfids):
        sizeratio90[vf] = ratio90[i]
    
    for vf in np.arange(len(v.main))[groupMembs]:
        if vf == 5851:
            print(vf,v.main['VFID'][vf])
            xp = v.magphys['logMstar_best'][int(vf)]
            yp = v.magphys['logSFR_best'][int(vf)]


        else:
            print(vf,v.main['VFID'][vf])
            xp, yp = x[int(vf)],y[int(vf)]
        if v.main['HAobsflag'][vf]:
            if nha == 0:
                label = r"$\rm Group \ Members \ In \ H\alpha \ FOV$"
                nha += 1
            else:
                label = "_nolegend_"
            #plt.plot(xp,yp,'bs',markersize=12,c='magenta',label=label)
            # add scatter for those with HIdef values
            #if HIflag[vf]:
            plotflag = HIflag[vf]
            #print(f"HIdef[vf] = {HIdef[vf]}")
            # use R90 as the color
            vfid = 'VFID'+str(vf)
            mycolor = allratio90[vfid]
            print(f"vfid = {vfid}, ratio = {mycolor}")
            #print(f"HIdef[vf] = {HIdef[vf]}")            
            #plt.scatter(xp,yp,c=HIdef[vf],s=120,marker='s',label=label,vmin=0,vmax=1.3)
            plt.scatter(xp,yp,c=mycolor,s=160,marker='s',label=label,vmin=0.1,vmax=1.2,alpha=.9)
            # use HI def as the color
            #plt.scatter(xp,yp,c=HIdef[vf],s=120,marker='s',label=label,vmin=-0.5,vmax=0.8)
            #else:
            #    #plt.plot(xp,yp,'bs',markersize=12,c='0.5')
                
            plt.text(xp-.02,yp+.2,v.main['VFID'][vf],horizontalalignment='center',fontsize=14)
        else:
            
            if nnha == 0:
                label = r"$All \ members$"
                nnha += 1
            else:
                label = "_nolegend_"
            if HIflag[vf]:
                plotflag = HIflag[vf]
                print(f"HIdef[vf] = {HIdef[vf]}")
                #plt.scatter(xp,yp,c=HIdef[vf],s=120,marker='8',label=label,vmin=-0.5,vmax=0.8)
            else:
                pass
                #plt.plot(xp,yp,'bo',markersize=12,c='0.5',label=label)

            #plt.plot(xp,yp,'bo',c=mycolors[0],label=label)

    xp = x[v.groupMembs & ~v.main['HAobsflag']]
    yp = y[v.groupMembs & ~v.main['HAobsflag']]    
    plt.scatter(xp,yp,c='0.55',s=160,marker='s',alpha=.9,label=r'$\rm Group \ Members /No \ H\alpha$')    
    cb = plt.colorbar()
    # cb.set_label('HI Def',fontsize=16)
    cb.set_label('$R_{90}(SFR)/R_{90}(M_\star)$',fontsize=16)
    
    plt.xlabel(r"$\rm \log_{10}(M_\star/M_\odot)$",fontsize=22)
    plt.ylabel(r"$\rm \log_{10}(SFR/(M_\odot/yr))$",fontsize=22)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    

    plt.xlim(xmin+.2,xmax)
    plt.ylim(-5.1,.2)
    plt.legend()
    plt.savefig(plotdir+'/NGC5364-sfr-mstar.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/NGC5364-sfr-mstar.pdf',bbox_inches="tight")        


def plot_NUVr_mstar():

    plt.figure(figsize=(8,6))
    plt.subplots_adjust(bottom=.15)

    xall = v.magphys['logMstar_med']

    ##################################################
    # get NUV-r mag from JM's elliptical photometry
    # 
    # NUV = 22.5 - 2.5*log10(flux_NUV)
    # r = 22.5 - 2.5*log10(flux_r)
    # NUV - r = 2.5*log10(flux_r/flux_NUV)
    #
    ##################################################    
    yall = 2.5*np.log10(v.ephot['FLUX_AP06_R']/v.ephot['FLUX_AP06_NUV'])

    ##################################################
    # cut low-mass crap out of sample
    ##################################################    
    magflag = xall > 4
    nuvflag = v.ephot['FLUX_AP06_NUV'] > .1
    flag = magflag    & nuvflag
    
    ##################################################
    # plot the full VFS sample
    ##################################################
    plt.plot(xall[flag],yall[flag],'k.',alpha=.05,label='All VFS Galaxies')

    ##################################
    # plot the main sequence
    ##################################    
    #xmin = 7
    #xmax = 11.1
    #xline = np.linspace(xmin,xmax,100)
    #yline = 0.8*xline-8.32
    #plt.plot(xline,yline,'k-',label='MS (Conger+2024)',alpha=.5)
    #plt.plot(xline,yline+0.3,'k--',alpha=.5)
    #plt.plot(xline,yline-0.3,'k--',alpha=.5)
    
    ##################################################    
    # plot the group members
    ##################################################    
    nha = 0
    nnha = 0

    sizeratio90 = np.zeros(len(v.a100))
    vfids = [5859,5892,5855,5842,5889,5851,5844,5879]
    ratio90 = [0.6,0.35,1.3,0.6,1.11,0.11,0,0]
                  

    for i,vf in enumerate(vfids):
        sizeratio90[vf] = ratio90[i]
    
    for vf in np.arange(len(v.main))[groupMembs]:
        if vf == 5851:
            print(vf,v.main['VFID'][vf])
            xp = v.magphys['logMstar_best'][int(vf)]
            yp = yall[int(vf)]


        else:
            print(vf,v.main['VFID'][vf])
            xp, yp = xall[int(vf)],yall[int(vf)]
        if v.main['HAobsflag'][vf]:
            if nha == 0:
                label = r"$\rm Group \ Members \ In \ H\alpha \ FOV$"
                nha += 1
            else:
                label = "_nolegend_"
            vfid = 'VFID'+str(vf)
            mycolor = allratio90[vfid]
            print(f"vfid = {vfid}, ratio = {mycolor}")
            plt.scatter(xp,yp,c=mycolor,s=160,marker='s',label=label,vmin=0.1,vmax=1.2,alpha=.9)
                
            plt.text(xp-.02,yp+.2,v.main['VFID'][vf],horizontalalignment='center',fontsize=14)
        else:
            pass

    xp = xall[v.groupMembs & ~v.main['HAobsflag'] & flag]
    yp = yall[v.groupMembs & ~v.main['HAobsflag'] & flag]    
    plt.scatter(xp,yp,c='0.55',s=160,marker='s',alpha=.9,label=r'$\rm Group \ Members /No \ H\alpha$')    
    cb = plt.colorbar()
    # cb.set_label('HI Def',fontsize=16)
    cb.set_label('$R_{90}(SFR)/R_{90}(M_\star)$',fontsize=16)
    
    plt.xlabel(r"$\rm \log_{10}(M_\star/M_\odot)$",fontsize=22)
    plt.ylabel(r"$\rm NUV-r$",fontsize=22)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    

    plt.xlim(7.2,11)
    plt.ylim(0,7)
    plt.legend()
    plt.savefig(plotdir+'/NGC5364-NUVr-mstar.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/NGC5364-NUVr-mstar.pdf',bbox_inches="tight")        


def plot_HIdef_sizeratio(paper1=False):
    """
    galaxies with HI detections: VFID5859, VFID5892, VFID5855, VFID5842, VFID5889

    #VFID5851: R90_mstar=81.2, R90_sfr=9.3, ratio=0.11 
    #VFID5851: R50_mstar=20.6, R50_sfr=6.9, ratio=0.33

    VFID5842: R90_mstar=80.8, R90_sfr=48.5, ratio=0.60 
    VFID5842: R50_mstar=36.6, R50_sfr=31.1, ratio=0.85
    
    VFID5859: R90_mstar=14.8, R90_sfr=8.9, ratio=0.60 
    VFID5859: R50_mstar=8.1, R50_sfr=4.8, ratio=0.59

    VFID5892: R90_mstar=59.3, R90_sfr=21.0, ratio=0.35 
    VFID5892: R50_mstar=25.8, R50_sfr=8.4, ratio=0.33

    VFID5855: R90_mstar=106.5, R90_sfr=138.4, ratio=1.30 
    VFID5855: R50_mstar=44.5, R50_sfr=50.1, ratio=1.12

    VFID5889: R90_mstar=163.8, R90_sfr=181.2, ratio=1.11 
    VFID5889: R50_mstar=51.3, R50_sfr=112.7, ratio=2.19

    #VFID5879: R90_mstar=15.2, R90_sfr=2.3, ratio=0.15 
    #VFID5879: R50_mstar=7.5, R50_sfr=1.6, ratio=0.22

    #VFID5844: R90_mstar=17.3, R90_sfr=9.9, ratio=0.57 
    #VFID5844: R50_mstar=8.1, R50_sfr=6.9, ratio=0.86 

    """

    # this is super clunky but going with it for now...
    sizeratio90 = np.zeros(len(v.a100))

    vfids = [5859,5892,5855,5842,5889,5851]
    ratio90 = [0.6,0.35,1.3,0.6,1.11,0.11]


    for i,vf in enumerate(vfids):
        sizeratio90[vf] = ratio90[i]


    vfids = [5859,5892,5855,5842,5889]
    plotflag = np.zeros(len(v.main),'bool')        
    for i,vf in enumerate(vfids):
        plotflag[vf] = True
    # now plot HIdef vs sizeratio90

    plt.figure(figsize=(8,6))
    plt.subplots_adjust(bottom=.15)
    
    x = sizeratio90[plotflag]

    
    #y = v.a100['HIdef_bos'][plotflag]
    
    y = v.a100['HIdef_bos'][plotflag]
    if paper1:
        y = v.paper1['HIdef'][plotflag]
    c = v.magphys['logMstar_med'][plotflag]
    plt.scatter(x,y,c=c,s=120,marker='s')
    plt.xlabel(r"$\rm R_{90}(SFR)/R_{90}(M_\star)$",fontsize=22)
    plt.ylabel(r"$\rm HI \ Deficiency$",fontsize=22)

    cb = plt.colorbar()
    cb.set_label("$\log(M_\star/M_\odot)$",fontsize=16)

    for vf in vfids:
        xp = sizeratio90[vf]
        yp = v.a100['HIdef_bos'][vf]
        if paper1:
            yp = v.paper1['HIdef'][vf] 
        plt.text(xp-.01,yp+.05,v.main['VFID'][vf],horizontalalignment='center',fontsize=14)

    plt.ylim(-.5,1.02)
    if paper1:
        plt.ylim(-.5,1.2)
    plt.xlim(0.22,1.39)
    plt.axhline(ls='-',c='k',alpha=.5)
    plt.axhline(y=0.25,ls='--',c='k',alpha=.5)
    plt.axhline(y=-0.25,ls='--',c='k',alpha=.5)
    
    #######################################
    # calculate spearman rank coeff
    #######################################
    r,pvalue = spearmanr(x,y)
    print(f"Spearman rank correlation coeff = {r:.2f}, pvalue = {pvalue:.3f}")


    plotflag[v.main['VFID']=='VFID5859'] = False
    x = sizeratio90[plotflag]
    y = v.a100['HIdef_bos'][plotflag]
    r,pvalue = spearmanr(x,y)
    print()
    print("removing VFID5859")
    print(f"Spearman rank correlation coeff = {r:.2f}, pvalue = {pvalue:.3f}")
    
    # now plot HIdef vs sizeratio90
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)    

    plt.savefig(plotdir+'/NGC5364-HIdef-sizeratio90.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/NGC5364-HIdef-sizeratio90.pdf',bbox_inches="tight")     


def plot_MHI_Mstar_sizeratio(paper1=False):
    """
    galaxies with HI detections: VFID5859, VFID5892, VFID5855, VFID5842, VFID5889

    #VFID5851: R90_mstar=81.2, R90_sfr=9.3, ratio=0.11 
    #VFID5851: R50_mstar=20.6, R50_sfr=6.9, ratio=0.33

    VFID5842: R90_mstar=80.8, R90_sfr=48.5, ratio=0.60 
    VFID5842: R50_mstar=36.6, R50_sfr=31.1, ratio=0.85
    
    VFID5859: R90_mstar=14.8, R90_sfr=8.9, ratio=0.60 
    VFID5859: R50_mstar=8.1, R50_sfr=4.8, ratio=0.59

    VFID5892: R90_mstar=59.3, R90_sfr=21.0, ratio=0.35 
    VFID5892: R50_mstar=25.8, R50_sfr=8.4, ratio=0.33

    VFID5855: R90_mstar=106.5, R90_sfr=138.4, ratio=1.30 
    VFID5855: R50_mstar=44.5, R50_sfr=50.1, ratio=1.12

    VFID5889: R90_mstar=163.8, R90_sfr=181.2, ratio=1.11 
    VFID5889: R50_mstar=51.3, R50_sfr=112.7, ratio=2.19

    #VFID5879: R90_mstar=15.2, R90_sfr=2.3, ratio=0.15 
    #VFID5879: R50_mstar=7.5, R50_sfr=1.6, ratio=0.22

    #VFID5844: R90_mstar=17.3, R90_sfr=9.9, ratio=0.57 
    #VFID5844: R50_mstar=8.1, R50_sfr=6.9, ratio=0.86 

    """

    # this is super clunky but going with it for now...
    sizeratio90 = np.zeros(len(v.a100))

    vfids = [5859,5892,5855,5842,5889,5851]
    ratio90 = [0.6,0.35,1.3,0.6,1.11,0.11]


    for i,vf in enumerate(vfids):
        sizeratio90[vf] = ratio90[i]


    vfids = [5859,5892,5855,5842,5889]
    plotflag = np.zeros(len(v.main),'bool')        
    for i,vf in enumerate(vfids):
        plotflag[vf] = True
    # now plot HIdef vs sizeratio90

    #plt.figure(figsize=(8,6))
    #plt.subplots_adjust(bottom=.15)
    
    plt.figure(figsize=(8,6))
    plt.subplots_adjust(bottom=.15,right=.9)
    
    x = sizeratio90[plotflag]

    

    # updating to plot HImass / stellar mass
    yall = np.zeros(len(v.main))
    for i,vf in enumerate(vfids):
        logMHI = meerkat_HI_mass['VFID'+str(vf)]
        # convert masses to Vcosmic
        logMstar = v.magphys['logMstar_best'][vf] + 2*np.log10(v.env['Vcosmic'][vf]/v.main['vr'][vf])
        
        yall[vf] = logMHI - logMstar
        print(vf,logMHI,v.magphys['logMstar_best'][vf],logMstar,yall[vf])

        
    y = yall[plotflag]
    c = v.magphys['logMstar_med'][plotflag]
    plt.scatter(x,y,c=c,s=120,marker='s')
    plt.xlabel(r"$\rm R_{90}(SFR)/R_{90}(M_\star)$",fontsize=22)
    plt.ylabel(r"$\rm \log_{10}(M_{HI}/M_\star) $",fontsize=22)

    cb = plt.colorbar()
    cb.set_label("$\log(M_\star/M_\odot)$",fontsize=16)

    for vf in vfids:
        xp = sizeratio90[vf]
        yp = yall[vf]
        plt.text(xp-.01,yp+.08,v.main['VFID'][vf],horizontalalignment='center',fontsize=14)

    plt.ylim(-2,1.2)
    #if paper1:
    #    plt.ylim(-.5,1.2)
    plt.xlim(0.19,1.41)
    #plt.axhline(ls='-',c='k',alpha=.5)
    #plt.axhline(y=0.25,ls='--',c='k',alpha=.5)
    #plt.axhline(y=-0.25,ls='--',c='k',alpha=.5)
    
    #######################################
    # calculate spearman rank coeff
    #######################################
    r,pvalue = spearmanr(x,y)
    print(f"Spearman rank correlation coeff = {r:.2f}, pvalue = {pvalue:.3f}")


    plotflag[v.main['VFID']=='VFID5859'] = False
    x = sizeratio90[plotflag]
    y = yall[plotflag]
    r,pvalue = spearmanr(x,y)
    print()
    print("removing VFID5859")
    print(f"Spearman rank correlation coeff = {r:.2e}, pvalue = {pvalue:.3e}")
    
    # now plot HIdef vs sizeratio90
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)    

    plt.savefig(plotdir+'/NGC5364-MHI-Mstar-sizeratio90.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/NGC5364-MHI-Mstar-sizeratio90.pdf',bbox_inches="tight")     

def plot_MH2_Mstar_sizeratio(paper1=False):
    """
    galaxies with HI detections: VFID5859, VFID5892, VFID5855, VFID5842, VFID5889


    """

    # this is super clunky but going with it for now...
    sizeratio90 = np.zeros(len(v.a100))

    vfids = [5859,5892,5855,5842,5889,5851]
    ratio90 = [0.6,0.35,1.3,0.6,1.11,0.11]


    for i,vf in enumerate(vfids):
        sizeratio90[vf] = ratio90[i]


    #vfids = [5859,5892,5855,5842,5889]
    # updating for H2 measurements
    vfids = [5892,5855,5842,5889,5851]
    plotflag = np.zeros(len(v.main),'bool')        
    for i,vf in enumerate(vfids):
        plotflag[vf] = True
    # now plot HIdef vs sizeratio90

    #plt.figure(figsize=(8,6))
    #plt.subplots_adjust(bottom=.15)
    
    plt.figure(figsize=(8,6))
    plt.subplots_adjust(bottom=.15,right=.9)
    
    x = sizeratio90[plotflag]

    

    # updating to plot HImass / stellar mass
    yall = np.zeros(len(v.main))
    for i,vf in enumerate(vfids):
        logMH2 = np.log10(v.paper1['Mmol'][vf])
        # convert masses to Vcosmic
        logMstar = v.magphys['logMstar_best'][vf] + 2*np.log10(v.env['Vcosmic'][vf]/v.main['vr'][vf])
        
        yall[vf] = logMH2 - logMstar
        print(vf,logMH2,v.magphys['logMstar_best'][vf],logMstar,yall[vf])

        
    y = yall[plotflag]
    c = v.magphys['logMstar_med'][plotflag]
    plt.scatter(x,y,c=c,s=120,marker='s')
    plt.xlabel(r"$\rm R_{90}(SFR)/R_{90}(M_\star)$",fontsize=22)
    plt.ylabel(r"$\rm \log_{10}(M_{H2}/M_\star) $",fontsize=22)

    cb = plt.colorbar()
    cb.set_label("$\log(M_\star/M_\odot)$",fontsize=16)

    for vf in vfids:
        xp = sizeratio90[vf]
        yp = yall[vf]
        plt.text(xp-.01,yp+.08,v.main['VFID'][vf],horizontalalignment='center',fontsize=14)

    #plt.ylim(-2,1.2)
    #if paper1:
    #    plt.ylim(-.5,1.2)
    #plt.xlim(0.19,1.41)
    #plt.axhline(ls='-',c='k',alpha=.5)
    #plt.axhline(y=0.25,ls='--',c='k',alpha=.5)
    #plt.axhline(y=-0.25,ls='--',c='k',alpha=.5)
    
    #######################################
    # calculate spearman rank coeff
    #######################################
    r,pvalue = spearmanr(x,y)
    print(f"Spearman rank correlation coeff = {r:.2f}, pvalue = {pvalue:.3f}")


    plotflag[v.main['VFID']=='VFID5859'] = False
    x = sizeratio90[plotflag]
    y = yall[plotflag]
    r,pvalue = spearmanr(x,y)
    print()
    print("removing VFID5859")
    print(f"Spearman rank correlation coeff = {r:.2e}, pvalue = {pvalue:.3e}")
    
    # now plot HIdef vs sizeratio90
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)    

    plt.savefig(plotdir+'/NGC5364-MH2-Mstar-sizeratio90.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/NGC5364-MH2-Mstar-sizeratio90.pdf',bbox_inches="tight")     
    
def get_rad_fluxfrac_sersfunc(pfit,frac=0.9,verbose=False):
    from scipy.interpolate import interp1d
    N = 10
    x = np.linspace(0,N*pfit.r_eff.value,100*N)
    flux = pfit(x)*2*np.pi*x*(x[1]-x[0])
        
    integral2 = np.cumsum(flux)
    
    interp_profile = interp1d(integral2,x)
    if verbose:
        print(f"cumulative sum of profile = {integral2[-1]:.2e}")
    return interp_profile(frac*integral2[-1])

def get_rad_fluxfrac(pfit,frac=0.9,verbose=False):
    from scipy.interpolate import interp1d
    N = 10
    x = np.linspace(0,N*pfit.r_eff.value,100*N)
    flux = pfit(x)*2*np.pi*x*(x[1]-x[0])
        
    integral2 = np.cumsum(flux)
    
    interp_profile = interp1d(integral2,x)
    if verbose:
        print(f"cumulative sum of profile = {integral2[-1]:.2e}")
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
                   clevels=[4],contourFlag=True,zoomflag=True,COlevels=None):
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
    if ('INT' in dirname) or ('BOK' in dirname):
        try:
            xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
        except KeyError:
            print("no xmin,xmax in dictionary")
            t = fits.getdata(massim)
            xmin=1
            ymin=1
            ymax,xmax = t.shape
            if zoomflag:
                # zoom by factor of two
                xcenter = xmax//2
                ycenter = ymax//2
                xmin = xcenter-xmax//4
                xmax = xcenter + xmax//4
                xmin = xcenter-xmax//4
                xmax = xcenter + xmax//4            
    else:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
        except KeyError:
            print("no xmin,xmax in dictionary")
            t = fits.getdata(massim)
            xmin=1
            ymin=1
            ymax,xmax = t.shape
            if zoomflag:
                # zoom by factor of two
                xcenter = xmax//2
                ycenter = ymax//2
                xmin = xcenter-xmax//4
                xmax = xcenter + xmax//4
                xmin = xcenter-xmax//4
                xmax = xcenter + xmax//4            
        
    

    xcoords = np.array([xmin,xmax])
    ycoords = np.array([ymin,ymax])

    header = fits.getheader(sfrim)
    wcs = WCS(header)
    sky = wcs.pixel_to_world(xcoords,ycoords)
    #print("skycoords = ",sky)

    try:
        myfigsize=afigsize[vfid]
    except KeyError:
        myfigsize = [14,5]
    try:
        mycbfrac=acbfrac[vfid]
    except KeyError:
        mycbfrac = 0.08

    try:
        myclevels=alevels[vfid]
    except KeyError:
        myclevels = [4]            

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
   
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-ssfr.png',dpi=150,bbox_inches="tight")


    os.chdir(cwd)


def plot_mstar_sfr_profiles(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,harmax=None,\
                            Re_mstar=None,Re_sfr=None,R90_mstar=None,R90_sfr=None,logMstar=None,\
                            vr=None,\
                            cmap='magma_r',markGroupCenter=False,beamCorrected=False):
    """
    same plot as mstar_sfr, but swap out ssfr for radial profiles in the 4th panel

    PARAMS:
    vr = recession velocity to use for DA when converting flux/pixel to flux/kpc^2

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
    #titles = [r'$\log_{10}(M_\star/M_\odot)$',r'$H\alpha \ SFR$','log sSFR']
    titles = [r'$\log_{10}(\Sigma_\star)$',r'$\Sigma_{SFR}$','log sSFR']
    vmin = [2,0e-5,-11.5]
    vmax = [6,.6e-5,-9]
    allim = [massim,sfrim,ssfrim]
    allim = [massim,sfrim]    
    #cblabels = [r'$\rm \log_{10}(M_\star/M_\odot)/pixel$',r'$\rm SFR(M_\star/yr)/pixel$']
    cblabels = [r'$\rm \log_{10}(M_\star/M_\odot/kpc^2)$',r'$\rm SFR(M_\star/yr)/kpc^2$']
    
    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)
    if 'INT' in dirname:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
        except KeyError:
            # use the full image
            # read in massim
            data = fits.getdata(massim)
            ymax,xmax = data.shape
            xmin=1
            ymin=1

    else:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
        except KeyError:
            # use the full image
            # read in massim
            data = fits.getdata(massim)
            ymax,xmax = data.shape
            xmin=1
            ymin=1
    

    try:
        myfigsize=afigsize[vfid]
    except KeyError:
        myfigsize = [14.5,5]
    try:
        mycbfrac=acbfrac[vfid]
    except KeyError:
        mycbfrac=.08
    try:
        myclevels=alevels[vfid]
    except KeyError:
        myclevels=3.75

        
    fig = plt.figure(figsize=(myfigsize[0],myfigsize[1]))

    plt.subplots_adjust(wspace=0.01,bottom=.15)

    if 'VFID5892' in dirname:
        cbaspect = 10
        plt.subplots_adjust(wspace=0.01,bottom=.2)
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
        ax2=plt.subplot(1,4,i+2,projection=imwcs)


        v1 = vmin[i] 
        v2 = vmax[i]       
        if vr is not None:
            # scale data to convert to flux/kpc^2
            kpc_per_pixel = get_kpc_per_pixel(imwcs,vr)
            if i == 0:
                dat = np.log10(10**dat/(kpc_per_pixel**2))
                v1 = vmin[i]-2.*np.log10(kpc_per_pixel) 
                v2 = vmax[i]-2.*np.log10(kpc_per_pixel)           
                #myclevels = myclevels-2.*np.log10(kpc_per_pixel)
                print("\nmass contour",myclevels-2.*np.log10(kpc_per_pixel))
            else:
                dat = dat/(kpc_per_pixel**2)
                v1 = vmin[i]/kpc_per_pixel**2 
                v2 = vmax[i]/kpc_per_pixel**2
                #myclevels = myclevels/kpc_per_pixel**2
                
        mdat = np.ma.array(dat,mask=maskdat)
        #if xmin is None:
        #    mdat = mdat
        #else:
        #    mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            plt.imshow(mdat,origin='lower',interpolation='nearest',cmap=cmap,vmin=v1,vmax=v2)

        else:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            plt.imshow(mdat,cmap=cmap,vmin=v1,vmax=v2)#cmap='viridis'
        
            lon = ax2.coords[0]
            lat = ax2.coords[1]
            lon.set_ticklabel_visible(False)
            lon.set_ticks_visible(False)            
            lat.set_ticklabel_visible(False)
            lat.set_ticks_visible(False)            
        
        #plt.colorbar(fraction=mycbfrac,aspect=cbaspect)
        # plot contours from mass

        #############################################################
        # add HI contour to halpha image
        #############################################################
        if i == 1:
            # check if HI moment zero map is available
            if beamCorrected:
                # beam corrected version
                HIfilename = homedir+'/research/Virgo/MeerKAT/2024Jul12-data/J1355_p_0512_pbcorr_mom0.fits'

                # lower resolution version
                #HIfilename = homedir+'/research/Virgo/MeerKAT/2024Jul12-data/J1355_p_0512_tp30_pbcorr_mom0.fits'

            else:
                HIfilename = HI_file[vfid]
            if HIfilename is not None:
                print("HIfilename = ",HIfilename)
                plot_HI_contours(ax2,HIfilename,color='steelblue')
                if vfid == 'VFID5859':
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                else:
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 

        #############################################################
        # add contours from stellar mass image
        #############################################################    
        if contourFlag:
            ax = plt.gca()
            ax.contour(mcontour_data,levels=myclevels, colors='k',linestyles='-',linewidths=1,transform=ax.get_transform(contour_WCS))

                
        plt.title(titles[i],fontsize=18)



        if xmin is not None:
            plt.axis([xmin,xmax,ymin,ymax])
        if i in [0,1]:
            plt.xticks([],[])
            plt.yticks([],[])

        elif not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        
        allax.append(plt.gca())

        #if i in [1,2]:
        cbaxes = inset_axes(ax2, width="80%", height="3%", loc=8) 
        cb = plt.colorbar(cax=cbaxes, orientation='horizontal')
        cb.set_label(label=cblabels[i],fontsize=12)

    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band

    
    ax1 = plt.subplot(1,4,1,projection=legwcs)
    plt.imshow(jpeg_data)

    if xmin is not None:
        # plot the legacy image in panel 1
        xcoords = np.array([xmin,xmax])
        ycoords = np.array([ymin,ymax])
    
        # get ramin,ramax and decmin,decmax from SFR image
        sfrim = dirname+"-sfr-vr.fits"
        header = fits.getheader(sfrim)
        #print(header)
        imwcs = WCS(header)
        sky = imwcs.pixel_to_world(xcoords,ycoords)
    
        # set limits in ra,dec
        x,y = legwcs.world_to_pixel(sky)
        # convert ramin,ramax and decmin,decmax to (x,y)
        #print(sky)
        plt.axis([x[0],x[1],y[0],y[1]])
    t = dirname.split('-')
    #plt.text(.05,.02,t[0],fontsize=20,transform=plt.gca().transAxes,horizontalalignment='left',color='white')
    

    #plt.title(f"{t[0]} Legacy grz",fontsize=18)
    plt.title(f"{t[0]}",fontsize=18)    
    #plt.sca(ax1)
    plt.xlabel("RA (hr)",fontsize=14)
    plt.ylabel("DEC (deg)",fontsize=14)


    ###############################################
    # add a physical scale
    ###############################################
    pscale = wcs.utils.proj_plane_pixel_scales(legwcs) # in deg -> arcsec
    pscale_arcsec = pscale[0]*3600
    add_scale(ax1,vr=v.mw_vrcenter,pscale=pscale_arcsec,color='w',barsize=2)    

    

    #############################################################
    # add arrow showing direction to group center
    # TODO : figure out how to make this smaller in smaller images
    #############################################################
    print("VFID = ",vfid)

    if (vfid == 'VFID5855') | (vfid == 'VFID5859'):
        print('hi')
        # get coords of galaxy
        vflag = v.main['VFID'] == vfid
        print(np.sum(vflag))
        x = v.main['RA'][vflag][0]
        y = v.main['DEC'][vflag][0]      
        dRA = v.mw_RAcenter - x
        dDEC = v.mw_DECcenter - y        
        norm = np.sqrt(dRA**2+dDEC**2)

        if  (vfid == 'VFID5855'):
            scale = 0.02
            mylw = 0.05
            dx = dRA/norm*scale
            dy = dDEC/norm*scale

            plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),lw=mylw,alpha=.7)
        elif (vfid == 'VFID5859'):
            scale = 0.0025
            mylw=.5
            dx = dRA/norm*scale
            dy = dDEC/norm*scale
            width=.001
            plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),linewidth=mylw,alpha=.7,head_width=2.5*width,head_length=2.5*width)
        print(x,y)
        #plt.plot(x*u.deg,y*u.deg,'bo',markersize=20,color='c',transform=ax1.get_transform('world'))
        

    #print(x[0],x[1],y[0],y[1])
    #############################################################
    # add HI contour to legacy image
    #############################################################    
    # check if HI moment zero map is available
    #HIfilename = HI_file[vfid]
    if HIfilename is not None:
        print("HIfilename = ",HIfilename)
        plot_HI_contours(plt.gca(),HIfilename,color='lightsteelblue')
        smallgals = ['VFID5859','VFID5879','VFID5844']
        if vfid in smallgals:
            plot_HI_beam(ax1,HIfilename,hdu.header,color='steelblue',expandBox=True)
        else:
            plot_HI_beam(ax1,HIfilename,hdu.header,color='steelblue') 
        
        
    #############################################################
    # add stellar mass from magphys to the legacy image
    #############################################################    
    if logMstar is not None:
        print("adding logMstar = ",logMstar)
        plt.text(0.05,0.05,logMstar,fontsize=16,color='white',transform=plt.gca().transAxes,horizontalalignment='left')
    
    

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
    labels = ['$M_\star$',r'$H\alpha \ SFR$']    
    colors = ['navy','mediumvioletred']
    markers = ['s','o']
    for i,t in enumerate(tables):
        x = t['sma_arcsec']
        y = t['sb']
        yerr = t['sb_err']
        if (i == 1) & (('VFID5844' in dirname) | ('VFID5879' in dirname)):
            continue
        if ('VFID5859' in dirname) | ('VFID5892' in dirname):
            norm_factor = np.median(y[0:5])
        else:
            norm_factor = np.median(y[0:15])

        pflag = x < rmax
        if i == 1:
            if harmax is not None:
                pflag = x < harmax
        plt.plot(x[pflag],y[pflag]/norm_factor,'bo',c=colors[i],label=labels[i],marker=markers[i])
        #plt.fill_between(x,(y+yerr)/norm_factor,y2=(y-yerr)/norm_factor,color=mycolors[i],alpha=.2)
    # add R50 if it's provided
    if Re_mstar is not None:
        plt.axvline(x=Re_mstar,ls='--',lw=2,color=colors[0],label='$R_e$')#,label='$R_e(M_\star)$')
    if R90_mstar is not None:
        plt.axvline(x=R90_mstar,ls='-',lw=2,color=colors[0],label='$R_{90}$')#,label='$R_{90}(M_\star)$')

    if (('VFID5844' in dirname) | ('VFID5879' in dirname)):
        print("not adding halpha radii")
    else:
        if Re_sfr is not None:
            plt.axvline(x=Re_sfr,ls='--',color=colors[1])#,label='$R_e(SFR)$')
            if R90_sfr is not None:
                plt.axvline(x=R90_sfr,ls='-',color=colors[1])#,label='$R_{90}(SFR)$')
        
    plt.legend()#bbox_to_anchor=(1.02,0.95))
    plt.gca().set_yscale('log')
    plt.gca().yaxis.tick_right()
    plt.ylim(.003,3.5)
    #print("rmax = ",rmax)    
    if rmax is not None:
        
        plt.xlim(-1,rmax)

    #plt.ylim(-.1,1.5)
    plt.xlabel("SMA (arcsec)",fontsize=16)

    
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-profiles-4panel.png',dpi=150,bbox_inches="tight")


    os.chdir(cwd)


    return ax1


def plot_mstar_sfr_CO(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,\
                            logMstar=None,cmap='magma_r',markGroupCenter=False,COcolor='white',COlevels=None):
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
    #titles = [r'$\log_{10}(M_\star/M_\odot)$',r'$H\alpha \ SFR$','log sSFR']
    titles = [r'$\log_{10}(\Sigma_\star)$',r'$\Sigma_{SFR} \ w/HI$',r'$\Sigma_{SFR} \ w/CO$']
    vmin = [2,0e-5,0e-5]
    vmax = [6,.6e-5,.6e-5]
    allim = [massim,sfrim,ssfrim]
    #allim = [massim,sfrim]
    allim = [massim,sfrim,sfrim] # show halpha with HI and CO    
    cblabels = [r'$\rm \log_{10}(M_\star/M_\odot)/pixel$',r'$\rm SFR(M_\star/yr)/pixel$',r'$\rm SFR(M_\star/yr)/pixel$']
    
    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)
    if 'INT' in dirname:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
        except KeyError:
            # use the full image
            # read in massim
            data = fits.getdata(massim)
            ymax,xmax = data.shape
            xmin=1
            ymin=1

    else:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
        except KeyError:
            # use the full image
            # read in massim
            data = fits.getdata(massim)
            ymax,xmax = data.shape
            xmin=1
            ymin=1
    

    try:
        myfigsize=afigsize[vfid]
    except KeyError:
        myfigsize = [14,5]
    try:
        mycbfrac=acbfrac[vfid]
    except KeyError:
        mycbfrac = 0.08

    try:
        myclevels=alevels[vfid]
    except KeyError:
        myclevels = [4]            

    #myfigsize=afigsize[vfid]
    #mycbfrac=acbfrac[vfid]
    #myclevels=alevels[vfid]

        
    fig = plt.figure(figsize=(myfigsize[0],myfigsize[1]))

    plt.subplots_adjust(wspace=0.01,bottom=.15)

    if 'VFID5892' in dirname:
        cbaspect = 10
        plt.subplots_adjust(wspace=0.01,bottom=.2)
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
        #if i == 2:
        #    ax2=plt.subplot(1,4,i+2,projection=imwcs,slices=('x','y',1))
        #else:
        ax2=plt.subplot(1,4,i+2,projection=imwcs)
        
                        


        mdat = np.ma.array(dat,mask=maskdat)
        #if xmin is None:
        #    mdat = mdat
        #else:
        #    mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],origin='lower',interpolation='nearest',cmap='gray_r')


        else:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],cmap=cmap)#cmap='viridis'
        
        lon = ax2.coords[0]
        lat = ax2.coords[1]
        lon.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)            
        lat.set_ticklabel_visible(False)
        lat.set_ticks_visible(False)            
        
        #plt.colorbar(fraction=mycbfrac,aspect=cbaspect)
        # plot contours from mass

        #############################################################
        # add HI contour to halpha image
        #############################################################
        if i == 1:
            # check if HI moment zero map is available
            try:
                HIfilename = HI_file[vfid]
            except KeyError:
                HIfilename = None
            if HIfilename is not None:
                print("HIfilename = ",HIfilename)
                plot_HI_contours(ax2,HIfilename,color='steelblue')
                if vfid == 'VFID5859':
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                else:
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 

        #############################################################
        # add CO contour to halpha image
        #############################################################
        if i == 2:
            # check if HI moment zero map is available
            COfilename = CO_file[vfid]
            COmask_filename = CO_mask[vfid]
            #if COmask_filename is not None:
            #print("COfilename = ",COfilename)

            plot_CO_contours(ax2,COfilename,color=COcolor,mask=COmask_filename,levels=COlevels)
                #if vfid == 'VFID5859':
                #    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                #else:
                #    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 


                    
        #############################################################
        # add contours from stellar mass image
        #############################################################    
        if contourFlag:
            ax = plt.gca()
            ax.contour(mcontour_data,levels=myclevels, colors='k',linestyles='-',linewidths=1,transform=ax.get_transform(contour_WCS))

                
        plt.title(titles[i],fontsize=18)



        if xmin is not None:
            plt.axis([xmin,xmax,ymin,ymax])
        if i in [0,1]:
            plt.xticks([],[])
            plt.yticks([],[])

        elif not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        
        allax.append(plt.gca())

        #if i in [1,2]:
        cbaxes = inset_axes(ax2, width="80%", height="3%", loc=8) 
        cb = plt.colorbar(cax=cbaxes, orientation='horizontal')
        cb.set_label(label=cblabels[i],fontsize=12)

    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band

    
    ax1 = plt.subplot(1,4,1,projection=legwcs)
    plt.imshow(jpeg_data)

    #ax1.set_xlabel("RA",fontsize=16)
    #ax1.set_ylabel("DEC",fontsize=16)    

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
    

    #plt.title(f"{t[0]} Legacy grz",fontsize=18)
    plt.title(f"{t[0]}-{t[1]}",fontsize=18)    

    #############################################################
    # add arrow showing direction to group center
    # TODO : figure out how to make this smaller in smaller images
    #############################################################
    print("VFID = ",vfid)

    if (vfid == 'VFID5855') | (vfid == 'VFID5859'):
        print('hi')
        # get coords of galaxy
        vflag = v.main['VFID'] == vfid
        print(np.sum(vflag))
        x = v.main['RA'][vflag][0]
        y = v.main['DEC'][vflag][0]      
        dRA = v.mw_RAcenter - x
        dDEC = v.mw_DECcenter - y        
        norm = np.sqrt(dRA**2+dDEC**2)

        if  (vfid == 'VFID5855'):
            scale = 0.02
            mylw = 0.05
            dx = dRA/norm*scale
            dy = dDEC/norm*scale

            #plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),lw=mylw,alpha=.7)
        elif (vfid == 'VFID5859'):
            scale = 0.0025
            mylw=.5
            dx = dRA/norm*scale
            dy = dDEC/norm*scale
            width=.001
            plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),linewidth=mylw,alpha=.7,head_width=2.5*width,head_length=2.5*width)
        print(x,y)
        #plt.plot(x*u.deg,y*u.deg,'bo',markersize=20,color='c',transform=ax1.get_transform('world'))
        

    #print(x[0],x[1],y[0],y[1])
    #############################################################
    # add HI contour to legacy image
    #############################################################    
    # check if HI moment zero map is available
    try:
        HIfilename = HI_file[vfid]
    except KeyError:
        HIfilename = None
    if HIfilename is not None:
        print("HIfilename = ",HIfilename)
        plot_HI_contours(plt.gca(),HIfilename,color='lightsteelblue')
        plot_HI_beam(plt.gca(),HIfilename,fits.getheader(legacyr),color='steelblue') 
        
    #############################################################
    # add stellar mass from magphys to the legacy image
    #############################################################    
    if logMstar is not None:
        print("adding logMstar = ",logMstar)
        #plt.text(0.05,0.05,logMstar,fontsize=16,color='white',transform=plt.gca().transAxes,horizontalalignment='left')
    
    
    ax1.update({'xlabel': 'RA', 'ylabel': 'DEC'})


    plt.sca(ax1)
    plt.xlabel("RA (hr)",fontsize=14)
    plt.ylabel("DEC (deg)",fontsize=14)
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-profiles-CO.png',dpi=150,bbox_inches="tight")


    os.chdir(cwd)


    return ax1


def plot_mstar_sfr_COall(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,4],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,vr=None,\
                            logMstar=None,cmap='magma_r',markGroupCenter=False,COcolor='white',COlevels=None):
    """
    same plot as mstar_sfr, but swap out ssfr for radial profiles in the 4th panel

    1 - legacy
    2 - stellar mass
    3 - SFR
    4 - SFR w/CO
    5 - CO image
    6 - CO vel
    7 - CO width

    """
    #%matplotlib inline
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name


    cwd = os.getcwd()
    os.chdir(dirname)

    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)

    ###########################################
    # mstar, sfr images
    ###########################################
    massim = dirname+"-logmstar-vr.fits"
    sfrim = dirname+"-sfr-vr.fits"

    
    ssfrim = dirname+"-ssfr.fits"
    mask = dirname+'-R-mask.fits'

    ###########################################
    # CO Images
    ###########################################
    coim = CO_file[vfid]
    #        COmask_filename = CO_mask[vfid]
    covelim = CO_vel[vfid]
    cowidthim = CO_width[vfid]

    ###########################################
    # other plot params
    ###########################################
    #titles = [r'$\log_{10}(M_\star/M_\odot)$',r'$H\alpha \ SFR$','log sSFR']
    titles = [r'$\log_{10}(\Sigma_\star)$',r'$\Sigma_{SFR} \ w/HI$',r'$\Sigma_{SFR} \ w/CO$',\
                 '$CO~ Intensity$','$CO ~velocity$',r'$CO ~ \Delta v $']
                 
    vmin = [2,0e-5,0e-5,0,-200,0]
    vmax = [6,.6e-5,.6e-5,2,200,100]
    allim = [massim,sfrim,ssfrim]
    #allim = [massim,sfrim]
    allim = [massim,sfrim,sfrim,coim,covelim,cowidthim] # show halpha with HI and CO

    cblabels = [r'$\rm \log_{10}(M_\star/M_\odot/kpc^2)$',\
                    r'$\rm SFR(M_\star/yr)/kpc^2$',\
                    r'$\rm SFR(M_\star/yr)/kpc^2$',\
                    r'$\rm Jy~beam^{-1}~ km~s^{-1}$',\
                    r'$\rm km~s^{-1}$',r'$\rm km~s^{-1}$']
    
    try:
        xmin,xmax,ymin,ymax = CO_zoom_coords_INT[vfid]
    except KeyError:
        data = fits.getdata(massim)
        ymax,xmax = data.shape
        xmin=1
        ymin=1

#        try:
#            xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
#        except KeyError:
#            # use the full image
#            # read in massim
#            data = fits.getdata(massim)
#            ymax,xmax = data.shape
#            xmin=1
#            ymin=1
    

    myfigsize=cofigsize[vfid]
    try:
        mycbfrac=acbfrac[vfid]
    except KeyError:
        mycbfrac = 0.08

    try:
        myclevels=alevels[vfid]
    except KeyError:
        myclevels = [4]            

    #myfigsize=afigsize[vfid]
    #mycbfrac=acbfrac[vfid]
    #myclevels=alevels[vfid]

        
    fig = plt.figure(figsize=(myfigsize[0],myfigsize[1]))

    plt.subplots_adjust(wspace=0.05,hspace=.3,bottom=.15)

    if 'VFID5892' in dirname:
        cbaspect = 10
        plt.subplots_adjust(wspace=0.01,bottom=.2)
    maskdat = fits.getdata(mask)

    if contourFlag:
        # get contours from logmstar image
        hdu = fits.open(massim)
        contour_data = hdu[0].data
        contour_header = hdu[0].header
        contour_WCS = WCS(contour_header)
        hdu.close()
        mcontour_data = np.ma.array(contour_data,mask=maskdat)        
        mcontour_data = convolution.convolve_fft(mcontour_data, convolution.Box2DKernel(20), allow_huge=True, nan_treatment='interpolate') 


    ncol = 4
    nrow = 2
    allax = []
    nplots = [2,3,4,6,7,8]
    #for i, im in enumerate(allim):
    for i,j in enumerate(nplots):
        im = allim[i]
        v1 = vmin[i]
        v2 = vmax[i]
        if i > 2:
            COfilename = CO_file[vfid]
            # need to get the header from Francoise's version
            cdat, cheader = convert_alma3d_alma2d(COfilename)
            imwcs = WCS(cheader)
            COheader = cheader
            if i == 3:
                COmask_filename = CO_mask[vfid]
            
                dat = fits.getdata(COmask_filename)
            elif i == 4:
                COmask_filename = CO_vel[vfid]
                dat = fits.getdata(COmask_filename)
            elif i== 5:
                COmask_filename = CO_width[vfid]
                dat = fits.getdata(COmask_filename)
        
        else:
            hdu = fits.open(im)[0]
            dat = hdu.data
            imwcs = WCS(hdu.header)

        if vr is not None:
            # scale data to convert to flux/kpc^2
            kpc_per_pixel = get_kpc_per_pixel(imwcs,vr)
            if i == 0:
                dat = np.log10(10**dat/(kpc_per_pixel**2))
                v1 = vmin[i]-2.*np.log10(kpc_per_pixel) 
                v2 = vmax[i]-2.*np.log10(kpc_per_pixel)           
                #myclevels = myclevels-2.*np.log10(kpc_per_pixel)
                print("\nmass contour",myclevels-2.*np.log10(kpc_per_pixel))
            elif i in [1,2]:
                dat = dat/(kpc_per_pixel**2)
                v1 = vmin[i]/kpc_per_pixel**2 
                v2 = vmax[i]/kpc_per_pixel**2
                #myclevels = myclevels/kpc_per_pixel**2
            
        if i < 3:
            mdat = np.ma.array(dat,mask=maskdat)            
        else:
            mdat = dat


        
        #if i == 2:
        #    ax2=plt.subplot(1,4,i+2,projection=imwcs,slices=('x','y',1))
        #else:
        ax2=plt.subplot(nrow,ncol,j,projection=imwcs)
        

        #if xmin is None:
        #    mdat = mdat
        #else:
        #    mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            try:
                #v1 = vmin[i]
                #v2 = vmax[i]
                plt.imshow(mdat,vmin=v1,vmax=v2,origin='lower',interpolation='nearest',cmap='gray_r')
            except IndexError:
                plt.imshow(mdat,origin='lower',interpolation='nearest',cmap='gray_r')
        elif i in [3,4,5]:
            #v1 = vmin[i]
            #v2 = vmax[i]
            plt.imshow(mdat,vmin=v1,vmax=v2,origin='lower',cmap=cmap,interpolation='nearest')

        else:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            try:
                #v1 = vmin[i]
                #v2 = vmax[i]
                #print(cmap)
                plt.imshow(mdat,vmin=v1,vmax=v2,cmap=cmap)
            except IndexError:
                plt.imshow(mdat)#,cmap=cmap)
            
                #plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],cmap=cmap)#cmap='viridis'

        
        #plt.colorbar(fraction=mycbfrac,aspect=cbaspect)
        # plot contours from mass

        #############################################################
        # add HI contour to halpha image
        #############################################################
        if i == 1:
            # check if HI moment zero map is available
            try:
                HIfilename = HI_file[vfid]
            except KeyError:
                HIfilename = None
            if HIfilename is not None:
                print("HIfilename = ",HIfilename)
                plot_HI_contours(ax2,HIfilename,color='steelblue')
                if vfid == 'VFID5859':
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                else:
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 

        #############################################################
        # add CO contour to halpha image
        #############################################################
        if i == 2:
            # check if HI moment zero map is available
            COfilename = CO_file[vfid]
            COmask_filename = CO_mask[vfid]
            #if COmask_filename is not None:
            #print("COfilename = ",COfilename)

            plot_CO_contours(ax2,COfilename,color=COcolor,mask=COmask_filename,levels=COlevels)
                #if vfid == 'VFID5859':
                #    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                #else:
                #    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 


                    
        #############################################################
        # add contours from stellar mass image
        #############################################################    
        if contourFlag:
            ax = plt.gca()
            ax.contour(mcontour_data,levels=myclevels, colors='k',linestyles='-',linewidths=1,transform=ax.get_transform(contour_WCS))

                
        plt.title(titles[i],fontsize=14)


        # plot CO beam
        if i >= 3:
            plot_HI_beam(plt.gca(),COfilename,COheader,color='steelblue') 

        if xmin is not None:
            
            plt.axis([xmin,xmax,ymin,ymax])

            # plot the legacy image in panel 1
            xcoords = np.array([xmin,xmax])
            ycoords = np.array([ymin,ymax])
    
            # get ramin,ramax and decmin,decmax from SFR image
            sfrim = dirname+"-sfr-vr.fits"
            header = fits.getheader(sfrim)
            #print(header)
            sfrwcs = WCS(header)
            sky = sfrwcs.pixel_to_world(xcoords,ycoords)
    
            # set limits in ra,dec
            x,y = imwcs.world_to_pixel(sky)
            # convert ramin,ramax and decmin,decmax to (x,y)
            #print(sky)
            plt.axis([x[0],x[1],y[0],y[1]])
            
        if i in [0,1]:
            plt.xticks([],[])
            plt.yticks([],[])

        elif not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        
        allax.append(plt.gca())
        myax = plt.gca()
        #if i in [1,2]:
        #
        if i == 23:
            cbaxes = inset_axes(ax2, width="5%", height="60%", loc=7)        
            cb = plt.colorbar(cax=cbaxes, orientation='vertical')
        else:
            cbaxes = inset_axes(ax2, width="80%", height="5%", loc=8)
            cb = plt.colorbar(cax=cbaxes, orientation='horizontal')            
        cb.set_label(label=cblabels[i],fontsize=12)
            
        lon = ax2.coords[0]
        lat = ax2.coords[1]
        if i == 23:
            plt.sca(myax)
            plt.xlabel("RA (hr)",fontsize=14)
            plt.ylabel("DEC (deg)",fontsize=14)
            
            continue
        else:
            lon.set_ticklabel_visible(False)
            lon.set_ticks_visible(False)            
            lat.set_ticklabel_visible(False)
            lat.set_ticks_visible(False)      
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band

    
    ax1 = plt.subplot(nrow,ncol,1,projection=legwcs)
    plt.imshow(jpeg_data)

    ax5 = plt.subplot(nrow,ncol,5,projection=legwcs)
    plt.imshow(jpeg_data)    

    #ax1.set_xlabel("RA",fontsize=16)
    #ax1.set_ylabel("DEC",fontsize=16)    
    leg_axes = [ax1,ax5]
    
    #ax1.set_xlabel("RA",fontsize=16)
    #ax1.set_ylabel("DEC",fontsize=16)    

    if xmin is not None:
        # plot the legacy image in panel 1
        xcoords = np.array([xmin,xmax])
        ycoords = np.array([ymin,ymax])
    
        # get ramin,ramax and decmin,decmax from SFR image
        sfrim = dirname+"-sfr-vr.fits"
        header = fits.getheader(sfrim)
        #print(header)
        sfrwcs = WCS(header)
        sky = sfrwcs.pixel_to_world(xcoords,ycoords)
    
        # set limits in ra,dec
        x,y = legwcs.world_to_pixel(sky)
        # convert ramin,ramax and decmin,decmax to (x,y)
        #print(sky)

        for ax in leg_axes:
            plt.sca(ax)
            plt.axis([x[0],x[1],y[0],y[1]])
            t = dirname.split('-')
            #plt.text(.05,.02,t[0],fontsize=20,transform=plt.gca().transAxes,horizontalalignment='left',color='white')
            #plt.title(f"{t[0]} Legacy grz",fontsize=18)

            # Manually set x-ticks at specific RA values (e.g., every 10 degrees)
            #ax.coords[0].set_ticks(spacing=40 * u.arcsec)
            #ax.xaxis.set_major_locator(ticker.MaxNLocator(2))  # Change 5 to the number of ticks you want
            lon = ax.coords[0]
            lat = ax.coords[1]

            #lon.set_major_formatter('hh:mm')
            if vfid in ['VFID6018']:
                dx = 0.04
            else:
                dx = .01
            #'VFID5855': ,'VFID5842': ,'VFID5709': ,
            lon.set_ticks(spacing=1/60 * u.deg)
            #lon.set_major_formatter('dd:mm')
            lat.set_major_formatter('dd:mm')

            ###############################################
            # add a physical scale
            ###############################################
            pscale = wcs.utils.proj_plane_pixel_scales(legwcs) # in deg -> arcsec
            pscale_arcsec = pscale[0]*3600
            add_scale(ax,vr=vr,pscale=pscale_arcsec,color='w',barsize=2,fontsize=12,yscale=.95)    

            
    plt.sca(ax1)
    #plt.xlabel("RA (hr)",fontsize=14)
    plt.xlabel("",fontsize=14)    
    plt.ylabel("DEC (deg)",fontsize=14)
    plt.title(f"{t[0]}-{t[1]}",fontsize=14)    

    plt.sca(ax5)
    plt.xlabel("RA (hr)",fontsize=14)
    #plt.xlabel("",fontsize=14)    
    plt.ylabel("DEC (deg)",fontsize=14)
    #plt.title(f"{t[0]}-{t[1]}",fontsize=14)    
    #############################################################
    # add arrow showing direction to group center
    # TODO : figure out how to make this smaller in smaller images
    #############################################################
    print("VFID = ",vfid)

    if (vfid == 'VFID5855') | (vfid == 'VFID5859'):
        print('hi')
        # get coords of galaxy
        vflag = v.main['VFID'] == vfid
        print(np.sum(vflag))
        x = v.main['RA'][vflag][0]
        y = v.main['DEC'][vflag][0]      
        dRA = v.mw_RAcenter - x
        dDEC = v.mw_DECcenter - y        
        norm = np.sqrt(dRA**2+dDEC**2)

        if  (vfid == 'VFID5855'):
            scale = 0.02
            mylw = 0.05
            dx = dRA/norm*scale
            dy = dDEC/norm*scale

            #plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),lw=mylw,alpha=.7)
        elif (vfid == 'VFID5859'):
            scale = 0.0025
            mylw=.5
            dx = dRA/norm*scale
            dy = dDEC/norm*scale
            width=.001
            plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),linewidth=mylw,alpha=.7,head_width=2.5*width,head_length=2.5*width)
        print(x,y)
        #plt.plot(x*u.deg,y*u.deg,'bo',markersize=20,color='c',transform=ax1.get_transform('world'))
        

    #print(x[0],x[1],y[0],y[1])
    #############################################################
    # add HI contour to legacy image
    #############################################################    
    # check if HI moment zero map is available
    try:
        HIfilename = HI_file[vfid]
    except KeyError:
        HIfilename = None
    if HIfilename is not None:
        print("HIfilename = ",HIfilename)
        plot_HI_contours(ax1,HIfilename,color='lightsteelblue')
        plot_HI_beam(ax1,HIfilename,fits.getheader(legacyr),color='steelblue') 
        
    #############################################################
    # add stellar mass from magphys to the legacy image
    #############################################################    
    if logMstar is not None:
        print("adding logMstar = ",logMstar)
        #plt.text(0.05,0.05,logMstar,fontsize=16,color='white',transform=plt.gca().transAxes,horizontalalignment='left')
    
    
    #ax1.update({'xlabel': 'RA', 'ylabel': 'DEC'})


    # check if HI moment zero map is available
    COfilename = CO_file[vfid]
    COmask_filename = CO_mask[vfid]
    plot_CO_contours(ax5,COfilename,color=COcolor,mask=COmask_filename,levels=COlevels)

    plot_HI_beam(ax5,COfilename,COheader,color='steelblue')

    
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-profiles-COall.png',dpi=150,bbox_inches="tight")
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-profiles-COall.pdf',dpi=150,bbox_inches="tight")    


    os.chdir(cwd)


    return ax1


def plot_mstar_sfr_COimages(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,\
                            logMstar=None,cmap='magma_r',markGroupCenter=False,COcolor='white'):
    """
    same plot as mstar_sfr, but swap out CO image radial profiles in the 4th panel
    add other moments as well to make 1x7 panel plot

    """
    #%matplotlib inline
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name

    vfid = dirname.split('-')[0]    

    cwd = os.getcwd()
    os.chdir(dirname)
    massim = dirname+"-logmstar-vr.fits"
    sfrim = dirname+"-sfr-vr.fits"

    
    ssfrim = dirname+"-ssfr.fits"
    mask = dirname+'-R-mask.fits'
    #titles = [r'$\log_{10}(M_\star/M_\odot)$',r'$H\alpha \ SFR$','log sSFR']
    titles = [r'$\log_{10}(\Sigma_\star)$',r'$\Sigma_{SFR} \ w/HI$',r'$CO$']
    vmin = [2,0e-5,0]
    vmax = [6,.6e-5,.1]
    allim = [massim,sfrim,ssfrim]
    #allim = [massim,sfrim]
    #allim = [massim,sfrim,sfrim] # show halpha with HI and CO
    allim = [massim,sfrim,CO_file[vfid]] # show halpha with HI and CO    
    cblabels = [r'$\rm \log_{10}(M_\star/M_\odot)/pixel$',r'$\rm SFR(M_\star/yr)/pixel$',r'$\rm CO/pixel$']
    

    print("VFID = ",vfid)
    if 'INT' in dirname:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
        except KeyError:
            # use the full image
            # read in massim
            data = fits.getdata(massim)
            ymax,xmax = data.shape
            xmin=1
            ymin=1

    else:
        try:
            xmin,xmax,ymin,ymax = zoom_coords_HDI[vfid]
        except KeyError:
            # use the full image
            # read in massim
            data = fits.getdata(massim)
            ymax,xmax = data.shape
            xmin=1
            ymin=1
    


    myfigsize=afigsize[vfid]
    mycbfrac=acbfrac[vfid]
    myclevels=alevels[vfid]

        
    fig = plt.figure(figsize=(myfigsize[0],myfigsize[1]))

    plt.subplots_adjust(wspace=0.01,bottom=.15)

    if 'VFID5892' in dirname:
        cbaspect = 10
        plt.subplots_adjust(wspace=0.01,bottom=.2)
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
        if i == 2:
            COfilename = CO_file[vfid]
            COmask_filename = CO_mask[vfid]
            dat, cheader = convert_alma3d_alma2d(COfilename)
            if COmask_filename is not None:
                maskdat = fits.getdata(COmask_filename)
                maskdat = maskdat == 1
                maskdat = ~maskdat
            else:
                maskdat = np.zeros_like(dat,'i')

            #maskdat = np.zeros_like(dat)
            imwcs = WCS(cheader)
            #plt.figure()
            #plt.imshow(dat)
            
            pass
        else:
            hdu = fits.open(im)[0]
            dat = fits.getdata(im)
            imwcs = WCS(fits.getheader(im))

            dat = hdu.data
            imwcs = WCS(hdu.header)
        #if i == 2:
        #    ax2=plt.subplot(1,4,i+2,projection=imwcs,slices=('x','y',1))
        #else:
        ax2=plt.subplot(1,4,i+2,projection=imwcs)
        
                        


        mdat = np.ma.array(dat,mask=maskdat)
        
        #if xmin is None:
        #    mdat = mdat
        #else:
        #    mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            plt.imshow(mdat,vmin=0,vmax=.2,interpolation='nearest')#,vmin=vmin[i],vmax=vmax[i],origin='lower',interpolation='nearest')

            continue
        else:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            plt.imshow(mdat,vmin=vmin[i],vmax=vmax[i],cmap=cmap)#cmap='viridis'
        
        lon = ax2.coords[0]
        lat = ax2.coords[1]
        #lon.set_ticklabel_visible(False)
        #lon.set_ticks_visible(False)            
        #lat.set_ticklabel_visible(False)
        #lat.set_ticks_visible(False)            
        
        #plt.colorbar(fraction=mycbfrac,aspect=cbaspect)
        # plot contours from mass

        #############################################################
        # add HI contour to halpha image
        #############################################################
        if i == 1:
            # check if HI moment zero map is available
            HIfilename = HI_file[vfid]
            if HIfilename is not None:
                print("HIfilename = ",HIfilename)
                plot_HI_contours(ax2,HIfilename,color='steelblue')
                if vfid == 'VFID5859':
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                else:
                    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 

        #############################################################
        # add CO contour to halpha image
        #############################################################
        #if i == 2:
            # check if HI moment zero map is available
        #    COfilename = CO_file[vfid]
        #    COmask_filename = CO_mask[vfid]
        #    if COfilename is not None:
        #        print("COfilename = ",COfilename)

                #plot_CO_image(ax2,COfilename,color=COcolor,mask=COmask_filename)
                #if vfid == 'VFID5859':
                #    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue',expandBox=True)
                #else:
                #    plot_HI_beam(ax2,HIfilename,hdu.header,color='steelblue') 


                    
        #############################################################
        # add contours from stellar mass image
        #############################################################    
        if contourFlag and (i < 2):
            ax = plt.gca()
            ax.contour(mcontour_data,levels=myclevels, colors='k',linestyles='-',linewidths=1,transform=ax.get_transform(contour_WCS))

                
        plt.title(titles[i],fontsize=18)



        if (xmin is not None) and (i < 2):
            plt.axis([xmin,xmax,ymin,ymax])
        if i in [0,1]:
            plt.xticks([],[])
            plt.yticks([],[])

        elif not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        
        allax.append(plt.gca())

        #if i in [1,2]:
        cbaxes = inset_axes(ax2, width="80%", height="3%", loc=8) 
        cb = plt.colorbar(cax=cbaxes, orientation='horizontal')
        cb.set_label(label=cblabels[i],fontsize=12)

    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band

    
    ax1 = plt.subplot(1,4,1,projection=legwcs)
    plt.imshow(jpeg_data)
    leg_axes = [ax1]
    if xmin is not None:
        # plot the legacy image in panel 1
        xcoords = np.array([xmin,xmax])
        ycoords = np.array([ymin,ymax])
    
        # get ramin,ramax and decmin,decmax from SFR image
        sfrim = dirname+"-sfr-vr.fits"
        header = fits.getheader(sfrim)
        #print(header)
        sfrwcs = WCS(header)
        sky = sfrwcs.pixel_to_world(xcoords,ycoords)
    
        # set limits in ra,dec
        x,y = legwcs.world_to_pixel(sky)
        # convert ramin,ramax and decmin,decmax to (x,y)
        #print(sky)
        for ax in leg_axes:
            plt.sca(ax)
            plt.axis([x[0],x[1],y[0],y[1]])
            t = dirname.split('-')
  

            plt.title(f"{t[0]} Legacy grz",fontsize=18)
            ax.xaxis.set_major_locator(ticker.MaxNLocator(2))  # Change 5 to the number of ticks you want
    #############################################################
    # add arrow showing direction to group center
    # TODO : figure out how to make this smaller in smaller images
    #############################################################
    print("VFID = ",vfid)

    if (vfid == 'VFID5855') | (vfid == 'VFID5859'):
        print('hi')
        # get coords of galaxy
        vflag = v.main['VFID'] == vfid
        print(np.sum(vflag))
        x = v.main['RA'][vflag][0]
        y = v.main['DEC'][vflag][0]      
        dRA = v.mw_RAcenter - x
        dDEC = v.mw_DECcenter - y        
        norm = np.sqrt(dRA**2+dDEC**2)

        if  (vfid == 'VFID5855'):
            scale = 0.02
            mylw = 0.05
            dx = dRA/norm*scale
            dy = dDEC/norm*scale

            plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),lw=mylw,alpha=.7)
        elif (vfid == 'VFID5859'):
            scale = 0.0025
            mylw=.5
            dx = dRA/norm*scale
            dy = dDEC/norm*scale
            width=.001
            plt.arrow(x,y,dx,dy,color='c',transform=ax1.get_transform('world'),linewidth=mylw,alpha=.7,head_width=2.5*width,head_length=2.5*width)
        print(x,y)
        #plt.plot(x*u.deg,y*u.deg,'bo',markersize=20,color='c',transform=ax1.get_transform('world'))
        

    #print(x[0],x[1],y[0],y[1])
    #############################################################
    # add HI contour to legacy image
    #############################################################    
    # check if HI moment zero map is available
    HIfilename = HI_file[vfid]
    if HIfilename is not None:
        print("HIfilename = ",HIfilename)
        plot_HI_contours(plt.gca(),HIfilename,color='lightsteelblue')
        plot_HI_beam(plt.gca(),HIfilename,fits.getheader(legacyr),color='steelblue') 
        
    #############################################################
    # add stellar mass from magphys to the legacy image
    #############################################################    
    if logMstar is not None:
        print("adding logMstar = ",logMstar)
        plt.text(0.05,0.05,logMstar,fontsize=16,color='white',transform=plt.gca().transAxes,horizontalalignment='left')
    
    
    ax1.update({'xlabel': 'RA', 'ylabel': 'DEC'})


    plt.sca(ax1)
    plt.xlabel("RA (hr)",fontsize=14)
    plt.ylabel("DEC (deg)",fontsize=14)
    plt.savefig(os.path.join(plotdir,dirname)+'-mstar-sfr-profiles-CO.png',dpi=150,bbox_inches="tight")


    os.chdir(cwd)


    return ax1


def fit_profiles(rp,hp,weights=None,rmax=None,fixN=False,labels = ['logMstar','logSFR'],log1Flag=False):
    """
    Fit a sersic model to a light profile

    PARAMS: 
    * rp = r-band photometry file
    * hp = halpha photometry file
    * rmax = max radius to use for fitting

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
    '''
    rmax in arcsec
    '''
    
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
    radii50 = []
    radii90 = []
    for i,f in enumerate(fits):
        plt.subplot(2,1,i+1)
        if rmax is not None:
            plt.xlim(0,rmax)
            flag = phots[i]['sma_arcsec'] < rmax
        else:
            flag = np.ones(len(phots),'bool')
        #plt.gca().set_yscale('log')

        plt.plot(phots[i]['sma_arcsec'][flag],phots[i]['flux'][flag],'bo',c=mycolors[0],label=labels[i])
        
        plt.axhline(y=0.5*np.max(phots[i]['flux'][flag]),ls=':',color='k')
        plt.axvline(x=fits[i].r_eff.value,ls='--',color='k')

    
        # plot R50, where the enclosed flux reaches the 0.5*total
        if rmax is not None:

            r50flag = phots[i]['flux'] > 0.5*np.max(phots[i]['flux'][flag])
            r50 = phots[i]['sma_arcsec'][r50flag][0]
            radii50.append(r50)
            
            r90flag = phots[i]['flux'] > 0.9*np.max(phots[i]['flux'][flag])
            r90 = phots[i]['sma_arcsec'][r90flag][0]
            radii90.append(r90)
        else:
            r50flag = phots[i]['flux'] > 0.5*np.max(phots[i]['flux'])
            r50 = phots[i]['sma_arcsec'][r50flag][0]
            radii50.append(r50)

            r90flag = phots[i]['flux'] > 0.9*np.max(phots[i]['flux'])
            r90 = phots[i]['sma_arcsec'][r90flag][0]
            radii90.append(r90)
            

        plt.axvline(x=r50,ls=':',color='k')                                     
        plt.text(0.05,0.8,f"{names[i]}: Re = {fits[i].r_eff.value:.1f}, n={f.n.value:.1f}, R50={r50:.1f}",transform=plt.gca().transAxes)
                                 
        if i == 1:

            plt.xlabel("SMA (arcsec)",fontsize=16)
            plt.text(-.15,1,"Enclosed Flux",fontsize=16,transform=plt.gca().transAxes,rotation=90,horizontalalignment='center')
    return radii50,radii90
    
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
    plt.savefig(vfid+'-mstar-sfr-sersic-profiles.png',dpi=150,bbox_inches="tight")
        
    radii50,radii90 = plot_cog(rp,hp,mfit,sfit,rmax=rmax,labels=['logMstar','logSFR'])
    plt.savefig(vfid+'-mstar-sfr-cog.png',dpi=150,bbox_inches="tight")    
    # use this to run on R and CS Halpha
    rphot = dirname+'-R_phot.fits'
    haphot = dirname+'-CS-gr_phot.fits'
    try:
        rp = Table.read(rphot)
        hp = Table.read(haphot)
    except FileNotFoundError:
        rphot = dirname+'-R-phot.fits'
        haphot = dirname+'-CS-phot.fits'
        rp = Table.read(rphot)
        hp = Table.read(haphot)
    rfit,hfit = fit_profiles(rp,hp,rmax=rmax,labels=['r','halpha'])
    plt.savefig(vfid+'-r-halpha-profiles.png',dpi=150,bbox_inches="tight")    
    return mfit,sfit,rfit,hfit,radii50,radii90


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
    vmins = [0.002,0,0]
    vmaxs = [.03,.00002,14000]    
    print(nuvim,w3im)
    
    if xmin is not None:
        # get ra  and dec corresponding to xmin,ymin xmax,ymax from sfr image
        # plot the legacy image in panel 1
        xcoords = np.array([xmin,xmax])
        ycoords = np.array([ymin,ymax])

        header = fits.getheader(sfrim)
        sfrwcs = WCS(header)
        sky = sfrwcs.pixel_to_world(xcoords,ycoords)
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
        #display_image(mdat,percent=99.9,cmap='magma_r')#,vmin=vmin[i],vmax=vmax[i])
        if i == 1:
            plt.imshow(mdat,cmap='magma_r',origin='lower',vmin=vmins[i],vmax=vmaxs[i])
        else:
            plt.imshow(mdat,vmin=vmins[i],vmax=vmaxs[i],cmap='magma_r',origin='lower')
        if not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        plt.title(titles[i],fontsize=20)
        plt.colorbar(fraction=cbfrac,aspect=cbaspect)
        # plot contours from mass
        allax.append(plt.gca())
    # plot the legacy image in panel 1
    xcoords = np.array([xmin,xmax])
    ycoords = np.array([ymin,ymax])
    
    # get ramin,ramax and decmin,decmax from SFR image
    sfrim = dirname+"-sfr-vr.fits"
    header = fits.getheader(sfrim)
    sfrwcs = WCS(header)
    sky = sfrwcs.pixel_to_world(xcoords,ycoords)
    
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
    #plt.title("Legacy grz",fontsize=20)
    t = dirname.split('-')
    plt.title(t[0],fontsize=20)
    plt.xlabel('RA (deg)',fontsize=20)
    plt.ylabel('DEC (deg)',fontsize=20)
    
    # read in the
    plt.savefig(os.path.join(plotdir,dirname)+'-sfr-indicators.png',dpi=150,bbox_inches="tight")
    plt.show()

    os.chdir(cwd)


def plot_sfr_indicators_nohalpha(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],cbfrac=.08,cbaspect=20,clevels=[4,5],contourFlag=False,logMstar=None):
    #%matplotlib inline
    #os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name

    print()
    cwd = os.getcwd()
    os.chdir(dirname)
    massim = dirname+"-logmstar-vr.fits"
    #sfrim = dirname+"-sfr-vr.fits"
    #ssfrim = dirname+"-ssfr.fits"
    #mask = dirname+'-R-mask.fits'
    vfid = dirname.split('-')[0]
    print(dirname+'/galex/*nuv*.fits')
    #print("sfrim = ",sfrim)
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
        t = glob.glob('unwise/*w3-img-m.fits')
        w3im = t[0]
        
        
    titles = ['NUV',r'$unWISE \ 12\mu m$']
    #vmin = [2,0,-11.5]
    #vmax = [6,1.e-6,-9]
    allim = [nuvim,w3im]
    vmins = [0.002,0]
    vmaxs = [.03,14000]    
    #print(nuvim,w3im)
    
        
    plt.figure(figsize=(figsize[0],figsize[1]))

    plt.subplots_adjust(wspace=0.05)
    allax = []
    for i, im in enumerate(allim):
        plt.subplot(1,3,i+2)
        dat = fits.getdata(im)
        
            
        if i == 1:
            plt.imshow(dat,cmap='magma_r',origin='lower',vmin=vmins[i],vmax=vmaxs[i])
        else:
            plt.imshow(dat,vmin=vmins[i],vmax=vmaxs[i],cmap='magma_r',origin='lower')
        if not xticks: 
            plt.xticks([],[])
            plt.yticks([],[])
        plt.title(titles[i],fontsize=20)
        #plt.colorbar(fraction=cbfrac,aspect=cbaspect)
        if i == 1:
            t = dirname.split('-')
            if len(t) == 1:
                s = t[0]
            elif len(t) == 2:
                s = t[0]+' '+t[1]
            #plt.xlabel(s,fontsize=20)
        allax.append(plt.gca())
    # plot the legacy image in panel 1
    
    # read in header from legacy r-band image
    legacyr = glob.glob("legacy/*r.fits")[0]
    #print(legacyr)
    legacy_jpg = legacyr.replace('-r.fits','.jpg')
    jpeg_data = Image.open(legacy_jpg)
    legwcs = WCS(fits.getheader(legacyr))
    # plot jpg as projection of legacy r-band
    plt.subplot(1,3,1,projection=legwcs)
    ax1 = plt.gca()
    plt.imshow(jpeg_data)
    # set limits in ra,dec
    #x,y = legwcs.world_to_pixel(sky)
    # convert ramin,ramax and decmin,decmax to (x,y)
    #plt.title("Legacy grz",fontsize=20)
    plt.title(vfid,fontsize=20)
    # read in the
    plt.xlabel("RA (hr)",fontsize=18)
    plt.ylabel("DEC (deg)",fontsize=18)

    # add stellar mass
    #############################################################
    # add stellar mass from magphys to the legacy image
    #############################################################    
    if logMstar is not None:
        print("adding logMstar = ",logMstar)
        plt.text(0.05,0.05,f"{logMstar:.1f}",fontsize=18,color='white',transform=plt.gca().transAxes,horizontalalignment='left')

    #############################################################
    # add stellar mass contours
    #############################################################    
    if contourFlag:
        # get contours from logmstar image
        hdu = fits.open(massim)
        contour_data = hdu[0].data
        contour_header = hdu[0].header
        contour_WCS = WCS(contour_header)
        hdu.close()
        mcontour_data = np.ma.array(contour_data,mask=maskdat)
        
        ax = plt.gca()
        myclevels = [4]
        ax.contour(mcontour_data,levels=myclevels, colors='k',linestyles='-',linewidths=1,transform=ax.get_transform(contour_WCS))
    
    ###############################################
    # add a physical scale
    ###############################################
    pscale = wcs.utils.proj_plane_pixel_scales(legwcs) # in deg -> arcsec
    pscale_arcsec = pscale[0]*3600
    add_scale(ax1,vr=v.mw_vrcenter,pscale=pscale_arcsec,color='w',barsize=2)    
    
    
    plt.savefig(os.path.join(plotdir,dirname)+'-sfr-indicators-nohalpha.png',dpi=150,bbox_inches="tight")
    plt.savefig(os.path.join(plotdir,dirname)+'-sfr-indicators-nohalpha.pdf',dpi=150,bbox_inches="tight")    
    plt.show()

    os.chdir(cwd)
    
def Jaffe15_phase_space_region(c=2.8,sigma=1.,M200=1e+14,R200=1.55):
    """                                                                                                                                               
    PURPOSE:                                                                                                                                          
    macro that plots the virialized region in the phase space diagram following Jaffe et al. 2015, MNRAS 448.1715J                                     
    
    INPUT:                                                                                                                                            
    c : r200/rs ,  concentration parameter                                                                                                            
    sigma : velocity dispersion (km/s)                                                                                                                
    M200 : M200 mass (Msun)                                                                                                                          
    R200 : R200 (Mpc)                                                                                                                                
    """
    G = 6.67260e-8 #(cm3 g-1 s-2) Gravitational constant                                                                                           
    Msun = 1.99e+33 #(g) solar mass                                                                                                                  
    Mpc = 3.086e+24 #(cm) Mpc                                                                                                                        
    km = 1e+5 #(cm)                                                                                                                                   
    ##########                                                                                                                                        
    M200 *= Msun #M200 in g                                                                                                                           
    R200 *= Mpc #R200 in cm                                                                                                                          
    sigma *= km #sigma in cm/s                                                                                                                        
    ##################                                                                                                                               
    r_to_r200 = np.arange(0.,3.,0.01) #r/r200 ratio                                                                                                   
    vesc = 0.*r_to_r200-1000.
    ###                                                                                                                                              
    s = r_to_r200                                                                                                                          
    gc = 1./(np.log(1.+c)-c/(1.+c))
    K = gc*(np.log(1+c*s)/s-np.log(1+c))+1.
    ######################                                                                                                                            
    cond = (r_to_r200<1.)
    vesc[cond] = np.sqrt(2*G*M200*K[cond]/(3.*R200))
    vesc[~cond] = np.sqrt(2*G*M200/(3.*(R200*s[~cond])))
    return r_to_r200,vesc/sigma

def plot_gianluca_caustics(sepmax=15):
    """ code from Gianluca to plot caustics for Virgo """
    ################
    #I plot an hyperbolic caustic derived empirically from Fig 4 (Kim_14)
    #y= a/x+b, two points in the fig are A = (1.7,3000); B = (20,1000)
    a = 3968.0
    b = 1016. #cosmic velocity of Virgo from Mould+2000 (see NED database: "V (Virgo + GA only) : Quantities derived from redshift)    
    x = np.arange(0.0001,sepmax,0.01)
    y_up = a/x+b
    y_down = -1*a/x+b
    ################
    """
    phase space model (Jaffe+15)
    """
    G = 6.67260e-8 #(cm3 g-1 s-2) Gravitational constant                                                                                              
    Msun = 1.99e+33 #(g) solar mass                                                                                                                   
    Mpc = 3.086e+24 #(cm) Mpc                                                                                                                         
    km = 1e+5 #(cm)                                                                                                                                   
    H0 = 74.
    #############
    r200 = 1.55*70./H0 #Mpc (McLaughlin 1999)
    concentration = 2.8 #concentration parameter (McLaughlin 1999)
    rhoc = 3./(8.*np.pi*G)*((H0/Mpc*1e+5)**2) #g/cm3 - critical density
    M200 = 4./3*np.pi*((r200*Mpc)**3)*200.*rhoc/Msun #Mo - M200 mass
    phase_space = Jaffe15_phase_space_region(c=concentration,sigma=1.,M200=M200,R200=r200)
    separation = phase_space[0]*r200*H0/b*180./np.pi
    #separation = phase_space[0]*r200*H0/b*180./np.pi
    vel_up = phase_space[1]+b
    vel_down = b-phase_space[1]
    
    plt.plot(separation,vel_up,'k-')
    plt.plot(separation,vel_down,'k-')    

def plot_sky_positions():
    plt.figure(figsize=(8,6))

    x = v.main['RA']
    y = v.main['DEC']
    c = v.dist3dVirgo
    plt.scatter(x,y,c='k',marker='.',alpha=.1,vmin=1,vmax=10,s=8)
    plt.scatter(x[v.groupMembs],y[v.groupMembs],c=c[v.groupMembs],alpha=1,vmin=1,vmax=10,label='NGC5364 Group')
    plt.scatter(x[v.env['cluster_member']],y[v.env['cluster_member']],alpha=.4,vmin=1,vmax=10,s=50,marker='*',label='Virgo Cluster Members')
    #plt.colorbar(label="3D Distance (Mpc)")
    plt.axis([170,214,-3,28])
    plt.gca().invert_xaxis()
    plt.ylabel("RA (deg)",fontsize=16)
    plt.xlabel("DEC (deg)",fontsize=16)



    spinedir = '/home/rfinn/research/Virgo/tables-north/v2/spines/'
    spine  = Table.read(spinedir+'filament_spine_VirgoIII.fits')
    plt.plot(spine['ra'],spine['dec'],'k--',label='Virgo III Filament Spine',alpha=.6)

    plt.subplots_adjust(right=.9)
    plt.legend(loc='upper left')
    mainax = plt.gca()
    x1,x2 = plt.xlim()
    y1,y2 = plt.ylim()
    # add projected separation on other axes

    # virgo cluster
    RAvirgo,DECvirgo = 187.697083, 12.336944


    from astropy.cosmology import WMAP9 as cosmo

    H0 = cosmo.H0.value

    DA = cosmo.angular_diameter_distance(1100/3.e5)
    DAperdeg = DA*np.pi/180


    ax2 =mainax.twiny()

    ax2.set_xlim((x1-RAvirgo)*DAperdeg.value,(x2-RAvirgo)*DAperdeg.value)
    ax2.set_xlabel("Angular Separation (Mpc)",fontsize=16)
    plt.grid()

    ax2 =mainax.twinx()

    ax2.set_ylim((y1-DECvirgo)*DAperdeg.value,(y2-DECvirgo)*DAperdeg.value)
    ax2.set_ylabel("Angular Separation (Mpc)",fontsize=16,rotation=270,labelpad=20)
    plt.grid()
    plt.savefig(plotdir+'/NGC5364_Virgo_sky_positions.png',dpi=150,bbox_inches="tight")

def plot_phase_space(sepmax=15):
    # like Fig 3 in Castignani+2022b

    # virgo cluster
    RAvirgo,DECvirgo = 187.697083, 12.336944
    vr = 1150
    sigma=600
    # phase space

    dv = (v.main['vr'])

    dr = np.sqrt((RAvirgo-v.main['RA'])**2 + (DECvirgo - v.main['DEC'])**2)
    
    plt.figure(figsize=(8,6))

    x = dr
    y = dv
    c = v.dist3dVirgo

    plt.scatter(x,y,c='k',alpha=.1,vmin=1,vmax=10,s=6)

    #plt.scatter(x[v.env['cluster_member']],y[v.env['cluster_member']],c=c[v.env['cluster_member']],alpha=.4,vmin=1,vmax=10,s=50,marker='*',label='Cluster member')
    plt.scatter(x[v.env['cluster_member']],y[v.env['cluster_member']],alpha=.4,vmin=1,vmax=10,s=50,marker='*',label='Virgo Cluster members')
    plt.scatter(x[v.groupMembs],y[v.groupMembs],c=c[v.groupMembs],alpha=1,vmin=1,vmax=10,s=100,label='NGC 5364 Group')
    #plt.scatter(x[v.groupMembs],y[v.groupMembs],c='c',alpha=1,vmin=1,vmax=10,s=100,label='NGC 5364 Group')

    plot_gianluca_caustics(sepmax=sepmax)
    
    plt.colorbar(label="3D Distance (Mpc)")
    plt.xlabel("Angular Separation (deg)",fontsize=16)
    plt.ylabel(r"$\rm v_r \ (km/s)$",fontsize=16)
    plt.axis([-.2,26,-1000,3000])
    plt.legend()
    #plt.axhline(y=1100,ls='--',color='k',alpha=.5)



    from astropy.cosmology import WMAP9 as cosmo

    H0 = cosmo.H0.value

    DA = cosmo.angular_diameter_distance(1100/3.e5)
    DAperdeg = DA*np.pi/180
    x1,x2 = plt.xlim()
    print(x1,x2)

    ax2 = plt.gca().twiny()

    ax2.set_xlim(x1*DAperdeg.value,x2*DAperdeg.value)
    ax2.set_xlabel("Angular Separation (Mpc)",fontsize=16)
    
    plt.savefig(plotdir+'/NGC5364_Virgo_phasespace.png',dpi=150,bbox_inches="tight")


def plot_phase_space_normalized(sepmax=15):
    # like Fig 3 in Castignani+2022b

    # virgo cluster
    RAvirgo,DECvirgo = 187.697083, 12.336944
    vr = 1150
    sigma=600
    # phase space

    dv = (v.main['vr'] - v.mw_vrcenter)/v.bw_scale

    dr = np.sqrt((v.mw_RAcenter-v.main['RA'])**2 + (v.mw_DECcenter - v.main['DEC'])**2)/v.rvirial_deg
    
    plt.figure(figsize=(8,6))

    x = dr
    y = dv
    c = v.dist3dVirgo
    c = v.combinedMass
    plt.scatter(x[v.groupMembs],y[v.groupMembs],c='k',alpha=.3,vmin=1,vmax=10,s=20)

    #plt.scatter(x[v.env['cluster_member']],y[v.env['cluster_member']],c=c[v.env['cluster_member']],alpha=.4,vmin=1,vmax=10,s=50,marker='*',label='Cluster member')
    #plt.scatter(x[v.env['cluster_member']],y[v.env['cluster_member']],alpha=.4,vmin=1,vmax=10,s=50,marker='*',label='Virgo Cluster members')
    #flag = v.groupMembs & ~v.groupMembsHA
    #plt.scatter(x[v.groupMembsHA],y[v.groupMembsHA],c=c[v.groupMembsHA],marker='s',alpha=1,vmin=8,vmax=10,s=100,label=r'$\rm On H\alpha \ FOV$')    
    #plt.scatter(x[flag],y[flag],c=c[flag],marker='o',alpha=1,vmin=8,vmax=10,s=50,label='Other Members')

    #plt.scatter(x[v.groupMembs],y[v.groupMembs],c='c',alpha=1,vmin=1,vmax=10,s=100,label='NGC 5364 Group')


    vfids = ['VFID5842','VFID5892','VFID5889','VFID5855','VFID5859']#,'VFID5851','VFID5844','VFID5879']
    for i,vf in enumerate(vfids):
        flag = v.main['VFID'] == vf
        plt.plot(x[flag],y[flag],'bs',c=mycolors[i],label=vf,markersize=12)#marker='s',alpha=1,vmin=8,vmax=10,s=100,label=r'$\rm On H\alpha \ FOV$')


    vfids = ['VFID5851','VFID5844','VFID5879']
    for i,vf in enumerate(vfids):
        flag = v.main['VFID'] == vf
        plt.plot(x[flag],y[flag],'bv',c='k',alpha=.6,label=vf,markersize=12)#marker='s',alpha=1,vmin=8,vmax=10,s=100,label=r'$\rm On H\alpha \ FOV$')
        
    #plot_gianluca_caustics(sepmax=sepmax)
    ids =  np.arange(len(v.main))[v.groupMembsHA]
    #for i in ids:
    #    plt.text(dr[i]-.01,dv[i]+.1,v.main['VFID'][i])    
    #plt.colorbar(label="$\log (M/M_\odot)$")
    plt.xlabel(r"$\rm \Delta r/R_{virial}$",fontsize=16)
    plt.ylabel(r"$\rm \Delta v/\sigma_{bw}$",fontsize=16)
    plt.axis([0,1,-1.4,2.1])
    plt.legend()#loc='center right')
    plt.axhline(y=0,ls='--',color='k',alpha=.5)



    from astropy.cosmology import WMAP9 as cosmo

    H0 = cosmo.H0.value

    DA = cosmo.angular_diameter_distance(1100/3.e5)
    DAperdeg = DA*np.pi/180
    x1,x2 = plt.xlim()
    print(x1,x2)

    #ax2 = plt.gca().twiny()

    #ax2.set_xlim(x1*DAperdeg.value,x2*DAperdeg.value)
    #ax2.set_xlabel("Angular Separation (Mpc)",fontsize=16)
    
    plt.savefig(plotdir+'/NGC5364_Virgo_phasespace_normalized.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/NGC5364_Virgo_phasespace_normalized.pdf')    


def plot_sky_positions_with_filaments(multicolor=True,plotlegend=True):

    plt.figure(figsize=(11,5))
    plt.subplots_adjust(right=.95,top=.95)
    plot_spines(multicolor=multicolor,colorone='VirgoIII',color='m',legend=plotlegend)
    plt.plot(v.main['RA'],v.main['DEC'],'k.',alpha=.1,label='Virgo Filament Catalog')
    # plot CO sources
    flag = v.main['HAobsflag'] #& ~v.main['COflag']
    #plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'cs',alpha=.5,mfc='None',label=r'$H\alpha$',markersize=8)
    flag = v.groupMembs
    plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'co',alpha=.8,label=r'$\rm NGC5364~Group$',markersize=5)    


    

    #flag = v.main['COflag'] & ~v.main['HAobsflag']
    #plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'bo',color='purple',markersize=8,alpha=.7,label=r'$CO, No \ H\alpha$',mec='0.5')
    #flag = v.main['COflag'] & v.main['HAobsflag']
    #plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'co',markersize=5,alpha=.7,label=r'$CO + H\alpha$',mec='0.5')
    plt.axis([118,262,-1,65])
    plt.legend(loc='lower right')#bbox_to_anchor=(1.01,1),loc='upper left')



    plt.gca().invert_xaxis()
    plt.xlabel('RA (deg)',fontsize=20)
    plt.ylabel('DEC (deg)',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.legend(loc='upper right')
    #plt.title('Filamentary Structures Surrounding the Virgo Cluster',fontsize=18)
    plt.savefig(plotdir+'/ngc5364_positions_filaments.png',dpi=150,bbox_inches="tight")
    plt.savefig(plotdir+'/ngc5364_positions_filaments.pdf')
class grouptables(vtables):
    def get_group_members(self):
        self.groupMembs = self.kt['PGC1'] == 49547
        print(f"number of Kourkchi group members = {np.sum(self.groupMembs)}")
        self.groupMembsHA = (self.kt['PGC1'] == 49547) & (self.main['HAobsflag'])
    def get_mass_weighted_center(self):
        """get RA and DEC of mass weighted center  """

        # create a combined mass that uses the best fit magphys mass for those with bad med values
        badMassFlag = v.magphys['logMstar_med'] < 2
        combinedMass = v.magphys['logMstar_med'] * ~badMassFlag + v.magphys['logMstar_best'] * badMassFlag
        self.combinedMass = combinedMass
        # calculate the mass-weighted center
        self.mw_RAcenter = np.sum(v.main['RA'][self.groupMembs]*10.**combinedMass[self.groupMembs])/np.sum(10.**combinedMass[self.groupMembs])

        self.mw_DECcenter = np.sum(v.main['DEC'][self.groupMembs]*10.**combinedMass[self.groupMembs])/np.sum(10.**combinedMass[self.groupMembs])

        self.mw_vrcenter = np.sum(v.main['vr'][self.groupMembs]*10.**combinedMass[self.groupMembs])/np.sum(10.**combinedMass[self.groupMembs])
        self.mw_Vcosmic_center = np.sum(v.env['Vcosmic'][self.groupMembs]*10.**combinedMass[self.groupMembs])/np.sum(10.**combinedMass[self.groupMembs])

    def get_dr_rvirial(self):
        self.rvirial = 222. # kpc

        # projected distance from group center
        dRA = np.abs(v.main['RA'][self.groupMembsHA] - self.mw_RAcenter)
        dDEC = np.abs(v.main['DEC'][self.groupMembsHA] - self.mw_DECcenter)
        dtheta_deg = np.sqrt(dRA**2 + dDEC**2)

        # convert to kpc
        from astropy.cosmology import WMAP9 as cosmo
        H0 = cosmo.H0.value

        redshift = self.mw_Vcosmic_center/3.e5
        
        DA = cosmo.angular_diameter_distance(redshift) # Mpc per radian
        self.rvirial_deg = np.degrees(self.rvirial/1000/DA.value)
        
        dr_kpc = np.radians(np.mean(dtheta_deg))*DA*1000
        print(f"average separation of Halpha FOV galaxies from center is {dr_kpc.value/self.rvirial:.3f} Rvirial")
        
    def write_latex_table(self,flag=None,outfile=None):
        self.get_distance_Virgo()
        if flag is None:
            flag = self.groupMembs & self.main['HAobsflag']
        #col_names = ['VFID','NED Name','vr','Virgo $d_{3d}$','logMstar','logSFR','logsSFR','H2 def','HI def']#, 'HI def Bos']
        #col_names = ['VFID','NED Name','$v_r$','Virgo $d_{3d}$','$\log(M_\star)$','$\log(SFR)$','$\log(sSFR)$','H2 def','HI def']#, 'HI def Bos']
        tc = [r'\end{center} \n \tablecomments{$^a$ Heliocentric recession velocity. \n']
        tc.append(r'$^b$ Flow-corrected, 3D distance from center of Virgo Cluster from \citet{Castignani2022a} \\')
        tc.append(r'$^c$ Stellar mass from SED fitting. \\')
        tc.append(r'$^d$ SFR from SED fitting. \\')
        tc.append(r'$^e$ sSFR from SED fitting. \\')
        tc.append(r'$^f$ $H_2$ deficiency from \citet{Castignani2022a}.\\')

        # end table comments
        tc.append(r'}\\')

        print(tc)
        col_names = ['VFID','NED Name','$v_r^{a}$','Virgo $d_{3d}^{b}$','$\log(M_\star)^c$',\
                     '$\log(SFR)^d$','$\log(sSFR)^e$','H2 def$^f$']
        col_formats={'$\log(M_\star)^c$': '%5.2f',\
                     '$\log(SFR)^d$': '%5.2f',\
                     '$\log(sSFR)^e$': '%5.2f',\
                     'H2 def$^f$': '%5.2f',\
                     'HI def': '%5.2f',\
                     #'HI def Bos': '%5.2f',\
                     'Virgo $d_{3d}^{b}$':'%5.1f'}
        latexdict={'preamble': r'\begin{center}',\
                   #'tablefoot': tc.append(r'\end{center}'),\
                   'tablefoot': tc,\
                   'tabletype': 'table*',\
                   'header_start': '\\hline \\hline',\
                   'header_end': '\\hline',\
                   'data_end': '\\hline',\
                   'caption': 'NGC~5364 Group Members within the \\ha \\ footprint \\label{tab:sample}',\
                   }
        if outfile is not None:
            latexdict['caption'] = 'NGC~5364 Group Members outside the \\ha \\ footprint \\label{tab:nohalpha}'
        #paperTab = Table([self.main['VFID'],self.main['NEDname'],self.main['vr'],self.dist3dVirgo,self.magphys['logMstar_med'],self.magphys['logSFR_med'],self.magphys['logsSFR_med'],self.paper1['H2def'],self.paper1['HIdef']])[flag]#,self.a100['HIdef_bos']])[flag]

        # test if there is a meerkat HI mass
        HIdef_column = self.paper1['HIdef']
        vindex = np.arange(len(self.main))[flag]

        for i in vindex:
            try:
                logMH = meerkat_HI_mass[v.main['VFID'][i]]
                # get diff b/w 
            except KeyError:
                print(f"WARNING: no meerkat HI mass for {v.main['VFID'][i]}")
                pass

        ################################################
        # replace mass value for VFID5851 with the best value
        ################################################
        mass_column = self.magphys['logMstar_med']
        swap_flag = self.main['VFID'] == 'VFID5851'
        mass_column[swap_flag] = self.magphys['logMstar_best'][swap_flag]

        ################################################
        # construct table
        ################################################

        #paperTab = Table([self.main['VFID'],self.main['NEDname'],self.main['vr'],self.dist3dVirgo,self.magphys['logMstar_med'],self.magphys['logSFR_med'],self.magphys['logsSFR_med'],self.paper1['H2defM'],self.paper1['HIdef']],names=col_names)[flag]
        paperTab = Table([self.main['VFID'],self.main['NEDname'],self.main['vr'],self.dist3dVirgo,mass_column,self.magphys['logSFR_med'],self.magphys['logsSFR_med'],self.paper1['H2def']],names=col_names)[flag]        


        
        ################################################
        # specify units
        ################################################        
        paperTab['$v_r^{a}$'].unit = u.km/u.s
        paperTab['Virgo $d_{3d}^{b}$'].unit = u.Mpc
        paperTab['$\log(M_\star)^c$'].unit = u.Msun
        paperTab['$\log(SFR)^d$'].unit = u.Msun/u.yr
        paperTab['$\log(sSFR)^e$'].unit = u.Msun/u.yr/u.Msun        


        ################################################
        # write table
        ################################################        
        
        if outfile is None:
            outfile = plotdir+'/NGC5364_tab1.tex'
        paperTab.write(outfile,format='latex',names=col_names,formats=col_formats,latexdict=latexdict,\
                       fill_values=[(ascii.masked,'\\nodata')])#,(np.ma.masked,'no data')])
        #paperTab.write(plotdir+'/NGC5364_tab1.tex',format='latex',names=col_names,formats=col_formats,latexdict=ascii.latex.latexdicts['ApJ'])

        if outfile is None:
            self.paperTab = paperTab
        else:
            return paperTab
        pass


    def write_latex_table_nohalpha(self):
        self.get_distance_Virgo()
        flag = self.groupMembs & ~self.main['HAobsflag']
        self.paperTab_nohalpha = self.write_latex_table(flag=flag,outfile=plotdir+'/NGC5364_tab2.tex')

    
    def get_distance_Virgo(self):
        """ calculate 3D distance to Virgo for group members """
        # environment table has distSGX_Virgo, SGY, SGZ
        self.dist3dVirgo = np.sqrt(self.env['distSGX_Virgo']**2+self.env['distSGY_Virgo']**2+self.env['distSGZ_Virgo']**2)
        dist3dVirgo = np.sqrt(self.env['distSGX_Virgo']**2+self.env['distSGY_Virgo']**2+self.env['distSGZ_Virgo']**2)[self.groupMembs]        
        self.dist3dVirgo_mean = np.mean(dist3dVirgo)
        self.dist3dVirgo_med = np.median(dist3dVirgo)
        self.dist3dVirgo_std = np.std(dist3dVirgo)        
        print(f"distance to Virgo = {self.dist3dVirgo_mean:.1f} ({self.dist3dVirgo_med:.1f}) +/- {self.dist3dVirgo_std:.1f}")

    def get_biweight_location_scale(self):
        """ calculate 3D distance to Virgo for group members """
        # environment table has distSGX_Virgo, SGY, SGZ

        self.bw_location = apstats.biweight_location(self.env['Vcosmic'][self.groupMembs])
        self.bw_scale = apstats.biweight.biweight_scale(self.env['Vcosmic'][self.groupMembs])
        print(f"Vcosmic: biweight location (scale) = {self.bw_location:.1f} ({self.bw_scale:.1f})")


        self.bw_location = apstats.biweight_location(self.main['vr'][self.groupMembs])
        self.bw_scale = apstats.biweight.biweight_scale(self.main['vr'][self.groupMembs])
        print(f"vr: biweight location (scale) = {self.bw_location:.1f} ({self.bw_scale:.1f})")
        
    def fit_profiles(self):
        """ fit profile of each group member """
        pass


class fullsample(vtables):
    def plot_HImass_mstar(self):
        '''
        show 
        '''
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
    v.get_group_members()
    v.get_mass_weighted_center()    
    v.get_distance_Virgo()
    v.get_biweight_location_scale()
    v.get_dr_rvirial()
    groupMembs = v.kt['PGC1'] == 49547
