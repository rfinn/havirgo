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
import matplotlib.transforms as transforms
from matplotlib.offsetbox import AnchoredOffsetbox, AuxTransformBox
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
from astropy.coordinates import SkyCoord

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

# rms in the halpha images, measuring this using a box
halpha_rms = {
            'VFID5855': (1.29e-6),
            'VFID5842': (1.15e-6),\
            'VFID5709': (6.50e-7),\
            'VFID6018': (1.80e-6),\
            'VFID6033': (1.60e-6),\
            'VFID6091': (1.35e-6),\
            'VFID6362': (1.25e-6),\
            'VFID2822': 2.03e-6,\
            'VFID2140': 4e-6,\
            'VFID3574': 1.15e-5,
            
           }

# first estimate of CO rms
# 'VFID2822': 0.06234,\
co_rms = {
            'VFID5855': 0.025,
            'VFID5842': 0.025,\
            'VFID5709': 0.025,\
            'VFID6018': 0.025,\
            'VFID6033': 0.025,\
            'VFID6091': 0.025,\
            'VFID6362': 0.025,\
            #0.0623,\ - prior value for VFID2140
            # PGC025063
            #'VFID2822': 0.070, \
            #'VFID2822': 0.0623,\
            'VFID2822': 0.04,\
            #'VFID2822': 0.05,\
            # .082 - prior value for VFID2140
            # rms in Jy/beam/km/s +/- 200 km/s
            #'VFID2140': 0.064,\ # IC4336
            #'VFID2140': 0.1020,\
            'VFID2140': 0.04,\
            #'VFID2140': 0.04,\
            # NGC6186
            #'VFID3574': 0.0479
            #'VFID3574': 0.0694
            'VFID3574': 0.05
           }

rms_dir = homedir+'/research/Virgo/alma/rms_maps/'
co_rms_map = {
            'VFID5855': f"{rms_dir}rms_map_n5348.fits",
            'VFID5842': f"{rms_dir}rms_map_n5356.fits",\
            'VFID5709': f"{rms_dir}rms_map_n5470.fits",\
            'VFID6018': f"{rms_dir}rms_map_n5560.fits",\
            'VFID6033': f"{rms_dir}rms_map_n5566.fits",\
            'VFID6091': f"{rms_dir}rms_map_n5577.fits",\
            'VFID6362': f"{rms_dir}rms_map_u9661.fits",\
           }
    


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
                   'VFID2822':[50,135,50,130],\
                   'VFID2140':[225,390,180,440],\
                       }
zoom_coords_HDI = {'VFID5889':[250,1200,750,1700],\
                   'VFID5842':[400,1300,550,1130],\
                   'VFID5892':[200,600,190,550],\
                   'VFID5859':[35,115,25,115],\
                   'VFID3574':[175,445,175,420]
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
            'VFID2140':[10,8],\
            'VFID2822':[10,8],\
            'VFID3574':[10,8],\
            
            
           }
                     
codepletion_figsize = {
            'VFID5855':[6,8],\
            'VFID5842':[6,10],\
            'VFID5709':[10,5],\
            'VFID6018':[9,5.2],\
            'VFID6033':[9,8],\
            'VFID6091':[9,5.25],\
            'VFID6362':[9,6.5],\
            'VFID2140':[10,8],\
            'VFID2822':[10,8],\
            'VFID3574':[10,8],\
            
           }

codepletion_levels = {
            'VFID5855':[7.,9],\
            'VFID5842':[8.5,10.5],\
            'VFID5709':[9.,11],\
            'VFID6018':[7.5,11],\
            'VFID6033':[8.,9.75],\
            'VFID6091':[7.,9.75],\
            'VFID6362':[7.,9.75],\
            'VFID2140':[8.75,10.],\
            'VFID2822':[8.75, 10.],\
            'VFID3574':[8.75, 10.],\
           }

codepletion_zoom = {
            'VFID5855':[0.55,0.6],\
            'VFID5842':[0.5,0.6],\
            'VFID5709':[0.65,0.6],\
            'VFID6018':[0.5,0.4],\
            'VFID6033':[0.6,0.65],\
            'VFID6091':[0.6,0.6],\
            'VFID6362':[0.5,0.5],\
            'VFID2140':[0.5,0.5],\
            'VFID2822':[.3,.3],\
            'VFID3574':[.65,.6],\
            
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

# width of CO single dish beam
# NGC 5356, NGC 5560, and NGC 5577 >> Boselli+14: Kitt Peak HPBW = 55 arcsec
# otherwise from IRAM - same as for NOEMA paper
CO_single_dish_beam_arcsec = {
                'VFID5842':55,
                'VFID6018':55,
    'VFID6091':55}

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
               'VFID5844':None,\
               'VFID2140':None,\
               'VFID2822':None,\
               'VFID3574':None,\
               
               
}


COdir = homedir+'/research/Virgo/alma/combes_data/'
COveldir = homedir+'/research/Virgo/alma/combes_data/masked_velocity_maps/'
COwidthdir = homedir+'/research/Virgo/alma/combes_data/masked_velocity_width_maps/'
COmaskdir = homedir+'/research/Virgo/alma/combes_data/masked/'
noema_base = homedir+'/research/Virgo/noema/NOEMA_maps/'
noema_base_mask = homedir+'/research/Virgo/noema/NOEMA_maps_masked/'
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
           #'VFID2140':noema_base+'moment_maps_fits_s18cd002/s18cd002-usb-bas-CO10-mean.fits',\
           'VFID2140':noema_base+'moment_maps_fits_s18cd002/s18cd002.fits',\
           'VFID2822':noema_base+'moment_maps_fits_s18cd001/s18cd001.fits',\
           'VFID3574':noema_base+'moment_maps_fits_s18cd003/s18cd003-usb-full-bas-CO10-mean.fits',
           
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
           #'VFID2140':noema_base+'moment_maps_fits_s18cd002/s18cd002_COmap_masked.fits',\
           #'VFID2822':noema_base+'moment_maps_fits_s18cd001/s18cd001_COmap_masked.fits',\
           #'VFID3574':noema_base+'moment_maps_fits_s18cd003/s18cd003_COmap_masked.fits',
           'VFID2140':None,\
           'VFID2822':None,\
           'VFID3574':None

           
}

CO_vel = { 'VFID5709':COveldir+'n5470_velocity_COmap_masked.fits',\
           'VFID5855':COveldir+'n5348_velocity_COmap_masked.fits',\
           'VFID5842':COveldir+'n5356_velocity_COmap_masked.fits',\
           'VFID6018':COveldir+'n5560_velocity_COmap_masked.fits',\
           'VFID6033':COveldir+'n5566_velocity_COmap_masked.fits',\
           'VFID6091':COveldir+'n5577_velocity_COmap_masked.fits',\
           'VFID6362':COveldir+'u9661_velocity_COmap_masked.fits',\
           'VFID2140':noema_base+'moment_maps_fits_s18cd002/s18cd002-usb-bas-CO10-velo.fits',\
           'VFID2822':noema_base+'moment_maps_fits_s18cd001/s18cd001-usb-bas-CO10-velo.fits',\
           'VFID3574':noema_base+'moment_maps_fits_s18cd003/s18cd003-usb-full-bas-CO10-velo.fits',
           
}
CO_width = {'VFID5709':COwidthdir+'n5470_width_COmap_masked.fits',\
           'VFID5855':COwidthdir+'n5348_width_COmap_masked.fits',\
           'VFID5842':COwidthdir+'n5356_width_COmap_masked.fits',\
           'VFID6018':COwidthdir+'n5560_width_COmap_masked.fits',\
           'VFID6033':COwidthdir+'n5566_width_COmap_masked.fits',\
           'VFID6091':COwidthdir+'n5577_width_COmap_masked.fits',\
           'VFID6362':COwidthdir+'u9661_width_COmap_masked.fits',\
           'VFID2140':noema_base+'moment_maps_fits_s18cd002/s18cd002-usb-bas-CO10-width.fits',\
           'VFID2822':noema_base+'moment_maps_fits_s18cd001/s18cd001-usb-bas-CO10-width.fits',\
           'VFID3574':noema_base+'moment_maps_fits_s18cd003/s18cd003-usb-full-bas-CO10-width.fits',
           
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
           'VFID2140':[12,5],\
           'VFID2822':[14,5],\
           'VFID3574':[14,5],\
            
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
           'VFID6362':[3.9],\
           'VFID2140':[4.5],\
           'VFID2822':[4.5],\
           'VFID3574':[5],\

           
           }



HIfiles = {'VFID5859':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_5_mom0.fits',\
           'VFID5855':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_6_mom0.fits',\
           'VFID5851':'research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_2_mom0.fits',\
           }

###########################################
## FUNCTIONS
###########################################
def add_scale(
    ax,
    vr=1100,
    pscale=.331,
    barsize=5.,
    color='k',
    fontsize=16,
    xscale=.09,
    yscale=.9,
    dytext=.08,
    capfrac=0.18,
    lw=2,
):
    '''
    ax = axis for drawing on
    vr = cosmic recession velocity
    pscale = pixel scale, in arcsec/pix
    barsize = size of marker in kpc
    fontsize = fontsize for the text
    capfrac = endcap height as a fraction of the bar length in pixels
    lw = line width
    '''

    from astropy.cosmology import WMAP9 as cosmo
    import numpy as np

    z = vr / 3.e5

    # angular diameter distance in Mpc
    add = cosmo.angular_diameter_distance(z)

    # convert to kpc/arcsec
    add_kpc_arcsec = add.value * 1000 * np.pi / (180 * 3600)

    # convert barsize from kpc to arcsec
    barsize_arcsec = barsize / add_kpc_arcsec

    # convert barsize from arcsec to pixels
    barsize_pixels = barsize_arcsec / pscale

    # get current image limits
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()

    xline1 = x1 + xscale * (x2 - x1)
    xline2 = xline1 + barsize_pixels

    yline = y1 + yscale * (y2 - y1)

    # endcap height in pixels
    capsize_pixels = capfrac * barsize_pixels
    ycap1 = yline - 0.5 * capsize_pixels
    ycap2 = yline + 0.5 * capsize_pixels

    # label position
    xtext = xline1 + 0.5 * barsize_pixels
    ytext = yline - dytext * (y2 - y1)

    # horizontal bar
    ax.plot(
        [xline1, xline2],
        [yline, yline],
        '-',
        lw=lw,
        color=color,
    )

    # left endcap
    ax.plot(
        [xline1, xline1],
        [ycap1, ycap2],
        '-',
        lw=lw,
        color=color,
    )

    # right endcap
    ax.plot(
        [xline2, xline2],
        [ycap1, ycap2],
        '-',
        lw=lw,
        color=color,
    )

    # label
    ax.text(
        xtext,
        ytext,
        '{:.0f} kpc'.format(barsize),
        horizontalalignment='center',
        fontsize=fontsize,
        color=color,
    )
    
def add_scale_original(ax,vr=1100,pscale=.331,barsize=5.,color='k',fontsize=16,xscale=.09,yscale=.9,dytext=.08):
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

    for c in inheader.cards:
        field = c[0]
        if field == 'NAXIS':
            outheader.set(field,2)
        if field.endswith('3'):
            continue
        else:
            try:
                outheader.set(field,inheader[field])
            except:
                print(f"WARNING: problem setting header {field}")# to {inheader[field]}")
    #for c in delcards:
    #    for i in range(2):
    #        ci = c+str(i+1)
    #        outheader.set(ci,inheader[ci])

    #addcards = ['BMAJ','BMIN','BPA','OBJECT','BUNIT']
    #for c in addcards:
    #    outheader.set(c,inheader[c])

    
    #outheader['NAXIS'] = 2
    #outheader['Number of WCS axes'] = 2
    return outheader


def convert_alma3d_alma2d(alma3d_image,nhdu=0):
    hdu = fits.open(alma3d_image)
    outdata = hdu[0].data[nhdu]

    # create a new HDU with input data and no header
    newhdu = fits.PrimaryHDU(outdata)
    # build new header
    outheader = convert_alma_header(hdu[0].header,newhdu.header)
    newhdu.header = outheader
    return outdata, outheader


def turn_off_ax_labels(ax):

    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_ticklabel_visible(False)
    lon.set_ticks_visible(False)            
    lat.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)      
    
def add_inset_colorbar(fig, ax, sc, orientation='horizontal',cblabel=None):
    # Use inset_axes to create a new axes for the colorbar
    # The parameters define the position [left, bottom, width, height] in axes coordinates (0 to 1)
    # For positioning like an AnchoredOffsetbox, you can use the 'loc' keyword
    cax = inset_axes(ax, 
                 width="80%",  # width of the colorbar in percentage of parent axes width
                 height="5%", # height of the colorbar in percentage of parent axes height
                 loc='lower center', # location (e.g., 'upper right', 'lower left')
                 borderpad=1 # padding from the border
                )

    # Add the colorbar to the created inset axes
    cb = fig.colorbar(sc, cax=cax, orientation=orientation)

    if cblabel is not None:
        cb.set_label(cblabel, fontsize=12)
        
def add_metallicity_markers(ax,fig):
    Zsolar = 8.69
    # read in file
    mtab = Table.read(homedir+'/research/Virgo/alma/legacy_metallicity/metallicity_subset.fits')

    # convert RA DEC to skycoord

    coords = SkyCoord(mtab['RA']*u.deg,mtab['DEC']*u.deg, frame='icrs')
    Z_keyword = 'METALLICITY_R23'
    Z_keyword = 'METALLICITY_O3N2'
    sc = ax.scatter(coords.ra, coords.dec, c=mtab[Z_keyword]-Zsolar, facecolors='none', edgecolors='face',s=40, transform=ax.get_transform('world'),vmin=-0.2, vmax=0.2)#, c=mtab['METALLICITY_R23'])
    #cb = fig.colorbar(sc, location='bottom', ax=ax)
    add_inset_colorbar(fig, ax, sc, cblabel=r"$\rm \log_{10}(Z/Z_\odot)$")

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
    dA = cosmo.angular_diameter_distance(z)

    # convert to kpc/arcsec
    dA_kpc_arcsec = dA.value*1000*np.pi/(180*3600)
    
    return dA_kpc_arcsec

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

    
def plot_HI_beam(ax,HIfilename,hostim_header,color='white',expandBox=False,halpha=False, fwhm=1.5):
    """
    DESCRIPTION:
    plot a beam on the provided axis

    PARAMS:
    * ax - axis to plot the beam on.
    * HIfilename - filename where we get the beam size using BPA, BMAJ, BMIN (except for halpha, see below)
    * hostim_header - this is the image where we are plotting the beam.  It is used to set the size of the beam.
    * halpha - set to use with halpha images.  hostim_header is ignored in this case.
    * fwhm - provide the fwhm in arcsec when running on halpha
    """
    ###
    # get size of the host image
    ###
    naxis1 = hostim_header['NAXIS1']
    naxis2 = hostim_header['NAXIS2']
    pscale = (wcs.utils.proj_plane_pixel_scales(WCS(hostim_header))*u.deg).to(u.arcsec) # in deg
    print(f"pscale = {pscale}")
    image_x_arcsec = naxis1*pscale[0]
    image_y_arcsec = naxis2*pscale[1] 

    ###
    # get PSF
    ###

    if halpha:
        PA = 90 # in degrees
        a_arcsec = fwhm*u.arcsec
        b_arcsec = a_arcsec
    
    else:
        hdu = fits.open(HIfilename)
        HI_WCS = WCS(hdu[0].header)
        PA = hdu[0].header['BPA']# in degrees
        a_arcsec = (hdu[0].header['BMAJ']*u.deg).to(u.arcsec)
        b_arcsec = (hdu[0].header['BMIN']*u.deg).to(u.arcsec)
        hdu.close()
    
    #patch_height = (cube_params['bmaj'] / opt_view).decompose()
    #patch_width = (cube_params['bmin'] / opt_view).decompose()
    patch_height = a_arcsec/pscale[0] # in pix
    patch_width = b_arcsec/pscale[1] # in pix
    patch = {'width': patch_width, 'height': patch_height}
    print(f"beam size info: a={a_arcsec:.3f},b={b_arcsec:.3f},image_x_arcsec={image_x_arcsec:.1f},pa={PA:.2f},patch_height={patch['height']:.4f}, patch_width={patch['width']:.4f}")
    #ax.add_patch(Ellipse((0.9, 0.9), height=patch['height'], width=patch['width'], angle=PA,transform=ax.transAxes, edgecolor=color,facecolor=color, linewidth=1,alpha=.8))    

    '''
    # add white box to highlight position of beam
    dx = 0.1#10*u.arcsec/image_x_arcsec
    dy = 0.1#10*u.arcsec/image_y_arcsec
    if expandBox:
        dx=.25
        dy=.25
    xc = .9
    yc = 0.9

    if xcenter_box is not None:
        dx = 10*u.arcsec/image_x_arcsec # box width in pixels
        dy = 10*u.arcsec/image_y_arcsec # box height in pixels
        # describe the box in pixel coordinates
        rect = Rectangle((xcenter_box,ycenter_box),dx,dy,edgecolor=color,facecolor='None',lw=1,alpha=.8)
    else:
        rect = Rectangle((xc-0.5*dx,yc-0.5*dy),dx,dy,transform=ax.transAxes,edgecolor=color,facecolor='None',lw=1,alpha=.8)
    #ax.add_patch(rect)
    '''

    # found this https://matplotlib.org/stable/gallery/misc/anchored_artists.html#sphx-glr-gallery-misc-anchored-artists-py
    aux_tr_box = AuxTransformBox(ax.transData)
    myellipse = Ellipse((0., 0.), height=patch['height'], width=patch['width'], angle=PA,\
                         edgecolor=color,facecolor=color, linewidth=1,alpha=.8)    

    aux_tr_box.add_artist(myellipse)
    box = AnchoredOffsetbox(child=aux_tr_box, loc="upper right", frameon=True)
    ax.add_artist(box)
    


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

def plot_noema_footprint(ax,COfilename, r=22.):
    """ add circle to axis  """
    header = fits.getheader(COfilename)
    
    # get RA and DEC from header
    ra = header['CRVAL1']
    dec = header['CRVAL2']

    # convert radius to deg
    rad = r/3600.

    # for testing
    plt.plot(ra,dec,'cX')
    # make a circle
    c = plt.Circle((ra,dec),rad)
    # add to axis
    ax.add_patch(c)
    
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

def plot_spines(multicolor=True,colorone=None,color=None,legend=True,plotone=None):
    import glob
    sfiles = glob.glob(homedir+'/research/Virgo/tables-north/spines/filament*.fits')
    ncolor = 0
    for i,f in enumerate(sfiles):
        if (plotone is not None) and (plotone not in f):
            continue
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


def get_image_footprint(image_name):
    hdu = fits.open(image_name)[0]
    wcs = WCS(hdu.header)
    # Calculate the celestial footprint (4 corners)
    footprint = wcs.calc_footprint() # Returns (4, 2) array of RA/Dec
    return footprint

def plot_footprint(ax, footprint, color='c', label="footprint"):
    ax.plot(footprint[:, 0], footprint[:, 1], 'r-', color=color, label=label)


def get_CO_footprint(ax, COfilename):
    try:
        codata, coheader = convert_alma3d_alma2d(COfilename)
    except KeyError:
        codata, coheader = fits.getdata(COfilename, header=True)
    #print(coheader)
    CO_WCS = WCS(coheader)
    footprint = CO_WCS.calc_footprint()
    print("DEBUG: footprint = ",footprint)
    return footprint

def plot_single_dish_CO_beam(ax, vfid,jpeg_data, legwcs):

    from astropy.visualization.wcsaxes import SphericalCircle
    from astropy.coordinates import SkyCoord
    import astropy.units as u


    try:
        co_single_radius_arcsec = CO_single_dish_beam_arcsec[vfid]
    except KeyError:
        co_single_radius_arcsec = 22.
    rad = co_single_radius_arcsec * u.arcsec
    # center of image from WCS
    #ny, nx = np.shape(jpeg_data)[0], np.shape(jpeg_data)[1]
    #xc, yc = nx/2, ny/2
    #ra_c, dec_c = legwcs.wcs_pix2world(xc, yc, 0)

    vflag = v.main['VFID'] == vfid

    ra_c = v.main['RA'][vflag][0]
    dec_c = v.main['DEC'][vflag][0]      
    
    center = SkyCoord(ra_c*u.deg, dec_c*u.deg, frame='icrs')
    
    # 22 arcsec circle
    circ = SphericalCircle(
        center,
        rad,
        edgecolor='cyan',
        facecolor='none',
        linewidth=2,
        transform=ax.get_transform('world')
        )
    ax.add_patch(circ)
    
def plot_mstar_sfr_CO(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,6],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,\
                            logMstar=None,cmap='magma_r',markGroupCenter=False,COcolor='white',COlevels=None,\
                          noHI=False,alt_title=None,vr=None, sfr_limits=None,noema=False):
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
    if noHI:
        titles = [r'$\log_{10}(\Sigma_\star)$',r'$\Sigma_{SFR}$',r'$\Sigma_{SFR} \ w/CO$']
    else:
        titles = [r'$\log_{10}(\Sigma_\star)$',r'$\Sigma_{SFR} \ w/HI$',r'$\Sigma_{SFR} \ w/CO$']

    #vmin = [2,0e-5,0e-5,0,-200,0]
    #vmax = [6,.6e-5,.6e-5,2,200,100]

    vmin = [4,0e-5,0e-5]
    vmax = [7,3e-5,3e-5]


    
    allim = [massim,sfrim,ssfrim]
    #allim = [massim,sfrim]
    allim = [massim,sfrim,sfrim] # show halpha with HI and CO    
    cblabels = [r'$\rm \log_{10}(M_\star/M_\odot)/pixel$',r'$\rm SFR(M_\star/yr)/pixel$',r'$\rm SFR(M_\star/yr)/pixel$']


    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)
    if ('INT' in dirname) | ('BOK' in dirname):
        try:
            xmin,xmax,ymin,ymax = zoom_coords_INT[vfid]
        except KeyError:
            # use the full image
            # read in massim
            print("WARNING: no zoom coords found ",dirname)
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
        
        if i > 0:              
            if vr is not None:
                # scale data to convert to flux/kpc^2
                kpc_per_pixel = get_kpc_per_pixel(imwcs,vr)
                dat = dat/(kpc_per_pixel**2)
                vmin[i] = 0 # vmin[i]/kpc_per_pixel**2 
                vmax[i] = 0.02#1e-5/kpc_per_pixel**2
                cblabels[i] = r'$\rm SFR(M_\star/yr)/kpc^2$'


        mdat = np.ma.array(dat,mask=maskdat)
        #if xmin is None:
        #    mdat = mdat
        #else:
        #    mdat = mdat[ymin:ymax,xmin:xmax]
        if i == 2:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            if sfr_limits is not None:
                v1, v2 = sfr_limits
            else:
                v1, v2 = vmin[i], vmax[i]
            plt.imshow(mdat,vmin=v1,vmax=v2,origin='lower',interpolation='nearest',cmap='gray_r')


        elif i == 1:
            #display_image(mdat,percent=99.5,cmap='viridis')#,vmin=vmin[i],vmax=vmax[i])
            if sfr_limits is not None:
                v1, v2 = sfr_limits
            else:
                v1, v2 = vmin[i], vmax[i]
            plt.imshow(mdat,vmin=v1,vmax=v2,cmap=cmap)#cmap='viridis'
        else:
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

    #add_metallicity_markers(ax1, fig)

    # check if HI moment zero map is available
    COfilename = CO_file[vfid]
    COmask_filename = CO_mask[vfid]
    plot_CO_contours(ax1,COfilename,color=COcolor,mask=COmask_filename,levels=COlevels)

    # plot CO footprint
    #footprint = get_CO_footprint(ax1, COfilename)
    #footcoords = SkyCoord(footprint, unit=u.deg)
    #print("DEBUG: footcoords.ra",footcoords.ra.hms)
    #plt.scatter(footcoords.ra.deg, footcoords.dec.deg)
    #ax1.set_xlabel("RA",fontsize=16)
    #ax1.set_ylabel("DEC",fontsize=16)    

    # add noema primary beam of 22" circle
    #plot_single_dish_CO_beam(ax1, vfid, jpeg_data, legwcs)
    
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
    if alt_title is not None:
        plt.title(alt_title,fontsize=18)
    else:
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
            if CO_mask[vfid] is not None:
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
                dat = cdat
        
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
            cbaxes = inset_axes(ax2, width="80%", height="5%", loc='lower center')#loc=8)
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
    turn_off_ax_labels(ax1)
    add_metallicity_markers(ax1, fig)

    ax5 = plt.subplot(nrow,ncol,5,projection=legwcs)
    plt.imshow(jpeg_data)    


    plot_single_dish_CO_beam(ax5, vfid, jpeg_data, legwcs)    
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

def hide_xyticks_wcs(ax):
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_ticklabel_visible(False)
    lon.set_ticks_visible(False)            
    lat.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)         
    
from scipy import ndimage as ndi

def clean_limit_mask_for_contours(mask, min_pixels=10, do_closing=True):
    """
    Clean a boolean mask for display contours only.

    Removes small isolated regions and optionally smooths/fills small gaps.
    This should not be used for the actual depletion map calculation.
    """
    mask = np.asarray(mask, dtype=bool)

    # Label connected regions
    structure = ndi.generate_binary_structure(2, 2)  # 8-connected
    labels, nlab = ndi.label(mask, structure=structure)

    if nlab == 0:
        return mask

    # Remove small islands
    counts = np.bincount(labels.ravel())
    keep = counts >= min_pixels
    keep[0] = False  # background

    cleaned = keep[labels]

    # Optional cosmetic cleanup of jagged holes/gaps
    if do_closing:
        cleaned = ndi.binary_closing(cleaned, structure=structure, iterations=1)

    return cleaned


def add_limit_contours(
    ax,
    upper_mask,
    lower_mask,
    x1,
    x2,
    y1,
    y2,
    show_upper=True,
    show_lower=True,
    upper_color="magenta",
    lower_color="cyan",
    upper_ls="solid",
    lower_ls="solid",
    lw=1.2,
    alpha=0.9,
):
    """
    Add depletion-time upper/lower-limit contours to an image panel.

    upper_mask:
        Halpha detected, CO not detected.
    lower_mask:
        CO detected, Halpha not detected.
    """

    if show_upper:
        upper_zoom = upper_mask[y1:y2, x1:x2].astype(float)

        if np.any(upper_zoom):
            ax.contour(
                upper_zoom,
                levels=[0.5],
                colors=upper_color,
                linestyles=upper_ls,
                linewidths=lw,
                alpha=alpha,
            )

    if show_lower:
        lower_zoom = lower_mask[y1:y2, x1:x2].astype(float)

        if np.any(lower_zoom):
            ax.contour(
                lower_zoom,
                levels=[0.5],
                colors=lower_color,
                linestyles=lower_ls,
                linewidths=lw,
                alpha=alpha,
            )

def get_ks_points_one_galaxy(
    dirname,
    vr,
    nsigma=3,
    alma=False,
    rootdir=None,
    H0=74.,
    return_mask=False,
    verbose=True,
):
    """
    Extract resolved Kennicutt-Schmidt points for one galaxy.

    Returns pixels/regions where:
      - CO is detected with SNR > nsigma
      - Halpha/SFR is detected with SNR > nsigma
      - both Sigma_H2 and Sigma_SFR are finite and positive

    Assumes get_depletion_map has already created:
      dirname-Sigma-H2-kpc2.fits
      dirname-Sigma-sfr-kpc2.fits

    Parameters
    ----------
    dirname : str
        Galaxy cutout directory name.
    vr : float
        Recession velocity used for distance/kpc conversion.
    nsigma : float
        Detection threshold. Default is 3.
    alma : bool
        If True, use co_rms_map[vfid]. If False, use scalar co_rms[vfid].
    rootdir : str or None
        Parent directory containing dirname. If None, use current directory.
    H0 : float
        Hubble constant used for D_Mpc = vr/H0.
    return_mask : bool
        If True, also return the boolean detected-both mask.
    verbose : bool
        Print summary information.

    Returns
    -------
    result : dict
        Dictionary containing sigma_H2, sigma_sfr, log_sigma_H2,
        log_sigma_sfr, galaxy, vfid, and optionally detected_both_mask.
    """

    import os
    import numpy as np
    from astropy.io import fits
    from astropy import wcs
    from astropy.wcs import WCS

    vfid = dirname.split("-")[0]

    if rootdir is None:
        galpath = dirname
    else:
        galpath = os.path.join(rootdir, dirname)

    sigma_H2_file = os.path.join(galpath, f"{dirname}-Sigma-H2-kpc2.fits")
    sigma_sfr_file = os.path.join(galpath, f"{dirname}-Sigma-sfr-kpc2.fits")

    if not os.path.exists(sigma_H2_file):
        raise FileNotFoundError(f"Missing Sigma_H2 file: {sigma_H2_file}")

    if not os.path.exists(sigma_sfr_file):
        raise FileNotFoundError(f"Missing Sigma_SFR file: {sigma_sfr_file}")

    sigma_H2, h2_header = fits.getdata(sigma_H2_file, header=True)
    sigma_sfr, sfr_header = fits.getdata(sigma_sfr_file, header=True)

    if sigma_H2.shape != sigma_sfr.shape:
        raise ValueError(
            f"Shape mismatch: Sigma_H2 has shape {sigma_H2.shape}, "
            f"but Sigma_SFR has shape {sigma_sfr.shape}"
        )

    co_imwcs = WCS(h2_header)
    kpc_per_pixel = get_kpc_per_pixel(co_imwcs, vr)

    ##################################################
    # CO/H2 SNR threshold in Msun/kpc^2
    ##################################################

    D_Mpc = vr / H0
    alpha_CO = 4.3

    pixscale_deg = wcs.utils.proj_plane_pixel_scales(co_imwcs)[0]
    pixel_area_sq_deg = pixscale_deg**2

    beam_area_sq_deg = np.pi / 4 * float(h2_header["BMAJ"]) * float(h2_header["BMIN"])
    scale_factor = beam_area_sq_deg / pixel_area_sq_deg

    if alma:
        # rms map is assumed to be Jy/beam km/s on the CO grid
        rms_map = fits.getdata(co_rms_map[vfid])

        co_rms_Msun_beam = (
            rms_map
            * 1.05e4
            * alpha_CO / 4.3
            * D_Mpc**2
        )

        co_rms_Msun_pixel = co_rms_Msun_beam / np.sqrt(scale_factor)
        sigma_H2_rms = co_rms_Msun_pixel / kpc_per_pixel**2

    else:
        co_rms_Msun_beam = (
            co_rms[vfid]
            * 1.05e4
            * alpha_CO / 4.3
            * D_Mpc**2
        )

        co_rms_Msun_pixel = co_rms_Msun_beam / np.sqrt(scale_factor)
        sigma_H2_rms = co_rms_Msun_pixel / kpc_per_pixel**2

    co_detected = sigma_H2 > nsigma * sigma_H2_rms

    ##################################################
    # Halpha/SFR SNR threshold in Msun/yr/kpc^2
    ##################################################

    sigma_sfr_rms = halpha_rms[vfid] / kpc_per_pixel**2
    ha_detected = sigma_sfr > nsigma * sigma_sfr_rms

    ##################################################
    # Combined valid KS mask
    ##################################################

    detected_both = (
        co_detected
        & ha_detected
        & np.isfinite(sigma_H2)
        & np.isfinite(sigma_sfr)
        & (sigma_H2 > 0)
        & (sigma_sfr > 0)
    )

    sigma_H2_out = sigma_H2[detected_both]/1e6
    sigma_sfr_out = sigma_sfr[detected_both]


    result = {
        "galaxy": dirname,
        "vfid": vfid,

        # original internal units
        "sigma_H2_kpc2": sigma_H2_out,
        "sigma_sfr_kpc2": sigma_sfr_out,

        # plotting units
        "sigma_H2_pc2": sigma_H2_pc2_out,
        "sigma_sfr": sigma_sfr_out,

        "log_sigma_H2_pc2": np.log10(sigma_H2_pc2_out),
        "log_sigma_sfr": np.log10(sigma_sfr_out),
        "n_points": len(sigma_H2_out),
        }

    if return_mask:
        result["detected_both_mask"] = detected_both

    if verbose:
        print(
            f"{dirname}: {len(sigma_H2_out)} pixels with "
            f"CO SNR>{nsigma} and Halpha SNR>{nsigma}"
        )

        if len(sigma_H2_out) > 0:
            print(
                f"  log Sigma_H2 range: "
                f"{np.nanmin(result['log_sigma_H2']):.2f} to "
                f"{np.nanmax(result['log_sigma_H2']):.2f}"
            )
            print(
                f"  log Sigma_SFR range: "
                f"{np.nanmin(result['log_sigma_sfr']):.2f} to "
                f"{np.nanmax(result['log_sigma_sfr']):.2f}"
            )

    return result

def _trim_for_blocks(arr, block_size, yoff=0, xoff=0):
    """
    Trim a 2D array so it can be reshaped into non-overlapping square blocks.
    """

    arr = arr[yoff:, xoff:]

    ny, nx = arr.shape
    ny_trim = (ny // block_size) * block_size
    nx_trim = (nx // block_size) * block_size

    return arr[:ny_trim, :nx_trim]


def _block_sum(arr, block_size, good=None, yoff=0, xoff=0):
    """
    Sum a 2D array in non-overlapping block_size x block_size regions.

    If good is supplied, values outside good are treated as zero.
    """

    import numpy as np

    arr_trim = _trim_for_blocks(arr, block_size, yoff=yoff, xoff=xoff)

    if good is not None:
        good_trim = _trim_for_blocks(good, block_size, yoff=yoff, xoff=xoff)
        arr_trim = np.where(good_trim, arr_trim, 0.0)

    ny, nx = arr_trim.shape
    nyb = ny // block_size
    nxb = nx // block_size

    return arr_trim.reshape(nyb, block_size, nxb, block_size).sum(axis=(1, 3))


def _block_count(good, block_size, yoff=0, xoff=0):
    """
    Count True pixels in non-overlapping block_size x block_size regions.
    """

    return _block_sum(
        good.astype(float),
        block_size,
        good=None,
        yoff=yoff,
        xoff=xoff,
    )

def get_ks_points_one_galaxy(
    dirname,
    vr,
    nsigma=3,
    alma=False,
    rootdir=None,
    H0=74.,
    return_mask=False,
    verbose=True,
):
    """
    Extract resolved Kennicutt-Schmidt points for one galaxy.

    Returns pixels/regions where:
      - CO is detected with SNR > nsigma
      - Halpha/SFR is detected with SNR > nsigma
      - both Sigma_H2 and Sigma_SFR are finite and positive

    Assumes get_depletion_map has already created:
      dirname-Sigma-H2-kpc2.fits
      dirname-Sigma-sfr-kpc2.fits

    Parameters
    ----------
    dirname : str
        Galaxy cutout directory name.
    vr : float
        Recession velocity used for distance/kpc conversion.
    nsigma : float
        Detection threshold. Default is 3.
    alma : bool
        If True, use co_rms_map[vfid]. If False, use scalar co_rms[vfid].
    rootdir : str or None
        Parent directory containing dirname. If None, use current directory.
    H0 : float
        Hubble constant used for D_Mpc = vr/H0.
    return_mask : bool
        If True, also return the boolean detected-both mask.
    verbose : bool
        Print summary information.

    Returns
    -------
    result : dict
        Dictionary containing sigma_H2, sigma_sfr, log_sigma_H2,
        log_sigma_sfr, galaxy, vfid, and optionally detected_both_mask.
    """

    import os
    import numpy as np
    from astropy.io import fits
    from astropy import wcs
    from astropy.wcs import WCS

    vfid = dirname.split("-")[0]

    if rootdir is None:
        galpath = dirname
    else:
        galpath = os.path.join(rootdir, dirname)

    sigma_H2_file = os.path.join(galpath, f"{dirname}-Sigma-H2-kpc2.fits")
    sigma_sfr_file = os.path.join(galpath, f"{dirname}-Sigma-sfr-kpc2.fits")

    if not os.path.exists(sigma_H2_file):
        raise FileNotFoundError(f"Missing Sigma_H2 file: {sigma_H2_file}")

    if not os.path.exists(sigma_sfr_file):
        raise FileNotFoundError(f"Missing Sigma_SFR file: {sigma_sfr_file}")

    sigma_H2, h2_header = fits.getdata(sigma_H2_file, header=True)
    sigma_sfr, sfr_header = fits.getdata(sigma_sfr_file, header=True)

    if sigma_H2.shape != sigma_sfr.shape:
        raise ValueError(
            f"Shape mismatch: Sigma_H2 has shape {sigma_H2.shape}, "
            f"but Sigma_SFR has shape {sigma_sfr.shape}"
        )

    co_imwcs = WCS(h2_header)
    kpc_per_pixel = get_kpc_per_pixel(co_imwcs, vr)

    ##################################################
    # CO/H2 SNR threshold in Msun/kpc^2
    ##################################################

    D_Mpc = vr / H0
    alpha_CO = 4.3

    pixscale_deg = wcs.utils.proj_plane_pixel_scales(co_imwcs)[0]
    pixel_area_sq_deg = pixscale_deg**2

    beam_area_sq_deg = np.pi / 4 * float(h2_header["BMAJ"]) * float(h2_header["BMIN"])
    scale_factor = beam_area_sq_deg / pixel_area_sq_deg

    if alma:
        # rms map is assumed to be Jy/beam km/s on the CO grid
        rms_map = fits.getdata(co_rms_map[vfid])

        co_rms_Msun_beam = (
            rms_map
            * 1.05e4
            * alpha_CO / 4.3
            * D_Mpc**2
        )

        co_rms_Msun_pixel = co_rms_Msun_beam / np.sqrt(scale_factor)
        sigma_H2_rms = co_rms_Msun_pixel / kpc_per_pixel**2

    else:
        co_rms_Msun_beam = (
            co_rms[vfid]
            * 1.05e4
            * alpha_CO / 4.3
            * D_Mpc**2
        )

        co_rms_Msun_pixel = co_rms_Msun_beam / np.sqrt(scale_factor)
        sigma_H2_rms = co_rms_Msun_pixel / kpc_per_pixel**2

    co_detected = sigma_H2 > nsigma * sigma_H2_rms

    ##################################################
    # Halpha/SFR SNR threshold in Msun/yr/kpc^2
    ##################################################

    sigma_sfr_rms = halpha_rms[vfid] / kpc_per_pixel**2
    ha_detected = sigma_sfr > nsigma * sigma_sfr_rms

    ##################################################
    # Combined valid KS mask
    ##################################################

    detected_both = (
        co_detected
        & ha_detected
        & np.isfinite(sigma_H2)
        & np.isfinite(sigma_sfr)
        & (sigma_H2 > 0)
        & (sigma_sfr > 0)
    )

    sigma_H2_out = sigma_H2[detected_both]
    sigma_sfr_out = sigma_sfr[detected_both]

    result = {
        "galaxy": dirname,
        "vfid": vfid,
        "sigma_H2": sigma_H2_out,
        "sigma_sfr": sigma_sfr_out,
        "log_sigma_H2": np.log10(sigma_H2_out),
        "log_sigma_sfr": np.log10(sigma_sfr_out),
        "n_points": len(sigma_H2_out),
    }

    if return_mask:
        result["detected_both_mask"] = detected_both

    if verbose:
        print(
            f"{dirname}: {len(sigma_H2_out)} pixels with "
            f"CO SNR>{nsigma} and Halpha SNR>{nsigma}"
        )

        if len(sigma_H2_out) > 0:
            print(
                f"  log Sigma_H2 range: "
                f"{np.nanmin(result['log_sigma_H2']):.2f} to "
                f"{np.nanmax(result['log_sigma_H2']):.2f}"
            )
            print(
                f"  log Sigma_SFR range: "
                f"{np.nanmin(result['log_sigma_sfr']):.2f} to "
                f"{np.nanmax(result['log_sigma_sfr']):.2f}"
            )

    return result

def get_ks_points_one_galaxy_binned(
    dirname,
    vr,
    nsigma=3,
    alma=False,
    rootdir=None,
    H0=74.,
    binsize_pix=None,
    binsize_kpc=None,
    min_valid_frac=0.5,
    yoff=0,
    xoff=0,
    return_all_bins=False,
    verbose=True,
):
    """
    Extract resolved KS points using approximately beam-sized independent regions.

    The function bins the Sigma_H2 and Sigma_SFR maps into non-overlapping
    square regions whose area is approximately one CO beam area.

    Selection:
      - binned CO SNR > nsigma
      - binned Halpha/SFR SNR > nsigma
      - binned Sigma_H2 > 0
      - binned Sigma_SFR > 0
      - at least min_valid_frac of the pixels in the bin are valid

    Parameters
    ----------
    dirname : str
        Galaxy cutout directory name.
    vr : float
        Cosmic recession velocity used for kpc conversion.
    nsigma : float
        Detection threshold.
    alma : bool
        If True, use co_rms_map[vfid]. If False, use scalar co_rms[vfid].
    rootdir : str or None
        Parent directory containing dirname. If None, use current directory.
    H0 : float
        Hubble constant for D_Mpc = vr/H0.
    binsize_pix : int or None
        Size of square bin in pixels. If None, use sqrt(beam area / pixel area).
    min_valid_frac : float
        Minimum fraction of pixels in a bin that must be finite/usable.
    yoff, xoff : int
        Optional grid offsets in pixels. Useful for checking sensitivity to bin placement.
    return_all_bins : bool
        If True, include full binned maps and SNR maps in the output.
    verbose : bool
        Print summary information.

    Returns
    -------
    result : dict
        Contains sigma_H2, sigma_sfr, log_sigma_H2, log_sigma_sfr,
        co_snr, ha_snr, and metadata.
    """

    import os
    import numpy as np
    from astropy.io import fits
    from astropy import wcs
    from astropy.wcs import WCS

    vfid = dirname.split("-")[0]
    t = dirname.split("-")
    galtitle = t[0]+'-'+t[1]
    if rootdir is None:
        galpath = dirname
    else:
        galpath = os.path.join(rootdir, dirname)

    sigma_H2_file = os.path.join(galpath, f"{dirname}-Sigma-H2-kpc2.fits")
    sigma_sfr_file = os.path.join(galpath, f"{dirname}-Sigma-sfr-kpc2.fits")

    if not os.path.exists(sigma_H2_file):
        raise FileNotFoundError(f"Missing Sigma_H2 file: {sigma_H2_file}")

    if not os.path.exists(sigma_sfr_file):
        raise FileNotFoundError(f"Missing Sigma_SFR file: {sigma_sfr_file}")

    sigma_H2, h2_header = fits.getdata(sigma_H2_file, header=True)
    sigma_sfr, sfr_header = fits.getdata(sigma_sfr_file, header=True)

    if sigma_H2.shape != sigma_sfr.shape:
        raise ValueError(
            f"Shape mismatch: Sigma_H2 has shape {sigma_H2.shape}, "
            f"but Sigma_SFR has shape {sigma_sfr.shape}"
        )

    co_imwcs = WCS(h2_header)
    kpc_per_pixel = get_kpc_per_pixel(co_imwcs, vr)
    pixel_area_kpc2 = kpc_per_pixel**2

    ##################################################
    # Beam area in pixels
    ##################################################

    pixscale_deg = wcs.utils.proj_plane_pixel_scales(co_imwcs)[0]
    pixel_area_sq_deg = pixscale_deg**2

    beam_area_sq_deg = np.pi / 4 * float(h2_header["BMAJ"]) * float(h2_header["BMIN"])

    beam_area_pix = beam_area_sq_deg / pixel_area_sq_deg

    print(f"{dirname}: bin size = {np.sqrt(beam_area_pix)*kpc_per_pixel:.3f} kpc")    

    if binsize_kpc is not None:
        # Use square bins with side length approximately binsize_kpc
        binsize_pix = max(1, int(np.round(binsize_kpc / kpc_per_pixel)))

    elif binsize_pix is None:
        # Default: square bins with area approximately equal to one CO beam
        binsize_pix = max(1, int(np.round(np.sqrt(beam_area_pix))))


    bin_npix = binsize_pix**2
    min_valid_pixels = int(np.ceil(min_valid_frac * bin_npix))

    ##################################################
    # CO/H2 RMS in Sigma units
    ##################################################

    D_Mpc = vr / H0
    alpha_CO = 4.3

    if alma:
        rms_map = fits.getdata(co_rms_map[vfid])

        if rms_map.shape != sigma_H2.shape:
            raise ValueError(
                f"CO rms map shape {rms_map.shape} does not match "
                f"Sigma_H2 shape {sigma_H2.shape}"
            )

        co_rms_Msun_beam = (
            rms_map
            * 1.05e4
            * alpha_CO / 4.3
            * D_Mpc**2
        )

        # Same convention as get_depletion_map:
        # pixel rms = beam rms / sqrt(beam_area_pix)
        co_rms_Msun_pixel = co_rms_Msun_beam / np.sqrt(beam_area_pix)
        sigma_H2_rms = co_rms_Msun_pixel / pixel_area_kpc2

    else:
        co_rms_Msun_beam = (
            co_rms[vfid]
            * 1.05e4
            * alpha_CO / 4.3
            * D_Mpc**2
        )

        co_rms_Msun_pixel = co_rms_Msun_beam / np.sqrt(beam_area_pix)
        sigma_H2_rms = co_rms_Msun_pixel / pixel_area_kpc2

    ##################################################
    # Halpha/SFR RMS
    ##################################################

    # In get_depletion_map, halpha_rms[vfid] is used as SFR/pixel
    # on the reprojected SFR map. Convert the binned sum to SNR using
    # SFR/pixel noise directly.
    sfr_rms_pixel = halpha_rms[vfid]

    ##################################################
    # Valid pixels
    ##################################################

    valid = (
        np.isfinite(sigma_H2)
        & np.isfinite(sigma_sfr)
    )

    ##################################################
    # Convert surface densities to per-pixel mass/SFR
    ##################################################

    H2_mass_pix = sigma_H2 * pixel_area_kpc2
    sfr_pix = sigma_sfr * pixel_area_kpc2

    ##################################################
    # Bin signal
    ##################################################

    H2_mass_bin = _block_sum(
        H2_mass_pix,
        binsize_pix,
        good=valid,
        yoff=yoff,
        xoff=xoff,
    )

    sfr_bin = _block_sum(
        sfr_pix,
        binsize_pix,
        good=valid,
        yoff=yoff,
        xoff=xoff,
    )

    n_valid_bin = _block_count(
        valid,
        binsize_pix,
        yoff=yoff,
        xoff=xoff,
    )

    area_bin_kpc2 = n_valid_bin * pixel_area_kpc2

    with np.errstate(invalid="ignore", divide="ignore"):
        sigma_H2_bin = H2_mass_bin / area_bin_kpc2
        sigma_sfr_bin = sfr_bin / area_bin_kpc2

    ##################################################
    # Bin noise
    ##################################################

    if np.ndim(sigma_H2_rms) == 0:
        H2_mass_err_pix = sigma_H2_rms * pixel_area_kpc2
        H2_mass_err_bin = np.sqrt(n_valid_bin) * H2_mass_err_pix
    else:
        H2_mass_err_pix = sigma_H2_rms * pixel_area_kpc2
        H2_mass_err2_bin = _block_sum(
            H2_mass_err_pix**2,
            binsize_pix,
            good=valid & np.isfinite(H2_mass_err_pix),
            yoff=yoff,
            xoff=xoff,
        )
        H2_mass_err_bin = np.sqrt(H2_mass_err2_bin)

    sfr_err_bin = np.sqrt(n_valid_bin) * sfr_rms_pixel

    with np.errstate(invalid="ignore", divide="ignore"):
        co_snr_bin = H2_mass_bin / H2_mass_err_bin
        ha_snr_bin = sfr_bin / sfr_err_bin

    ##################################################
    # Final detected-both mask
    ##################################################

    detected_both_bin = (
        (n_valid_bin >= min_valid_pixels)
        & np.isfinite(sigma_H2_bin)
        & np.isfinite(sigma_sfr_bin)
        & np.isfinite(co_snr_bin)
        & np.isfinite(ha_snr_bin)
        & (sigma_H2_bin > 0)
        & (sigma_sfr_bin > 0)
        & (co_snr_bin > nsigma)
        & (ha_snr_bin > nsigma)
    )

    ##################################################
    # Bin centers in original pixel coordinates
    ##################################################

    nyb, nxb = sigma_H2_bin.shape

    y_centers = yoff + np.arange(nyb) * binsize_pix + 0.5 * (binsize_pix - 1)
    x_centers = xoff + np.arange(nxb) * binsize_pix + 0.5 * (binsize_pix - 1)
    xbin, ybin = np.meshgrid(x_centers, y_centers)

    ##################################################
    # Package output
    ##################################################
    sigma_H2_out = sigma_H2_bin[detected_both_bin]
    sigma_sfr_out = sigma_sfr_bin[detected_both_bin]

    sigma_H2_pc2_out = sigma_H2_out / 1e6

    result = {
        "galaxy": dirname,
        "vfid": vfid,
        "galtitle": galtitle,
        "sigma_H2_kpc2": sigma_H2_out,
        "sigma_H2_pc2": sigma_H2_pc2_out,
        "sigma_sfr_kpc2": sigma_sfr_out,

        "log_sigma_H2_pc2": np.log10(sigma_H2_pc2_out),
        "log_sigma_sfr": np.log10(sigma_sfr_out),

        "co_snr": co_snr_bin[detected_both_bin],
        "ha_snr": ha_snr_bin[detected_both_bin],
        "xbin": xbin[detected_both_bin],
        "ybin": ybin[detected_both_bin],
        "n_valid_pix": n_valid_bin[detected_both_bin],
        "n_points": len(sigma_H2_out),
        "binsize_pix": binsize_pix,
        "bin_npix": bin_npix,
        "beam_area_pix": beam_area_pix,
        "min_valid_pixels": min_valid_pixels,
        "kpc_per_pixel": kpc_per_pixel,
        }
    

    if return_all_bins:
        result["sigma_H2_bin_map"] = sigma_H2_bin
        result["sigma_sfr_bin_map"] = sigma_sfr_bin
        result["co_snr_bin_map"] = co_snr_bin
        result["ha_snr_bin_map"] = ha_snr_bin
        result["detected_both_bin_map"] = detected_both_bin
        result["n_valid_bin_map"] = n_valid_bin

    if verbose:
        print("########################################")
        print(dirname)
        print("########################################")
        print(f"beam area = {beam_area_pix:.1f} pixels")
        print(f"using binsize_pix = {binsize_pix}")
        print(f"bin area = {bin_npix} pixels")
        print(f"minimum valid pixels per bin = {min_valid_pixels}")
        print(f"binned CO detections:      {np.sum(co_snr_bin > nsigma)}")
        print(f"binned Halpha detections:  {np.sum(ha_snr_bin > nsigma)}")
        print(f"binned detected in both:   {result['n_points']}")

        if result["n_points"] > 0:
            print(
                f"log Sigma_H2 range: "
                f"{np.nanmin(result['log_sigma_H2_pc2']):.2f} to "
                f"{np.nanmax(result['log_sigma_H2_pc2']):.2f}"
            )
            print(
                f"log Sigma_SFR range: "
                f"{np.nanmin(result['log_sigma_sfr']):.2f} to "
                f"{np.nanmax(result['log_sigma_sfr']):.2f}"
            )

    return result

def get_ks_points_many_galaxies(
    dirnames,
    vr_dict,
    nsigma=3,
    alma=False,
    rootdir=None,
    H0=74.,
    verbose=True,
):
    """
    Collect resolved KS points for several galaxies.

    Parameters
    ----------
    dirnames : list of str
        Galaxy directory names.
    vr_dict : dict
        Dictionary keyed by VFID or dirname with recession velocity.
        Example keys can be 'VFID5855' or full dirname.
    nsigma : float
        Detection threshold.
    alma : bool or dict
        If bool, use same alma flag for all galaxies.
        If dict, keyed by VFID or dirname.
    rootdir : str or None
        Parent directory containing galaxy directories.

    Returns
    -------
    ks_df : pandas.DataFrame
        One row per detected CO+Halpha pixel.
    """

    import numpy as np
    import pandas as pd

    rows = []

    for dirname in dirnames:
        vfid = dirname.split("-")[0]

        if dirname in vr_dict:
            vr = vr_dict[dirname]
        elif vfid in vr_dict:
            vr = vr_dict[vfid]
        else:
            raise KeyError(f"No vr found for {dirname} or {vfid}")

        if isinstance(alma, dict):
            if dirname in alma:
                alma_flag = alma[dirname]
            elif vfid in alma:
                alma_flag = alma[vfid]
            else:
                alma_flag = False
        else:
            alma_flag = alma

        result = get_ks_points_one_galaxy(
            dirname,
            vr=vr,
            nsigma=nsigma,
            alma=alma_flag,
            rootdir=rootdir,
            H0=H0,
            return_mask=False,
            verbose=verbose,
        )

        n = result["n_points"]

        if n == 0:
            continue

        rows.append(
            pd.DataFrame(
                {
                    "dirname": dirname,
                    "vfid": vfid,
                    "galtitle":result["galtitle"],
                    "sigma_H2_kpc2": result["sigma_H2_kpc2"],
                    "sigma_H2_pc2": result["sigma_H2_pc2"],
                    "sigma_sfr_kpc2": result["sigma_sfr_kpc2"],
                    "log_sigma_H2_pc2": result["log_sigma_H2_pc2"],
                    "log_sigma_sfr": result["log_sigma_sfr"],
                    "co_snr": result["co_snr"],
                    "ha_snr": result["ha_snr"],
                }
            )
        )

        # rows.append(
        #     pd.DataFrame(
        #         {
        #             "dirname": dirname,
        #             "vfid": vfid,
        #             "sigma_H2": result["sigma_H2"],
        #             "sigma_sfr": result["sigma_sfr"],
        #             "log_sigma_H2": result["log_sigma_H2"],
        #             "log_sigma_sfr": result["log_sigma_sfr"],
        #         }
        #     )
        # )

    if len(rows) == 0:
        return pd.DataFrame(
            columns=[
                "dirname",
                "vfid",
                "sigma_H2",
                "sigma_sfr",
                "log_sigma_H2",
                "log_sigma_sfr",
            ]
        )

    ks_df = pd.concat(rows, ignore_index=True)

    if verbose:
        print()
        print(f"Total KS points: {len(ks_df)}")
        print(f"Number of galaxies: {ks_df['vfid'].nunique()}")

    return ks_df

def _binned_percentiles_by_x(
    x,
    y,
    bins,
    min_points_per_bin=5,
):
    """
    Calculate median x, median y, and 16/84 percentile y values
    in bins of x.
    """

    import numpy as np

    x = np.asarray(x)
    y = np.asarray(y)

    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]

    x_med = []
    y_med = []
    y_p16 = []
    y_p84 = []
    n_bin = []

    for i in range(len(bins) - 1):
        if i == len(bins) - 2:
            in_bin = (x >= bins[i]) & (x <= bins[i + 1])
        else:
            in_bin = (x >= bins[i]) & (x < bins[i + 1])

        n = np.sum(in_bin)

        if n < min_points_per_bin:
            continue

        x_med.append(np.nanmedian(x[in_bin]))
        y_med.append(np.nanmedian(y[in_bin]))
        y_p16.append(np.nanpercentile(y[in_bin], 16))
        y_p84.append(np.nanpercentile(y[in_bin], 84))
        n_bin.append(n)

    return {
        "x_med": np.asarray(x_med),
        "y_med": np.asarray(y_med),
        "y_p16": np.asarray(y_p16),
        "y_p84": np.asarray(y_p84),
        "n_bin": np.asarray(n_bin),
    }


def plot_resolved_ks_relation(
    ks_df,
    ax=None,
    x_col="log_sigma_H2_pc2",
    y_col="log_sigma_sfr_kpc2",
    group_col="vfid",
    alpha=0.2,
    s=18,
    nbins=10,
    bins=None,
    bins_by="global",
    min_points_per_bin=8,
    show_points=True,
    show_running_median=True,
    show_percentile_region=True,
    median_lw=2.5,
    percentile_alpha=0.18,
    show_legend=True,
    plot_constant_tdep=True,
    tdep_yr=(1e8, 1e9, 1e10),
    plot_comparisons=False,
):
    """
    Plot resolved Kennicutt-Schmidt relation for multiple galaxies.

    For each galaxy:
      - plot individual points with alpha=0.2 by default
      - calculate running median in bins of log Sigma_H2
      - shade 16-84 percentile range in each bin

    Expected default columns:
      log_sigma_H2_pc2
      log_sigma_sfr_kpc2
      vfid
    """

    import numpy as np
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5.2))
    else:
        fig = ax.figure

    df = ks_df.copy()

    good = (
        np.isfinite(df[x_col])
        & np.isfinite(df[y_col])
    )
    df = df.loc[good]

    if len(df) == 0:
        raise ValueError("No finite KS points found.")

    ##################################################
    # Set up global bins
    ##################################################

    if bins is None and bins_by == "global":
        xmin = np.nanmin(df[x_col])
        xmax = np.nanmax(df[x_col])
        bins_global = np.linspace(xmin, xmax, nbins + 1)
    elif bins is not None:
        bins_global = np.asarray(bins)
    else:
        bins_global = None

    ##################################################
    # Stable matplotlib color cycle
    ##################################################

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    ##################################################
    # Plot each galaxy separately
    ##################################################

    for j, (gal, g) in enumerate(df.groupby(group_col)):
        x = np.asarray(g[x_col])
        y = np.asarray(g[y_col])

        color = color_cycle[j % len(color_cycle)]

        if show_points:
            ax.scatter(
                x,
                y,
                s=s,
                alpha=alpha,
                color=color,
                label=str(gal),
                edgecolors="none",
            )

        if not show_running_median:
            continue

        if bins_by == "galaxy" and bins is None:
            if len(g) < min_points_per_bin:
                continue

            xmin_g = np.nanmin(x)
            xmax_g = np.nanmax(x)

            if not np.isfinite(xmin_g) or not np.isfinite(xmax_g) or xmin_g == xmax_g:
                continue

            bins_use = np.linspace(xmin_g, xmax_g, nbins + 1)
        else:
            bins_use = bins_global

        stats = _binned_percentiles_by_x(
            x,
            y,
            bins_use,
            min_points_per_bin=min_points_per_bin,
        )

        x_med = stats["x_med"]
        y_med = stats["y_med"]
        y_p16 = stats["y_p16"]
        y_p84 = stats["y_p84"]

        if len(x_med) == 0:
            continue

        if show_percentile_region:
            ax.fill_between(
                x_med,
                y_p16,
                y_p84,
                color=color,
                alpha=percentile_alpha,
                linewidth=0,
                zorder=2,
            )

        ax.plot(
            x_med,
            y_med,
            color=color,
            lw=median_lw,
            zorder=3,
        )

    ##################################################
    # Constant depletion-time guide lines
    ##################################################

    if plot_constant_tdep:
        xgrid = np.linspace(
            np.nanmin(df[x_col]) - 0.2,
            np.nanmax(df[x_col]) + 0.2,
            100,
        )

        for tdep in tdep_yr:
            log_tdep = np.log10(tdep)

            # x is log Sigma_H2 in Msun/pc^2.
            # Convert Sigma_H2 to Msun/kpc^2 by adding 6 dex.
            # log Sigma_SFR = log Sigma_H2(pc^-2) + 6 - log tau_dep(yr)
            ygrid = xgrid + 6 - log_tdep

            if log_tdep == 9:
                thiscolor = 'k'
            else:
                thiscolor = "0.6"
            ax.plot(
                xgrid,
                ygrid,
                ls="--",
                lw=1,
                color=thiscolor,
                alpha=0.8,
                zorder=1,
            )

            if log_tdep == 9:            

                label = fr"COLIBRE / {tdep / 1e9:g} Gyr"
                thiscolor='k' 
                thisfontsize=11               
            else:
                label = fr"{tdep / 1e9:g} Gyr"
                thiscolor = "0.35"
                thisfontsize=9
            ax.text(
                xgrid[-1],
                ygrid[-1],
                label,
                fontsize=thisfontsize,
                color=thiscolor,
                ha="right",
                va="bottom",
            )

    ##################################################
    # Plot comparison data
    ##################################################
    if plot_comparisons:
        # Leroy+2013
        log10_a = -2.40 # pm .09
        n = 0.95 #pm 0.15
        scatter = 0.3

        xgrid = np.linspace(0,2.5,
            
            100,
        )
        
        ygrid = n * (log10_a + xgrid - np.log10(10))
        ax.plot(xgrid, ygrid, ls="-", lw='2',color='0.35',label=r"Local disks: H$\alpha$+24$\mu$m (Leroy+2013)")
        plt.fill_between(xgrid, y1=ygrid + scatter, y2=ygrid - scatter, color='k',alpha=.15)

        label = f"local disks"
        ax.text(
                xgrid[-1],
                ygrid[-1],
                label,
                fontsize=11,
                color="0.35",
                ha="left",
                va="bottom",
            )
    ##################################################
    # Labels and formatting
    ##################################################

    ax.set_xlabel(
        r"$\log_{10}\,\Sigma_{\rm H_2}\ [{\rm M_\odot\,pc^{-2}}]$",fontsize=16
    )
    ax.set_ylabel(
        r"$\log_{10}\,\Sigma_{\rm SFR_{H\alpha}}\ [{\rm M_\odot\,yr^{-1}\,kpc^{-2}}]$",fontsize=16
    )

    ax.grid(alpha=0.25)

    if show_legend:
        leg = ax.legend(
            fontsize=8,
            frameon=False,
            loc="best",
            markerscale=1.8,
            
        )
        for lh in leg.legend_handles[0:-1]: 
            lh.set_alpha(1)

    fig.tight_layout()

    return fig, ax




# def plot_resolved_ks_relation(
#     ks_df,
#     ax=None,
#     alpha=0.2,
#     s=18,
#     show_legend=True,
#     plot_constant_tdep=True,
#     plot_running_median=True,
#     nbins=12,
#     min_points_per_bin=10,
#     median_color="k",
#     median_lw=2.5,
#     percentile_alpha=0.25,
# ):
#     """
#     Plot resolved Kennicutt-Schmidt relation with optional running median
#     and 16-84 percentile range.

#     Expected dataframe columns:
#         log_sigma_H2_pc2
#         log_sigma_sfr
#         vfid

#     Units:
#         x-axis: log10 Sigma_H2 [Msun pc^-2]
#         y-axis: log10 Sigma_SFR [Msun yr^-1 kpc^-2]
#     """

#     import numpy as np
#     import matplotlib.pyplot as plt

#     if ax is None:
#         fig, ax = plt.subplots(figsize=(6, 5))
#     else:
#         fig = ax.figure

#     ##################################################
#     # Plot individual points
#     ##################################################

#     for vfid, g in ks_df.groupby("vfid"):
#         ax.scatter(
#             g["log_sigma_H2_pc2"],
#             g["log_sigma_sfr_kpc2"],
#             s=s,
#             alpha=alpha,
#             label=vfid,
#         )

#     ##################################################
#     # Running median and 16-84 percentile envelope
#     ##################################################

#     if plot_running_median:
#         x = np.asarray(ks_df["log_sigma_H2_pc2"])
#         y = np.asarray(ks_df["log_sigma_sfr_kpc2"])

#         good = np.isfinite(x) & np.isfinite(y)
#         x = x[good]
#         y = y[good]

#         bins = np.linspace(np.nanmin(x), np.nanmax(x), nbins + 1)
#         bin_centers = 0.5 * (bins[:-1] + bins[1:])

#         x_med = []
#         y_med = []
#         y_p16 = []
#         y_p84 = []
#         n_bin = []

#         for i in range(nbins):
#             in_bin = (x >= bins[i]) & (x < bins[i + 1])

#             # include the right edge in the last bin
#             if i == nbins - 1:
#                 in_bin = (x >= bins[i]) & (x <= bins[i + 1])

#             if np.sum(in_bin) < min_points_per_bin:
#                 continue

#             x_med.append(np.nanmedian(x[in_bin]))
#             y_med.append(np.nanmedian(y[in_bin]))
#             y_p16.append(np.nanpercentile(y[in_bin], 16))
#             y_p84.append(np.nanpercentile(y[in_bin], 84))
#             n_bin.append(np.sum(in_bin))

#         x_med = np.asarray(x_med)
#         y_med = np.asarray(y_med)
#         y_p16 = np.asarray(y_p16)
#         y_p84 = np.asarray(y_p84)

#         if len(x_med) > 0:
#             ax.fill_between(
#                 x_med,
#                 y_p16,
#                 y_p84,
#                 color=median_color,
#                 alpha=percentile_alpha,
#                 linewidth=0,
#                 label="16-84 percentile",
#             )

#             ax.plot(
#                 x_med,
#                 y_med,
#                 color=median_color,
#                 lw=median_lw,
#                 label="Running median",
#             )

#     ##################################################
#     # Constant depletion-time guide lines
#     ##################################################

#     if plot_constant_tdep:
#         xgrid = np.linspace(
#             np.nanmin(ks_df["log_sigma_H2_pc2"]) - 0.2,
#             np.nanmax(ks_df["log_sigma_H2_pc2"]) + 0.2,
#             100,
#         )

#         # Sigma_H2 is in Msun/pc^2.
#         # Convert to Msun/kpc^2 by adding 6 dex:
#         # log Sigma_SFR = log Sigma_H2(pc^-2) + 6 - log tau_dep(yr)
#         for log_tdep in [8, 9, 10]:
#             ygrid = xgrid + 6 - log_tdep

#             ax.plot(
#                 xgrid,
#                 ygrid,
#                 ls="--",
#                 lw=1,
#                 color="0.6",
#                 alpha=0.8,
#             )

#             ax.text(
#                 xgrid[-1],
#                 ygrid[-1],
#                 rf"$10^{{{log_tdep}}}$ yr",
#                 fontsize=9,
#                 color="0.35",
#                 ha="right",
#                 va="bottom",
#             )

#     ##################################################
#     # Labels and formatting
#     ##################################################

#     ax.set_xlabel(
#         r"$\log_{10}\,\Sigma_{\rm H_2}\ [{\rm M_\odot\,pc^{-2}}]$"
#     )
#     ax.set_ylabel(
#         r"$\log_{10}\,\Sigma_{\rm SFR}\ [{\rm M_\odot\,yr^{-1}\,kpc^{-2}}]$"
#     )

#     ax.grid(alpha=0.25)

#     if show_legend:
#         ax.legend(
#             fontsize=8,
#             frameon=False,
#             loc="best",
#             markerscale=1.5,
#         )

#     fig.tight_layout()

#     return fig, ax



def plot_resolved_ks_relation_kpc(
    ks_df,
    ax=None,
    alpha=0.2,
    s=20,
    show_legend=True,
    plot_constant_tdep=True,
):
    """
    Plot resolved Kennicutt-Schmidt relation for several galaxies.

    Expects a dataframe from get_ks_points_many_galaxies().
    """

    import numpy as np
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.figure

    for vfid, g in ks_df.groupby("vfid"):
        ax.scatter(
            g["log_sigma_H2"],
            g["log_sigma_sfr"],
            s=s,
            alpha=alpha,
            label=vfid,
        )

    if plot_constant_tdep:
        xgrid = np.linspace(
            np.nanmin(ks_df["log_sigma_H2"]) - 0.2,
            np.nanmax(ks_df["log_sigma_H2"]) + 0.2,
            100,
        )

        # Sigma_SFR = Sigma_H2 / tau_dep
        # Units:
        #   Sigma_H2: Msun/kpc^2
        #   Sigma_SFR: Msun/yr/kpc^2
        for log_tdep in [8, 9, 10]:
            ygrid = xgrid - log_tdep

            ax.plot(
                xgrid,
                ygrid,
                ls="--",
                lw=1,
                color="0.5",
                alpha=0.7,
            )

            ax.text(
                xgrid[-1],
                ygrid[-1],
                rf"$10^{{{log_tdep}}}$ yr",
                fontsize=9,
                color="0.35",
                ha="right",
                va="bottom",
            )

    ax.set_xlabel(r"$\log_{10}\,\Sigma_{\rm H_2}\ [{\rm M_\odot\,kpc^{-2}}]$")
    ax.set_ylabel(r"$\log_{10}\,\Sigma_{\rm SFR}\ [{\rm M_\odot\,yr^{-1}\,kpc^{-2}}]$")

    ax.grid(alpha=0.25)

    if show_legend:
        ax.legend(
            fontsize=8,
            frameon=False,
            loc="best",
            markerscale=1.5,
        )

    fig.tight_layout()

    return fig, ax



def get_ks_points_many_galaxies_binned(
    dirnames,
    vr_dict,
    nsigma=3,
    alma=False,
    rootdir=None,
    H0=74.,
    binsize_pix=None,
    binsize_kpc=None,
    min_valid_frac=0.5,
    verbose=True,
):
    """
    Collect beam-binned KS points for several galaxies.

    Returns one pandas DataFrame with one row per detected beam-sized region.
    """

    import pandas as pd

    rows = []

    for dirname in dirnames:
        t = dirname.split("-")
        vfid = t[0]
        
        if dirname in vr_dict:
            vr = vr_dict[dirname]
        elif vfid in vr_dict:
            vr = vr_dict[vfid]
        else:
            raise KeyError(f"No vr found for {dirname} or {vfid}")

        if isinstance(alma, dict):
            if dirname in alma:
                alma_flag = alma[dirname]
            elif vfid in alma:
                alma_flag = alma[vfid]
            else:
                alma_flag = False
        else:
            alma_flag = alma

        result = get_ks_points_one_galaxy_binned(
            dirname,
            vr=vr,
            nsigma=nsigma,
            alma=alma_flag,
            rootdir=rootdir,
            H0=H0,
            binsize_pix=binsize_pix,
            binsize_kpc=binsize_kpc,
            min_valid_frac=min_valid_frac,
            verbose=verbose,
            )


        if result["n_points"] == 0:
            continue

        # reference dictionary
        # result = {
        #     "galaxy": dirname,
        #     "vfid": vfid,

        #     "sigma_H2_kpc2": sigma_H2_out,
        #     "sigma_H2_pc2": sigma_H2_pc2_out,
        #     "sigma_sfr_kpc2": sigma_sfr_out,

        #     "log_sigma_H2_pc2": np.log10(sigma_H2_pc2_out),
        #     "log_sigma_sfr": np.log10(sigma_sfr_out),

        #     "co_snr": co_snr_bin[detected_both_bin],
        #     "ha_snr": ha_snr_bin[detected_both_bin],
        #     "xbin": xbin[detected_both_bin],
        #     "ybin": ybin[detected_both_bin],
        #     "n_valid_pix": n_valid_bin[detected_both_bin],
        #     "n_points": len(sigma_H2_out),
        #     "binsize_pix": binsize_pix,
        #     "bin_npix": bin_npix,
        #     "beam_area_pix": beam_area_pix,
        #     "min_valid_pixels": min_valid_pixels,
        #     "kpc_per_pixel": kpc_per_pixel,
        #     }
            

        rows.append(
            pd.DataFrame(
                {
                    "dirname": dirname,
                    "vfid": vfid,
                    "galtitle": result["galtitle"],                    
                    "sigma_H2_pc2": result["sigma_H2_pc2"],
                    "sigma_sfr_kpc2": result["sigma_sfr_kpc2"],
                    "log_sigma_H2_pc2": result["log_sigma_H2_pc2"],
                    "log_sigma_sfr_kpc2": result["log_sigma_sfr"],
                    "co_snr": result["co_snr"],
                    "ha_snr": result["ha_snr"],
                    "xbin": result["xbin"],
                    "ybin": result["ybin"],
                    "n_valid_pix": result["n_valid_pix"],
                    "binsize_pix": result["binsize_pix"],
                    "beam_area_pix": result["beam_area_pix"],
                }
            )
        )

    if len(rows) == 0:
        return pd.DataFrame(
            columns=[
                "dirname",
                "vfid",
                "sigma_H2_pc2",
                "sigma_sfr_kpc2",
                "log_sigma_H2_pc2",
                "log_sigma_sfr_kpc2",
                "co_snr",
                "ha_snr",
                "xbin",
                "ybin",
                "n_valid_pix",
                "binsize_pix",
                "beam_area_pix",
            ]
        )

    ks_df = pd.concat(rows, ignore_index=True)

    if verbose:
        print()
        print(f"Total beam-binned KS points: {len(ks_df)}")
        print(f"Number of galaxies: {ks_df['vfid'].nunique()}")

    return ks_df



def get_depletion_map(dirname,vr=None, H0=74., cmap='magma_r', verbose=False, sfr_limit=1e-5, alma=False,convolveflag=False,alt_title=None,
                          show_tdep_limits=True):
    """
    GOAL:
    * get depletion map for CO

    PARAMS:
    * dirname - root name on halpha cutout directory
    * vr - this should actually be Vcosmis
    * H0 - optional input - hubble constant.  use 74 in Castignani+2022


    assume the haimage is linear, in units of SFR/pix/1.e3

    Note in alma_sources.ipynb
    """
    from astropy.wcs import WCS
    vfid = dirname.split('-')[0]
    print("VFID = ",vfid)

    #%matplotlib inline
    os.chdir(homedir+'/research/Virgo-dev/cont-sub-gr')
    # add scale factor for continuue after the directory name

    cwd = os.getcwd()
    os.chdir(dirname)

    # get D_Mpc
    
    D_Mpc = vr / H0

    cblabels = [ r'$\rm Jy~beam^{-1}~ km~s^{-1}$',\
                    r'$\rm SFR(M_\star/yr)/kpc^2$',\
                    r'$\rm \log_{10}(\tau_{dep}/yr)$',]

    titles = ['$CO~ Intensity$',r'$\Sigma_{SFR}$',r'$Depletion \ time$']
    temp = dirname.split('-')
    galname = temp[0]+'-'+temp[1]
    if alt_title is not None:
        galname = alt_title
    ##################################################
    ## CO Image
    ##################################################
    COfilename = CO_file[vfid]
    # need to get the header from Francoise's version
    cdat_unmasked, cheader = convert_alma3d_alma2d(COfilename)

    # print CO header
    #print(cheader)
    
    co_imwcs = WCS(cheader)
    COheader = cheader
    
    COmask_filename = CO_mask[vfid]
    if COmask_filename is None:
        print("did not find masked CO file!")        
        cdat_masked = cdat_unmasked
    else:
        print("reading masked CO file!")
        cdat_masked = fits.getdata(COmask_filename)

    #co_mdat = np.ma.array(cdat,mask=maskdat)

    # cdat is in Jy/beam * km/s
    # convert CO flux to H2 mass based on eqn 1 from Castignani+2022
    alpha_CO = 4.3
    #cdat = cdat_unmasked * 1.05e4 * alpha_CO/4.3 * (D_Mpc)**2 # in Msun/beam
    cdat = cdat_masked * 1.05e4 * alpha_CO/4.3 * (D_Mpc)**2 # in Msun/beam

    # add correction to convert from Jy/beam to Jy/pixel
    pixscale_deg_out = wcs.utils.proj_plane_pixel_scales(WCS(cheader))[0] # pixel scale
    pixscale_area_sr_out = (pixscale_deg_out * np.pi/180)**2 # area per pixel
    pixel_area_sq_deg = pixscale_deg_out**2 # pixel area in sq deg
    
    # get beam area in sq deg
    beam_area_sq_deg = np.pi/4 * float(cheader['BMAJ']) * float(cheader['BMIN'])

    # scale factor between beam area / pixel area
    scale_factor = beam_area_sq_deg/pixel_area_sq_deg # beam area/pixel area

    cdat = cdat / scale_factor # this should now be Msun/(pixel area)

    # convert units for the CO rms value
    if alma: # use rms image for alma
        rms_map = fits.getdata(co_rms_map[vfid]) # rms per beam
        co_rms_Msun_beam = rms_map * 1.05e4 * alpha_CO/4.3 * (D_Mpc)**2 # Msun/beam
    else: # use the fixed rms value for noema
        co_rms_Msun_beam =co_rms[vfid] * 1.05e4 * alpha_CO/4.3 * (D_Mpc)**2

    # GL thinks we should scale by sqrt of beam area/pixel area
    # signal scales with scale_factor, but rms scales as sqrt(scale_factor)
    # so snr scales as sqrt(scale_factor)
    co_rms_Msun_pixel = (co_rms_Msun_beam ) / np.sqrt(scale_factor) # Msun/pixel
    #print("DEBUG: ",co_rms_Msun_pixel)
    try:
        print(f"CO rms Msun/pixel = {co_rms_Msun_pixel:.3e}")
    except TypeError:
        print(f"CO rms Msun/pixel = {np.nanmedian(co_rms_Msun_pixel):.3e}")
    if vr is not None:
        # scale data to convert to flux/kpc^2
        kpc_per_pixel = get_kpc_per_pixel(co_imwcs,vr)
        try:
            print(f"CO rms Msun/kpc^2 = {co_rms_Msun_pixel/kpc_per_pixel**2:.3e}, kpc/pixel={kpc_per_pixel:.2e}, vr={vr:.1f}")
        except TypeError:
            print(f"CO rms Msun/kpc^2 = {np.nanmedian(co_rms_Msun_pixel)/kpc_per_pixel**2:.3e}, kpc/pixel={kpc_per_pixel:.2e}, vr={vr:.1f}")

        print()
        
    ##################################################
    ## SFR Image
    ##################################################
    sfrim = dirname+"-sfr-vr.fits"
    hamask_filename = dirname+'-R-mask.fits'
    # sfr maps are scaled by a factor of 1e3 for convenience
    # divide sfr by 1e3 to get true values.
    hdu = fits.open(sfrim)[0]
    #hdu.data = hdu.data
    sfr_header = hdu.header

    try:
        ha_fwhm = sfr_header['SEFWHM']
    except KeyError:
        print("WARNING: no SEFWHM is header, assuming 1.5 arcsec")
        ha_fwhm = 1.5
    
    # convert to surface brightness
    # times by area in sr
    #pixscale_deg_in = np.abs(hdu.header['CDELT1'])
    pixscale_deg_in = wcs.utils.proj_plane_pixel_scales(WCS(hdu.header))[0]
    ha_imwcs = WCS(hdu.header)
    
    wcs_in = WCS(hdu.header)
    hdu.header = wcs_in.to_header()

    
    #print(WCS(hdu.header))
    print(f"pix scale in input image = {pixscale_deg_in*3600:.2f}")
    pixscale_area_sr_in = (pixscale_deg_in * np.pi/180)**2 # pixel area in steradians

    ###
    # CONVERT TO SB UNITS - DON'T NEED TO ANYMORE WITH reproject_adaptive, conserve_flux=True
    # input data is in flux/pixel -> sb, or flux/sr
    ###
    #hdu.data = hdu.data #/ pixscale_area_sr_in
    
    # replace nans with zero
    #hdu.data[hdu.data == np.nan] = 0
    
    hdat = hdu.data
    #hdat = hdu.data/1.e3
    hheader = hdu.header

    hmaskdat, hmaskheader = fits.getdata(hamask_filename, header=True)
    ha_mdat = np.ma.array(hdat, mask=hmaskdat)

    if verbose:
        print("CHECK")
        print(f"{vfid}: max value in image before reprojection = {np.nanmax(hdat)}, halpha_rms = {halpha_rms[vfid]}")


    # ##################################################
    # ## Convolve Halpha Image to match CO
    # ##################################################

    # if convolveflag:
    #     # convolve halpha to radio
    #     sys.path.append(homedir+'/github/havirgo/python/')
    #     import convolve_halpha_2_radio as cr
    #     radiobeam_PA = cheader['BPA']# in degrees
    #     radiobeam_a = (cheader['BMAJ']*u.deg).to(u.arcsec).value
    #     radiobeam_b = (cheader['BMIN']*u.deg).to(u.arcsec).value
        
    #     outname, csfr_dat = cr.convolve_image(sfrim, ha_fwhm, radiobeam_a, radiobeam_b, radiobeam_PA)


    ##################################################
    ## Convolve/Reproject Halpha Image to match CO grid
    ##################################################
    from reproject import reproject_adaptive, reproject_interp
    #from reproject import reproject_adaptive

    out_header = cheader

    if convolveflag:
        outfile = f"{vfid}_sfr_reproject_adaptive_True_convolved_alma.fits"
    else:
        outfile = f"{vfid}_sfr_reproject_adaptive_True_alma.fits"

    if os.path.exists(outfile):
        print(f"reading existing reprojected SFR image: {outfile}")

        if convolveflag:
            crsfr_dat, cached_header = fits.getdata(outfile, header=True)
        else:
            rsfr_dat, cached_header = fits.getdata(outfile, header=True)

    else:
        print(f"creating reprojected SFR image: {outfile}")

        if convolveflag:
            # convolve halpha to radio
            sys.path.append(homedir + '/github/havirgo/python/')
            import convolve_halpha_2_radio as cr

            radiobeam_PA = cheader['BPA']  # degrees
            radiobeam_a = (cheader['BMAJ'] * u.deg).to(u.arcsec).value
            radiobeam_b = (cheader['BMIN'] * u.deg).to(u.arcsec).value

            outname, csfr_dat = cr.convolve_image(
                sfrim,
                ha_fwhm,
                radiobeam_a,
                radiobeam_b,
                radiobeam_PA,
            )

            convolved_hdu = fits.ImageHDU(data=csfr_dat, header=hdu.header)

            crsfr_dat, crfootprint = reproject_adaptive(
                convolved_hdu,
                out_header,
                conserve_flux=True,
            )

            fits.writeto(outfile, crsfr_dat, header=out_header, overwrite=True)

        else:
            rsfr_dat, rfootprint = reproject_adaptive(
                hdu,
                out_header,
                conserve_flux=True,
            )

            fits.writeto(outfile, rsfr_dat, header=out_header, overwrite=True)


        

        
    # ##################################################
    # ## Reproject Image
    # ##################################################
    

    # ##
    # # reproject the halpha image onto the CO image
    # ##    
    # # TODONE - are we doing the reprojection right - yes, adaptive with conserve_flux=True is correct
    # out_header = cheader

    # if convolveflag:
    #     convolved_hdu = fits.ImageHDU(data=csfr_dat,header=hdu.header)
    #     crsfr_dat, crfootprint = reproject_adaptive(convolved_hdu, out_header, conserve_flux=True)
    # else:
    #     rsfr_dat, rfootprint = reproject_adaptive(hdu, out_header, conserve_flux=True)


        
    # # when using reproject_exact, need to be in sb units
    # #rsfr_sb, rfootprint = reproject_exact(hdu, out_header)#, conserve_flux=True)        
    # # need to convert back to flux
    # #rsfr_dat = rsfr_sb #* (pixscale_area_sr_out)#/pixscale_area_sr_in)
    
    # #print("ratio of out/in pixel areas = ",(pixscale_area_sr_out/pixscale_area_sr_in))

    # #print(f"ratio of rsfr_dat / rsfr_sb = {np.nanmean(rsfr_dat/rsfr_sb)}")
    
    # # write out reprojected halpha image
    # if convolveflag:
    #     outfile = f"{vfid}_sfr_reproject_adaptive_True_convolved_alma.fits"
    #     fits.writeto(outfile,crsfr_dat, header = out_header, overwrite=True)
    # else:
    #     outfile = f"{vfid}_sfr_reproject_adaptive_True_alma.fits"
    #     fits.writeto(outfile,rsfr_dat, header = out_header, overwrite=True)
    

    # Reproject the Halpha/R mask onto the CO grid.
    # Use nearest-neighbor because this is a mask, not flux.
    hmask_hdu = fits.ImageHDU(
        data=np.asarray(hmaskdat, dtype=float),
        header=hdu.header,
    )

    rhmask_dat, rhmask_footprint = reproject_interp(
        hmask_hdu,
        out_header,
        order="nearest-neighbor",
    )

    # True where the Halpha/SFR pixel should be ignored.
    ha_mask_on_co_grid = (
        (rhmask_dat > 0.5)
        | ~np.isfinite(rhmask_dat)
        | (rhmask_footprint == 0)
    )

    
        
    # compare the sum of the original and reprojected images
    #flux_in = np.nansum(hdu.data) * pixscale_area_sr_in

    #flux_out = np.nansum(rsfr_dat) * pixscale_area_sr_out
    #print("Output total flux:", flux_out)
    #print("Flux ratio (out/in):", flux_out/flux_in)

    #print(f"sum of original = {np.sum(hdat):.3e}, reprojected = {np.sum(rsfr_dat):.3e}")
    #print()

        
    ##################################################
    ## Get Depletion Map
    ##################################################

    #if vr is not None:
        # scale data to convert to flux/kpc^2
    #    kpc_per_pixel = get_kpc_per_pixel(co_imwcs,vr)
    #    rsfr_dat = rsfr_dat/(kpc_per_pixel**2)



    # get zoom in region
    ymax,xmax = cdat.shape
    xcenter = xmax/2
    ycenter = ymax/2
    xzoom, yzoom = codepletion_zoom[vfid]

    x1 = int(xcenter - (xmax/2)*xzoom)
    x2 = int(xcenter + (xmax/2)*xzoom)

    y1 = int(ycenter - (ymax/2)*yzoom)
    y2 = int(ycenter + (ymax/2)*yzoom)
    #print(f"Testing: xmax,ymax={xmax},{ymax}, x1,x2,y1,y2={x1,x2,y1,y2}")


    ################################################
    ## CALCULATE EXTINCTION MAP WITH UPPER LIMITS
    ################################################

    # one magnitude of extinction at halpha
    extinction_correction = 2.5 # assume one magnitude until we have something better
    #extinction_correction = 10 # assume one magnitude until we have something better
    if convolveflag:
        sfr_dep_dat = crsfr_dat
    else:
        sfr_dep_dat = rsfr_dat
    depletion = cdat / (sfr_dep_dat * extinction_correction)

    NSIGMA=3 # noise cut
    # cdat is nan but halpha is > 3*image_std
    #halpha_nan_flag = sfr_dep_dat < 3*halpha_rms[vfid] # halpha_rms measured from sfr image
    halpha_nan_flag = sfr_dep_dat < NSIGMA*halpha_rms[vfid] # halpha_rms measured from sfr image

    print(f"SFR rms Msun/yr/pixel = {halpha_rms[vfid]:.3e}")
    if vr is not None:
        # scale data to convert to flux/kpc^2
        sfr_wcs = WCS(sfr_header)
        kpc_per_pixel = get_kpc_per_pixel(sfr_wcs,vr)
        print(f"SFR rms Msun/yr/kpc^2 = {halpha_rms[vfid]/kpc_per_pixel**2:.3e}")
        print()

    
    #if verbose:
    #    print("CHECK")
    #    print(f"{vfid}: max value in image = {np.nanmax(rsfr_dat)}, max value in image before reprojection = {np.nanmax(hdat)}, halpha_rms = {halpha_rms[vfid]}")


    # don't need this b/c Gianluca redid the images with noise clipping
    #if vfid in ['VFID3757', 'VFID2822', 'VFID2140']: # noema galaxies
    #    # replace values in CO image with a nan if the value is < 2* rms
    #    co_lowflux = cdat < 2.5*co_rms[vfid]
    #    cdat[co_lowflux] = np.nan
    #    cdat_masked[co_lowflux] = np.nan        
    
    # Feb 13, 2026
    # adding case for noema to disregard regions where co flux is < 3sigma
    #co_nan_flag = ~(cdat == cdat) | (cdat < 1.*co_rms[vfid]) | (cdat == 0)
    #co_rms_Msun_pixel
    co_nan_flag = ~(cdat == cdat) | (cdat < NSIGMA*co_rms_Msun_pixel) | (cdat == 0)

    # testing a constant sfr map - should lead to continuous depletion time map WRT upper/lower limits
    #sfr_med = 3 * halpha_rms[vfid] * extinction_correction
    #sfr_dep_dat = sfr_med*np.ones_like(sfr_dep_dat)

    
    #plt.figure()
    #plt.imshow(sfr_dep_dat)
    #plt.colorbar()
    #plt.title("TEST SFR Image")
    
    #if verbose:
    #    plt.figure()
    #    plt.imshow(halpha_nan_flag, cmap='viridis', origin="lower")
    #    plt.colorbar()
    #    plt.title("halpha_nan_flag")
        
    #    plt.figure()
    #    plt.imshow(co_nan_flag, cmap='viridis', origin="lower")
    #    plt.colorbar()
    #    plt.title("co_nan_flag")
    
    # if CO is detected (not nan), but SFR isn't
    # 
    tdep_lower_flag = ~co_nan_flag & halpha_nan_flag

    # CO not detected (nan_flag), but SFR is (not nan)
    # CO_rms / SFR gives the upper limit b/c CO could be lower
    tdep_upper_flag = co_nan_flag & ~halpha_nan_flag


    # Do not identify upper/lower limits inside masked Halpha regions
    # such as stars or other bad pixels.
    ha_valid_on_co_grid = ~ha_mask_on_co_grid
    
    # CO detected, Halpha not detected, and Halpha pixel is usable
    tdep_lower_flag = ~co_nan_flag & halpha_nan_flag & ha_valid_on_co_grid

    # CO not detected, Halpha detected, and Halpha pixel is usable
    tdep_upper_flag = co_nan_flag & ~halpha_nan_flag & ha_valid_on_co_grid
    
    # create an image that just has the upper and lower limts so we can visualize
    tdep_limits_image = np.zeros_like(depletion)

    # set upper limits to CO upper limit / SFR
    # TODO: make sure halpha_rms is actually in SFR units
    if alma:
        tdep_limits_image[tdep_upper_flag] = (NSIGMA * co_rms_Msun_pixel[tdep_upper_flag]) / (sfr_dep_dat[tdep_upper_flag] * extinction_correction)
        # we detect CO but not SFR
        tdep_limits_image[tdep_lower_flag] = cdat[tdep_lower_flag]/ (NSIGMA * halpha_rms[vfid] * extinction_correction)    

    else:
        tdep_limits_image[tdep_upper_flag] = (NSIGMA * co_rms_Msun_pixel) / (sfr_dep_dat[tdep_upper_flag] * extinction_correction)

        upper_limit = tdep_limits_image[tdep_upper_flag]

        if upper_limit.size:
            print(f"tau_dep upper limits: min={upper_limit.min():.2e}, max={upper_limit.max():.2e}, med={np.nanmedian(upper_limit):.2e}")
        else:
            print("tau_dep upper limits: no upper-limit pixels")        
        #tdep_limits_image[tdep_upper_flag] = (co_rms_Msun_pixel) / (sfr_dep_dat[tdep_upper_flag] * extinction_correction)        
        ## we detect CO but not SFR
        #tdep_limits_image[tdep_lower_flag] = cdat[tdep_lower_flag]/ (3 * halpha_rms[vfid] * extinction_correction)    
        # print out min of lower limits
        #lower_limit = (cdat[tdep_lower_flag]/ (3 * halpha_rms[vfid] * extinction_correction)).flatten()
        #print(f"tau_dep lower limits: min={np.min(lower_limit)}, max={np.min(lower_limit)}")

        tdep_limits_image[tdep_lower_flag] = cdat[tdep_lower_flag] / (NSIGMA * halpha_rms[vfid] * extinction_correction)

        lower_limit = tdep_limits_image[tdep_lower_flag]

        if lower_limit.size:
            print(f"tau_dep lower limits: min={lower_limit.min():.2e}, max={lower_limit.max():.2e}, med={np.nanmedian(lower_limit):.2e}")
        else:
            print("tau_dep lower limits: no lower-limit pixels")
    
    # combine limits with depletion map
    limit_flag = tdep_lower_flag | tdep_upper_flag
    depletion_with_limits = np.zeros_like(depletion)
    depletion_with_limits[~limit_flag] = depletion[~limit_flag]

    depletion_with_limits[limit_flag] = tdep_limits_image[limit_flag]

    # save depletion map as a fits image

    dep_map_fits = fits.PrimaryHDU(depletion_with_limits)
    dep_map_fits.header = cheader
    dep_map_fits.header["BUNIT"] = ("yr", "Gas depletion time")
    dep_map_filename = sfrim.replace('-sfr-vr.fits','-tdepletion.fits')
    dep_map_fits.writeto(dep_map_filename, overwrite=True)
    print("saving depletion map to ",dep_map_filename)


    log_dep = np.log10(depletion_with_limits)
    dep_map_fits.header["BUNIT"] = ("log10(yr)", "Gas depletion time")    
    log_dep_map_filename = sfrim.replace('-sfr-vr.fits','-tdepletion-log.fits')
    dep_map_fits.writeto(log_dep_map_filename, overwrite=True)
    
    # testing if we remove limits


    # create smoothed version of limits for plotting contours
    # Use approximately one CO beam area as the minimum contourable region.
    # scale_factor = beam_area / pixel_area was already calculated above.
    limit_min_pixels = max(10, int(np.ceil(scale_factor)))

    tdep_upper_flag_contour = clean_limit_mask_for_contours(
        tdep_upper_flag,
        min_pixels=limit_min_pixels,
        do_closing=True,
    )

    tdep_lower_flag_contour = clean_limit_mask_for_contours(
        tdep_lower_flag,
        min_pixels=limit_min_pixels,
        do_closing=True,
    )
    
    #depletion_with_limits = cdat / (sfr_dep_dat * extinction_correction)
    if verbose:
        ####################################
        # plot image of upper/lower limits 
        ####################################   
        plt.figure(figsize=(12,4))
    
        plt.subplot(1,3,1)
        plt.imshow(tdep_lower_flag[y1:y2,x1:x2], cmap='viridis', origin="lower")#, vmin=0, vmax=1.5)
        cb = plt.colorbar()
        plt.title("lower limit flag")
        plt.subplot(1,3,2)
        plt.imshow(tdep_upper_flag[y1:y2,x1:x2], cmap='viridis', origin="lower")#, vmin=0, vmax=1.5)
        cb = plt.colorbar() 
        plt.title("upper limit flag")   
        plt.subplot(1,3,3)
        plt.imshow(np.log10(tdep_limits_image[y1:y2,x1:x2]), cmap='viridis', origin="lower")
        cb = plt.colorbar()
        plt.title("upper/lower limit image")    
        


    #################################################
    ## FIGURE - CO, SIGMA_SFR, DEPLETION MAP
    #################################################    

    # set up figure parameters
    myfigsize=codepletion_figsize[vfid]
    plt.figure(figsize=(myfigsize[0],myfigsize[1]))
    plt.subplots_adjust(wspace=0.05)

    size_ratio = 486/ymax
    size_ratio = 1
    #print(f"Testing: size_ratio={size_ratio:.2f}", codepletion_figsize['VFID5842'][1],myfigsize[1])


    # Limit contour styling
    upper_limit_contour_color = "magenta"  # Halpha detected, CO not detected
    lower_limit_contour_color = "deepskyblue" #"cyan"     # CO detected, Halpha not detected
    limit_contour_lw = 1.2
    limit_contour_alpha = 0.9

    #################################################
    ## PANEL 1 - CO
    #################################################    
    
    ax1 = plt.subplot(1,3,1, projection=co_imwcs)
    plt.imshow(cdat_masked[y1:y2,x1:x2], cmap='viridis')

    # save image of CO mass/kpc^2
    kpc_per_pixel = get_kpc_per_pixel(co_imwcs,vr)

    sigma_H2_kpc2 = cdat / kpc_per_pixel**2

    co_kpc_fits = fits.PrimaryHDU(sigma_H2_kpc2)
    co_kpc_fits.header = cheader
    co_kpc_fits.header["BUNIT"] = ("Msun/kpc2", "Molecular gas surface density")
    co_kpc_filename = sfrim.replace('-sfr-vr.fits', '-Sigma-H2-kpc2.fits')
    co_kpc_fits.writeto(co_kpc_filename, overwrite=True)


    # co_kpc_fits = fits.PrimaryHDU(cdat_masked/kpc_per_pixel**2)
    # co_kpc_fits.header = cheader
    # co_kpc_filename = sfrim.replace('-sfr-vr.fits','-Sigma-H2-kpc2.fits')
    # co_kpc_fits.writeto(co_kpc_filename, overwrite=True)
    print("saving Sigma H2 mass/kpc^2 ",co_kpc_filename)

    ###############################################
    # add a physical scale to the CO panel
    ###############################################
    pscale = wcs.utils.proj_plane_pixel_scales(co_imwcs)  # deg/pixel
    pscale_arcsec = np.abs(pscale[0]) * 3600              # arcsec/pixel


    # Contours showing upper limits:
    # Halpha detected, CO not detected
    if show_tdep_limits:
        show_both_limit_contours=False
        add_limit_contours(
            ax1,
            tdep_upper_flag_contour,
            tdep_lower_flag_contour,
            x1,
            x2,
            y1,
            y2,
            show_upper=True,
            show_lower=show_both_limit_contours,
            )


 
    
    plt.title(titles[0],fontsize=14)
    plt.text(-.1,0.5,galname,verticalalignment='center',transform=plt.gca().transAxes,fontsize=14,rotation=90)
    #plt.title("CO Intensity Map")
    #plt.colorbar(fraction=.08)
    #cbaxes = inset_axes(ax1, width="80%", height="5%", loc=8)
    # Create inset axes for colorbar below the x-axis
    cbaxes = inset_axes(
        ax1, 
        width="80%", # 80% of parent axes width
        height="5%", #f"{10*size_ratio:.0f}%", # 5% of axes height
        loc='lower center', 
        #bbox_to_anchor=(.1,-0.07*size_ratio,0.8,0.4*size_ratio), # shift below axes
        bbox_transform=ax1.transAxes
        )

    cb = plt.colorbar(cax=cbaxes, orientation='horizontal')            
    cb.set_label(label=cblabels[0],fontsize=12)
    hide_xyticks_wcs(ax1)

    # plot CO beam
    plot_HI_beam(ax1,COfilename,cheader,color='steelblue')


    add_scale(
        ax1,
        vr=vr,
        pscale=pscale_arcsec,
        color='k',
        barsize=1,
        fontsize=12,
        #xscale=.85,
        yscale=0.95,
    )    
    

    #################################################
    ## PANEL 2 - SIGMA_SFR
    #################################################    

    # plot 2
    ax2 = plt.subplot(1,3,2, projection=co_imwcs)

    # convert from SFR/pixel to SFR/kpc^2
    if vr is not None:
        # scale data to convert to flux/kpc^2
        kpc_per_pixel = get_kpc_per_pixel(co_imwcs,vr)
        
        if convolveflag:
            rsfr_dat_kpc2 = crsfr_dat/(kpc_per_pixel**2)
        else:
            rsfr_dat_kpc2 = rsfr_dat/(kpc_per_pixel**2)
        # save SFR/kpc^2 image
        # save image of CO mass/kpc^2
        # create a masked version
        rsfr_dat_kpc2_masked = rsfr_dat_kpc2.copy()
        rsfr_dat_kpc2_masked[ha_mask_on_co_grid] = np.nan
        sfr_kpc_fits = fits.PrimaryHDU(rsfr_dat_kpc2_masked)
        sfr_kpc_fits.header = cheader
        sfr_kpc_filename = sfrim.replace('-sfr-vr.fits','-Sigma-sfr-kpc2.fits')
        sfr_kpc_fits.writeto(sfr_kpc_filename, overwrite=True)
        print("saving sfr/kpc^2 ",sfr_kpc_filename)

    
        v1 = 0 # vmin[i]/kpc_per_pixel**2
        v2 = sfr_limit/kpc_per_pixel**2
        #myclevels = myclevels/kpc_per_pixel**2
        plt.imshow(rsfr_dat_kpc2[y1:y2,x1:x2], vmin=v1, vmax=v2, cmap=cmap)
    else:
        # plot result
        plt.imshow(rsfr_dat[y1:y2,x1:x2], vmin=0.0,vmax=1.e-5, cmap=cmap)

    # Contours showing lower limits:
    # CO detected, Halpha not detected

    if show_tdep_limits:
        add_limit_contours(
            ax2,
            tdep_upper_flag_contour,
            tdep_lower_flag_contour,
            x1,
            x2,
            y1,
            y2,
            show_upper=show_both_limit_contours,
            show_lower=True,
            )

    # lower_limit_mask_zoom = tdep_lower_flag_contour[y1:y2, x1:x2].astype(float)


    # if np.any(lower_limit_mask_zoom):
    #     ax2.contour(
    #         lower_limit_mask_zoom,
    #         levels=[0.5],
    #         colors=lower_limit_contour_color,
    #         linewidths=limit_contour_lw,
    #         alpha=limit_contour_alpha,
    #         )

        
    plt.title(titles[1],fontsize=14)
    #plt.title("SFR")
    #cbaxes = inset_axes(ax2, width="80%", height="5%", loc=8)
    # Create inset axes for colorbar below the x-axis
    cbaxes = inset_axes(
        ax2, 
        width="80%", # 80% of parent axes width
        height="5%", # 5% of axes height
        loc='lower center', 
        #bbox_to_anchor=(.1,-0.07*size_ratio,0.8,0.4*size_ratio), # shift below axes
        bbox_transform=ax2.transAxes
        )


    cb = plt.colorbar(cax=cbaxes, orientation='horizontal')            
    cb.set_label(label=cblabels[1],fontsize=12)                
    #plt.colorbar(fraction=.08)

    hide_xyticks_wcs(ax2)

    #print(hheader)
    # plot Halpha beam


    # plot halpha fwhm
    #if not convolveflag:
    plot_HI_beam(ax2,COfilename,cheader,color='steelblue', halpha=True, fwhm=ha_fwhm)
    
    #################################################
    ## PANEL 3- DEPLETION MAP
    #################################################    

    ax3 = plt.subplot(1,3,3, projection=co_imwcs)
    #plt.title("Depletion Map")
    v1, v2 = codepletion_levels[vfid]

    # adding upper/lower limits to depletion map
    #plt.imshow(np.log10(depletion[y1:y2,x1:x2]), cmap='viridis', vmin=v1, vmax=v2)
    plt.imshow(np.log10(depletion_with_limits[y1:y2,x1:x2]), cmap='viridis', vmin=v1, vmax=v2)    
    plt.title(titles[2],fontsize=14)

    #cbaxes = inset_axes(ax3, width="80%", height="5%", loc=8)
    # Create inset axes for colorbar below the x-axis
    cbaxes = inset_axes(
        ax3, 
        width="80%", # 80% of parent axes width
        height="5%", # 5% of axes height
        loc='lower center', 
        #bbox_to_anchor=(.1,-0.07*size_ratio**2,0.8,0.4*size_ratio), # shift below axes
        bbox_transform=ax3.transAxes
        )
    
    
    cb = plt.colorbar(cax=cbaxes, orientation='horizontal')            
    cb.set_label(label=cblabels[2],fontsize=12)        
    #plt.colorbar(fraction=.08, label='log10(t_depletion)')

    hide_xyticks_wcs(ax3)
    if show_tdep_limits:
        outname = os.path.join(plotdir,galname)+'-depletion-map-with-limits.png'
    else:
        outname = os.path.join(plotdir,galname)+'-depletion-map-no-limits.png'
        
    plt.savefig(outname,dpi=300,bbox_inches="tight")
        
    #plt.savefig(os.path.join(plotdir,galname)+'-depletion-map.png',dpi=150,bbox_inches="tight")        
    #plt.savefig(os.path.join(plotdir,dirname)+'-depletion-map.pdf',dpi=150,bbox_inches="tight")    

    
def plot_alma_tdepletion(dirname,xmin=None,xmax=None,ymin=None,ymax=None,xticks=True,figsize=[16,4],\
                            cbfrac=.08,cbaspect=20,clevels=[4],contourFlag=True,rmax=None,vr=None,\
                            logMstar=None,cmap='magma_r',markGroupCenter=False,COcolor='white',COlevels=None):
    """
   

    1 - legacy
    2 - CO image
    3 - SFR
    4 - CO/SFR

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
    nrow = 1
    allax = []
    nplots = [2,3,4,6,7,8]
    #nplots = [6,3,4]

    # legacy, integrated CO, sigma SFR, then depletion map
    indices = [2,3,4]

    # 2 = CO
    # 3 = SFR
    # 4 = depletion time
    
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
        #ax2=plt.subplot(nrow,ncol,j,projection=imwcs)
        ax2=plt.subplot(nrow,ncol,i+1,projection=imwcs)
        

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

    #ax5 = plt.subplot(nrow,ncol,5,projection=legwcs)
    #plt.imshow(jpeg_data)    

    leg_axes = [ax1,ax5]
    leg_axes = [ax1]    
    
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

    
    plt.savefig(os.path.join(plotdir,dirname)+'-CO-tdepletion.png',dpi=150,bbox_inches="tight")
    plt.savefig(os.path.join(plotdir,dirname)+'-CO-tdepletion.pdf',dpi=150,bbox_inches="tight")    


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


def plot_sky_positions_with_filaments(multicolor=True,plotlegend=True,plotN5566=False,zoom_virgo3=False):

    plt.figure(figsize=(11,5))
    plt.subplots_adjust(right=.95,top=.95)
    if zoom_virgo3:
        plotone="VirgoIII"
    else:
        plotone=None
    plot_spines(multicolor=multicolor,colorone='VirgoIII',color='m',legend=plotlegend,plotone=plotone)
    plt.plot(v.main['RA'],v.main['DEC'],'k.',alpha=.4,label='Virgo Filament Catalog')
    # plot CO sources
    flag = v.main['HAobsflag'] #& ~v.main['COflag']
    #plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'cs',alpha=.5,mfc='None',label=r'$H\alpha$',markersize=8)
    flag = v.groupMembs
    plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'cs',mec='k', color='#00FFFF',alpha=.8,label=r'$\rm NGC5364~Group$',markersize=10)    

    if plotN5566:
        flag = v.groupMembs_n5566
        plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'m^',color="chartreuse", mec='k',alpha=.8,label=r'$\rm NGC5566~Group$',markersize=10)    
    

    

    #flag = v.main['COflag'] & ~v.main['HAobsflag']
    #plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'bo',color='purple',markersize=8,alpha=.7,label=r'$CO, No \ H\alpha$',mec='0.5')
    #flag = v.main['COflag'] & v.main['HAobsflag']
    #plt.plot(v.main['RA'][flag],v.main['DEC'][flag],'co',markersize=5,alpha=.7,label=r'$CO + H\alpha$',mec='0.5')
    plt.axis([118,262,-1,65])
    plt.legend(loc='lower right',fontsize=14)#bbox_to_anchor=(1.01,1),loc='upper left')

    if zoom_virgo3:
        plt.axis([203,231,-3,10])


    plt.gca().invert_xaxis()
    plt.xlabel('RA (deg)',fontsize=20)
    plt.ylabel('DEC (deg)',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.legend(loc='upper right')
    #plt.title('Filamentary Structures Surrounding the Virgo Cluster',fontsize=18)
    if plotN5566:
        outfile = plotdir+'/ngc5364_ngc5566_positions_filaments.png'
    else:
        outfile = plotdir+'/ngc5364_positions_filaments.png'
    if zoom_virgo3:
        outfile = outfile.replace(".png","zoomv3.png")
    
    plt.savefig(outfile,dpi=300,bbox_inches="tight")
    plt.savefig(outfile.replace(".png",".pdf"))


class grouptables(vtables):
    def get_group_members(self):
        self.groupMembs = self.kt['PGC1'] == 49547
        print(f"number of Kourkchi group members = {np.sum(self.groupMembs)}")
        self.groupMembsHA = (self.kt['PGC1'] == 49547) & (self.main['HAobsflag'])

    def get_group_members_ngc5566(self):
        self.groupMembs_n5566 = self.kt['PGC1'] == 51233
        print(f"number of Kourkchi members in NGC5566 group = {np.sum(self.groupMembs_n5566)}")
        self.groupMembsHA_n5566 = (self.kt['PGC1'] == 51233) & (self.main['HAobsflag'])
        print(f"number of Kourkchi members in NGC5566 group w/Ha obs= {np.sum(self.groupMembsHA_n5566)}")
        
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
        tc = [r'\end{center} \n \tablecomments{$^a$ Hubble T-type from Hyperleda \citep{Makarov2014}. \\']
        tc.append(r'$^b$ Heliocentric recession velocity. \\')
        tc.append(r'$^d$ Flow-corrected, 3D distance from the center of the Virgo Cluster from \citet{Castignani2022a} \\')
        tc.append(r'$^d$ \ha \ flux correction, equal to the inverse of the filter transmission at the redshift of galaxy. \\')        
        tc.append(r'$^e$ Stellar mass from MAGPHYS SED fitting presented in \citet{Conger2025}. \\')
        tc.append(r'$^f$ SFR from MAGPHYS SED fitting presented in \citet{Conger2025}. \\')
        tc.append(r'$^g$ sSFR from MAGPHYS SED fitting presented in \citet{Conger2025}. \\')
        tc.append(r'$^h$ $H_2$ deficiency from \citet{Castignani2022a}.\\')

        # end table comments
        tc.append(r'}\\')

        print(tc)
        col_names = ['VFID','NED Name','$T^{a}$', '$v_r^{b}$','Virgo $d_{3d}^{c}$','$1/T(z)^{d}$', '$\log(M_\star)^e$',\
                     '$\log(SFR)^f$','$\log(sSFR)^g$','H2 def$^h$']
        col_formats={'$T^{a}$':'%5.1f', '$v_r^{b}$':'%.0f','Virgo $d_{3d}^{c}$':'%.1f' ,'$1/T(z)^{d}$': '%5.2f',\
                     '$\log(M_\star)^e$': '%5.2f',\
                     '$\log(SFR)^f$': '%5.2f',\
                     '$\log(sSFR)^g$': '%5.2f',\
                     'H2 def$^h$': '%5.2f',\
                     'HI def': '%5.2f',\
                     #'HI def Bos': '%5.2f',\
                     }
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
        paperTab = Table([self.main['VFID'],self.main['NEDname'],self.hl['t'],self.main['vr'],self.dist3dVirgo,self.halpha['FILT_COR'],mass_column,self.magphys['logSFR_med'],self.magphys['logsSFR_med'],self.paper1['H2def']],names=col_names)[flag]        


        
        ################################################
        # specify units
        ################################################        
        paperTab['$v_r^{b}$'].unit = u.km/u.s
        paperTab['Virgo $d_{3d}^{c}$'].unit = u.Mpc
        paperTab['$\log(M_\star)^e$'].unit = u.Msun
        paperTab['$\log(SFR)^f$'].unit = u.Msun/u.yr
        paperTab['$\log(sSFR)^g$'].unit = u.Msun/u.yr/u.Msun        


        ################################################
        # write table
        ################################################        
        
        if outfile is None:
            outfile = plotdir+'/NGC5364_tab1.tex'
        paperTab.write(outfile,format='latex',names=col_names,formats=col_formats,latexdict=latexdict,\
                       fill_values=[(ascii.masked,'\\nodata')],overwrite=True)#,(np.ma.masked,'no data')])
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
    v.get_group_members_ngc5566()        
    groupMembs = v.kt['PGC1'] == 49547
