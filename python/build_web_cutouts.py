#!/usr/bin/env python

'''
GOAL:
* create web page to inspect the cutouts

USAGE:
* run from cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)
* https://docs.astropy.org/en/stable/visualization/normalization.html#:~:text=The%20astropy.visualization%20module%20provides%20a%20framework%20for%20transforming,x%20represents%20the%20values%20in%20the%20original%20image%3A

TO DO:
* need to figure out how to handle repeated observations
  - don't overwrite directory

* fix how I combine unwise images when multiple images are returned
* same for galex
'''

import os
import numpy as np
import glob
from scipy.stats import scoreatpercentile
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch, PercentileInterval
from astropy.visualization import ImageNormalize
from astropy.visualization import LinearStretch,SinhStretch
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip

from PIL import Image

homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')
import rungalfit as rg #This code has galfit defined functions 

#from build_web_coadds import get_galaxies_fov, plot_vf_gals
from build_web_common import *
###########################################################
####  GLOBAL VARIABLES
###########################################################


VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_north_v2_main.fits'
haimaging_path = os.path.join(homedir,'github/HalphaImaging/python3/')
#sys.path.append(haimaging_path)

#vfmain = fits.getdata(VFMAIN_PATH)
residual_stretch = LinearStretch(slope=0.5, intercept=0.5) + SinhStretch() + \
    LinearStretch(slope=2, intercept=-1)
###########################################################
####  FUNCTIONS
###########################################################

def display_image(image,percentile1=.5,percentile2=99.5,stretch='asinh',mask=None,sigclip=True):
    lowrange=False
    # use inner 80% of image
    xdim,ydim = image.shape
    xmin = int(.1*xdim)
    xmax = int(.9*xdim)    
    ymin = int(.1*ydim)
    ymax = int(.9*ydim)
    if mask is not None:
        imdata = np.ma.array(image,mask=mask)
        
    else:
        imdata = image
    v1 = scoreatpercentile(imdata,percentile1)    
    v2 = scoreatpercentile(imdata,percentile2)
    
    if mask is not None:
        mask = mask[xmin:xmax,ymin:ymax]

    if sigclip:
        clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5)#,grow=3)
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]

    if mask is not None:
        norm = simple_norm(clipped_data[mask], stretch=stretch,max_percent=percentile2)
    else:
        norm = simple_norm(clipped_data, stretch=stretch,max_percent=percentile2)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower',vmin=v1,vmax=v2)
    

def make_png(fitsimage,outname):
    imdata,imheader = fits.getdata(fitsimage,header=True)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(projection=wcs.WCS(imheader))
    plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
    display_image(imdata,sigclip=False)
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('DEC (deg)',fontsize=16)        
    plt.savefig(outname)        
    plt.close(fig)

def display_galfit_model(galfile,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis',zoom=None,outdir=None):
      '''
      ARGS:
      galfile = galfit output image (with image, model, residual)
      percentile1 = min percentile for stretch of image and model
      percentile2 = max percentile for stretch of image and model
      p1residual = min percentile for stretch of residual
      p2residual = max percentile for stretch of residual
      cmap = colormap, default is viridis
      '''
      # model name


      image,h = fits.getdata(galfile,1,header=True)
      model = fits.getdata(galfile,2)
      residual = fits.getdata(galfile,3)

      if zoom is not None:
         print("who's zoomin' who?")
         # display central region of image

         # get image dimensions and center
         xmax,ymax = image.shape
         xcenter = int(xmax/2)
         ycenter = int(ymax/2)

         # calculate new size to display based on zoom factor
         new_xradius = int(xmax/2/(float(zoom)))
         new_yradius = int(ymax/2/(float(zoom)))

         # calculate pixels to keep based on zoom factor
         x1 = xcenter - new_xradius
         x2 = xcenter + new_xradius
         y1 = ycenter - new_yradius
         y2 = ycenter + new_yradius
         
         # check to make sure limits are not outsize image dimensions
         if (x1 < 1):
            x1 = 1
         if (y1 < 1):
            y1 = 1
         if (x2 > xmax):
            x2 = xmax
         if (y2 > ymax):
            y2 = ymax

         # cut images to new size
         image = image[x1:x2,y1:y2]
         model = model[x1:x2,y1:y2]
         residual = residual[x1:x2,y1:y2]         
         pass
      imwcs = wcs.WCS(h)
      images = [image,model,residual]
      titles = ['image','model','residual']
      v1 = [scoreatpercentile(image,percentile1),
            scoreatpercentile(image,percentile1),
            scoreatpercentile(residual,p1residual)]
      v2 = [scoreatpercentile(image,percentile2),
            scoreatpercentile(image,percentile2),
            scoreatpercentile(residual,p2residual)]
      norms = [simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(residual,'linear',max_percent=p2residual)]
      
      outim = ['galfit_image.png','galfit_model.png','galfit_residual.png']
      if outdir is not None:
          outim = [os.path.join(outdir,f) for f in outim]
      for i,im in enumerate(images):
          fig = plt.figure(figsize=(6,6))          
          plt.subplot(1,1,1,projection=imwcs)
          plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
          plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
          plt.xlabel('RA (deg)',fontsize=16)
          plt.ylabel('DEC (deg)',fontsize=16)
          #plt.title(titles[i],fontsize=16)
          plt.savefig(outim[i])
          plt.close(fig)
###########################################################
####  CLASSES
###########################################################

    
    
class cutout_dir():

    def __init__(self,cutoutdir=None,outdir=None):
        '''
        INPUT:
        * directory containing cutouts
        * output directory for png images

        This creates the png images for different cutouts
        '''
        self.gname = os.path.basename(os.path.dirname(cutoutdir))
        self.vfid = self.gname.split('-')[0]
        print('inside cutoutdir, gname = ',self.gname)
        print('cutoutdir = ',cutoutdir)
        print('outdir = ',outdir)        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir
        
        
        

        self.cutoutdir = cutoutdir
    def runall(self):
        self.get_halpha_names()
        try:
            self.get_legacy_names()
            self.legacy_flag = True
            self.get_legacy_jpg()                    
        except IndexError:
            print('WARNING: problem with legacy images')
            self.legacy_flag = False
        try:
            self.get_wise_names()
        except:
            print('WARNING: problem with wise images')
            self.wise_flag = False
        try:
            self.get_galex_names()
        except:
            print('WARNING: problem with wise images')
            self.galex_flag = False
        self.make_png_plots()
        self.make_cs_png()        
        self.get_galfit_model()
        try:
            self.get_phot_tables()
        except FileNotFoundError:
            print('WARNING: no phot files found - check this out')
    def get_halpha_names(self):
        search_string = os.path.join(self.cutoutdir,self.gname+'*-R.fits')
        print(search_string)
        t = glob.glob(search_string)
        print(t)
        
        self.rimage = t[0]
        self.haimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-Ha.fits'))[0]
        self.csimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-CS.fits'))[0]
        self.maskimage = self.rimage.replace('.fits','-mask.fits')
    def get_legacy_names(self):
        ''' get names of legacy images  '''
        legdir = os.path.join(self.cutoutdir,'legacy')
        self.legacy_g = glob.glob(legdir+'/*-g.fits')[0]
        self.legacy_r = glob.glob(legdir+'/*-r.fits')[0]
        self.legacy_z = glob.glob(legdir+'/*-z.fits')[0]
    def get_legacy_jpg(self):
        ''' copy jpg to local directory '''
        legdir = os.path.join(self.cutoutdir,'legacy')        
        legacy_jpg = glob.glob(legdir+'/*.jpg')[0]
        jpeg_data = Image.open(legacy_jpg)
        fig = plt.figure(figsize=(6,6))


        #hdu = fits.open(filename)[0]
        header = fits.getheader(self.legacy_g)
        imwcs = wcs.WCS(header)
        plt.subplot(projection=imwcs)
        plt.subplots_adjust(left=.2,bottom=.15,top=.95,right=.95)        
        plt.imshow(jpeg_data, origin='lower')
        plt.xlabel('RA (deg)',fontsize=16)
        plt.ylabel('Dec (deg)',fontsize=16)
        self.legacy_jpg = os.path.join(self.outdir,os.path.basename(legacy_jpg))        
        plt.savefig(self.legacy_jpg)
        plt.close(fig)

    def get_wise_names(self):
        ''' get names of wise images  '''
        # need to fix this to check for combined unwise images
        wisedir = os.path.join(self.cutoutdir,'unwise')
        self.w1 = glob.glob(wisedir+'/*-w1-img-m.fits')[0]
        self.w2 = glob.glob(wisedir+'/*-w2-img-m.fits')[0]
        self.w3 = glob.glob(wisedir+'/*-w3-img-m.fits')[0]
        self.w4 = glob.glob(wisedir+'/*-w4-img-m.fits')[0]                        
    def get_galex_names(self):
        ''' get names of galex images  '''        
        galdir = os.path.join(self.cutoutdir,'galex')
        galexfiles = glob.glob(galdir+'/*nuv*.fits')
        if len(galexfiles) > 0:
            for f in galexfiles:
                if f.find('nuv') > -1:
                    self.nuv = f
                    self.nuv_flag = True
        else:
            self.nuv_flag = False

    def define_png_names(self):
        pass
    def make_png_plots(self):
        self.fitsimages = [self.rimage,self.haimage,self.csimage,\
                           self.legacy_g,self.legacy_r,self.legacy_z,\
                           self.w1,self.w2,self.w3,self.w4,\
                           self.maskimage]
        if self.nuv_flag:
            self.fitsimages.append(self.nuv)
        self.pngimages=[]
        for f in self.fitsimages:
            pngfile = os.path.join(self.outdir,os.path.basename(f).replace('.fits','.png'))
            try:
                make_png(f,pngfile)
                self.pngimages.append(pngfile)                          
            except FileNotFoundError:
                print('WARNING: can not find ',f)
                self.pngimages.append(None)
            except TypeError:
                print('WARNING: problem making png for ',f)
                self.pngimages.append(None)

    def make_cs_png(self):
        csdata,csheader = fits.getdata(self.csimage,header=True)
        imx,imy,keepflag = get_galaxies_fov(self.csimage,vfmain['RA'],vfmain['DEC'])
        galsize=60/(abs(csheader['CD1_1'])*3600)        
        
        for i,p2 in enumerate([99.9,99.99]):
            fig = plt.figure(figsize=(6,6))
            plt.subplot(projection = wcs.WCS(csheader))
            plt.subplots_adjust(bottom=.15,left=.2,right=.95,top=.95)
            ax = plt.gca()
            display_image(csdata,stretch='asinh',percentile1=.5,percentile2=p2)
            # mark VF galaxies
            #plot_vf_gals(imx,imy,keepflag,vfmain,ax,galsize=galsize)
            suffix = "-{}.png".format(p2)
            pngfile = os.path.join(self.outdir,os.path.basename(self.csimage).replace('.fits',suffix))
            plt.xlabel('RA (deg)',fontsize=16)
            plt.ylabel('DEC (deg)',fontsize=16)        
            plt.savefig(pngfile)
            plt.close(fig)
            if i == 0:
                self.cs_png1 = pngfile
            elif i == 1:
                self.cs_png2 = pngfile                
    def get_galfit_model(self):
        ''' read in galfit model and make png '''
        self.galfit = self.rimage.replace('.fits','-1Comp-galfit-out.fits')
        if os.path.exists(self.galfit):
            # store fit results
            display_galfit_model(self.galfit,outdir=self.outdir)

            outim = ['galfit_image.png','galfit_model.png','galfit_residual.png']
        
            self.galimage = os.path.join(self.outdir,outim[0])
            self.galmodel = os.path.join(self.outdir,outim[1])
            self.galresidual = os.path.join(self.outdir,outim[2])        

            # store fitted parameters

            t = rg.parse_galfit_1comp(self.galfit)
        
            #header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
            self.xc, self.xc_err = t[0]
            self.yc, self.yc_err = t[1]
            self.mag, self.mag_err = t[2]
            self.re, self.re_err = t[3]
            self.nsersic, self.nsersic_err = t[4]
            self.BA, self.BA_err = t[5]
            self.PA, self.PA_err = t[6]
            self.sky, self.sky_err = t[7]
            self.error = t[8]
            self.chi2nu = t[9]
        else:
            self.galimage = None
        
        pass
    

    def get_phot_tables(self):
        ''' read in phot tables and make plot of flux and sb vs sma '''
        cs_galfit_phot = self.csimage.replace('.fits','-GAL-phot.fits')
        cs_gphot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','-GAL-phot.fits')
        r_gphot = fits.getdata(r_galfit_phot)

        # photutils flux
        cs_galfit_phot = self.csimage.replace('.fits','-phot.fits')
        cs_phot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','-phot.fits')
        r_phot = fits.getdata(r_galfit_phot)
        
        # plot enclosed flux
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
        alphas = [1,.4,.6,.4]
        for i,t in enumerate(tabs):
            y1 = t['flux_erg']+t['flux_erg_err']
            y2 = t['flux_erg']-t['flux_erg_err']
            if (i == 1) + (i == 3):
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend(loc='lower right')        
        self.efluxsma_png = os.path.join(self.outdir,self.gname+'-enclosed-flux.png')
        plt.savefig(self.efluxsma_png)
        plt.close(fig)
        
        # plot mag sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        #tabs = [r_gphot,cs_gphot]
        #labels = ['r','Halpha']
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halpha','photutil r','photutil Halpha']
        
        for i,t in enumerate(tabs):
            y1 = t['mag']+t['mag_err']
            y2 = t['mag']-t['mag_err']
            if (i == 1) + (i == 3):
                alpha=.4
            else:
                alpha=1
            
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alpha)
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('magnitude (AB)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()        
        self.emagsma_png = os.path.join(self.outdir,self.gname+'-mag-sma.png')
        plt.legend(loc='lower right')        
        plt.savefig(self.emagsma_png)
        plt.close(fig)
        
        # plot sb erg vs sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
        
        for i,t in enumerate(tabs):
            y1 = t['sb_erg_sqarcsec']+t['sb_erg_sqarcsec_err']
            y2 = t['sb_erg_sqarcsec']-t['sb_erg_sqarcsec_err']
            if (i == 1) + (i == 3):
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('SB (erg/s/cm^2/Hz/arcsec^2)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend()
        self.sbfluxsma_png = os.path.join(self.outdir,self.gname+'-sb-sma.png')
        plt.savefig(self.sbfluxsma_png)
        plt.close(fig)
                             
        # plot mag sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halpha','photutil r','photutil Halpha']
        
        for i,t in enumerate(tabs):
            y1 = t['sb_mag_sqarcsec']+t['sb_mag_sqarcsec_err']
            y2 = t['sb_mag_sqarcsec']-t['sb_mag_sqarcsec_err']
            
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()
        plt.legend()        
        self.sbmagsma_png = os.path.join(self.outdir,self.gname+'-sbmag-sma.png')
        plt.savefig(self.sbmagsma_png)
        plt.close(fig)
    
class build_html_cutout():

    def __init__(self,cutoutdir,outdir,previous=None,next=None,tel=None,run=None):
        ''' pass in instance of cutout_dir class and output directory '''
        self.cutout = cutoutdir


        outfile = os.path.join(outdir,self.cutout.gname+'.html')

        vfindices = np.arange(len(vfmain))
        self.vfindex = vfindices[vfmain['VFID'] == self.cutout.vfid]
        #print('inside build html')
        #print('coutdir = ',coutdir)
        print('outfile = ',outfile)        
        self.html = open(outfile,'w')
        self.htmlhome = 'index.html'
        self.next = next
        self.previous = previous

        self.telescope = tel
        self.run=run
        
        # for reference, this is the order of the png images
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\
        #              self.legacy_g,self.legacy_r,self.legacy_z,\
        #              self.w1,self.w2,self.w3,self.w4]

        
        #self.build_html()
    def build_html(self):
        self.write_header()
        self.write_navigation_links()
        self.write_image_stats()
        self.write_legacy_images()        
        self.write_sfr_images()
        self.write_wise_images()
        self.write_halpha_images()
        if self.cutout.galimage is not None:
            self.write_galfit_images()
            self.write_galfit_table()
        try:
            self.write_phot_profiles()
        except AttributeError:
            pass
        self.write_mag_table()
        self.write_morph_table()        
        self.write_navigation_links()
        self.close_html()
    def write_header(self):
        # title
        # home link
        # previous
        # next
        self.html.write('<html><body>\n')
        self.html.write('<style type="text/css">\n')
        self.html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        self.html.write('</style>\n')

    def write_navigation_links(self):
        # Top navigation menu--
        self.html.write('<h1>{}</h1>\n'.format(self.cutout.gname))

        self.html.write('<a href="../{}">Home</a>\n'.format(self.htmlhome))
        self.html.write('<br />\n')
        if self.previous is not None:
            previoushtml = "{}.html".format(self.previous)
            self.html.write('<a href="../{}/{}">Previous ({})</a>\n'.format(self.previous,previoushtml, self.previous))
            self.html.write('<br />\n')        
        if self.next is not None:
            nexthtml = "{}.html".format(self.next)
            self.html.write('<a href="../{}/{}">Next ({})</a>\n'.format(self.next,nexthtml,self.next))
            self.html.write('<br />\n')

    def write_image_stats(self):
        self.html.write('<h2>Image Statistics</h2>\n')        
        labels=['Telescope','Run','Pointing','R FWHM <br> (arcsec)','Halpha FWHM <br> (arcsec)','Filter Ratio','Filter Correction']
        myrow = vfha[self.vfindex]
        #print('HALPHA FILE')
        #print(self.vfindex)
        #print(myrow)
        colnames = ['POINTING','R_FWHM','H_FWHM','FILTER_RATIO','FILT_COR']
        
        data = [myrow['POINTING'][0], "{:.2f}".format(myrow['R_FWHM'][0]),\
                "{:.2f}".format(myrow['H_FWHM'][0]),\
                "{:.4f}".format(myrow['FILTER_RATIO'][0]),\
                "{:.2f}".format(myrow['FILT_COR'][0])]
        print('self.run = ',self.run)
        data.insert(0,self.run)
        data.insert(0,self.telescope)
        write_text_table(self.html,labels,data)        
    def write_sfr_images(self):
        ''' legacy jpeg, galex nuv, halpha, w4 '''
        self.html.write('<h2>Star-Formation Indicators</h2>\n')
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\ # 0,1,2
        #                   self.legacy_g,self.legacy_r,self.legacy_z,\ # 3,4,5
        #                   self.w1,self.w2,self.w3,self.w4] # 6,7,8,9

        if self.cutout.nuv_flag:
            images = [self.cutout.pngimages[-1],\
                      self.cutout.cs_png1,\
                      self.cutout.pngimages[8],self.cutout.pngimages[9]]
            labels = ['NUV','Halpha','W3','W4']
        else:
            images = [self.cutout.legacy_jpg,\
                      self.cutout.pngimages[2],\
                      self.cutout.pngimages[8],self.cutout.pngimages[9]]                      

            labels = ['Legacy','Halpha','W3','W4']
        images = [os.path.basename(i) for i in images]            
        write_table(self.html,images=images,labels=labels)
        pass

    
    def write_halpha_images(self):
        '''  r, halpha, cs, and mask images '''
        self.html.write('<h2>Halpha Images</h2>\n')        
        images = [self.cutout.pngimages[0],self.cutout.pngimages[1],self.cutout.cs_png1,self.cutout.cs_png2]
        images = [os.path.basename(i) for i in images]

        labels = ['R','Halpha+Continuum','CS, stretch 1','CS, stretch 2']
        write_table(self.html,images=images,labels=labels)
        pass

    def write_legacy_images(self):
        ''' jpg, g,r,z legacy images '''
        self.html.write('<h2>Legacy Images</h2>\n')

        images = self.cutout.pngimages[3:6]
        images = [os.path.basename(i) for i in images]        
        images.insert(0,os.path.basename(self.cutout.legacy_jpg))        
        labels = ['grz','g','r','z']
        write_table(self.html,images=images,labels=labels)

    def write_wise_images(self):
        ''' w1 - w4 images '''
        self.html.write('<h2>WISE Images</h2>\n')                
        pngimages = self.cutout.pngimages[6:10]
        wlabels = ['W1','W2','W3','W4']
        images=[]
        labels=[]
        for i,im in enumerate(pngimages):
            if im is not None:
                images.append(os.path.basename(im))
                labels.append(wlabels[i])

        write_table(self.html,images=images,labels=labels)
    
    def write_galex_images(self):
        ''' right now just nuv '''
        self.html.write('<h2>GALEX Images</h2>\n')                
        pass

    def write_galfit_images(self):
        ''' display galfit model and fit parameters for r-band image '''
        self.html.write('<h2>GALFIT r-band Modeling </h2>\n')                
        images = [self.cutout.galimage,self.cutout.galmodel,self.cutout.galresidual,\
                  self.cutout.pngimages[10]]
        images = [os.path.basename(i) for i in images]        
        labels = ['Image', 'Model', 'Residual','Mask']
        write_table(self.html,images=images,labels=labels)
        pass
    
    def write_galfit_table(self):
        ''' display galfit model and fit parameters for r-band image '''
        self.html.write('<h4>GALFIT Sersic Parameters</h4>\n')                
        labels=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        data = ['{:.1f}+/-{:.1f}'.format(self.cutout.xc, self.cutout.xc_err),
                '{:.1f}+/-{:.1f}'.format(self.cutout.yc, self.cutout.yc_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.mag, self.cutout.mag_err), 
                '{:.2f}+/-{:.2f}'.format(self.cutout.re, self.cutout.re_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.nsersic, self.cutout.nsersic_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.BA, self.cutout.BA_err),
                '{:.1f}+/-{:.1f}'.format(self.cutout.PA, self.cutout.PA_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.sky, self.cutout.sky_err),
                '{}'.format(self.cutout.error),
                '{:.2f}'.format(self.cutout.chi2nu)]
        
        #print(data)
        write_text_table(self.html,labels,data)
        pass

    def write_phot_profiles(self):
        ''' photometry table with galfit and photutil results '''
        self.html.write('<h2>Elliptical Photometry</h2>\n')
        self.html.write('<p>using galfit and photutil geometry</p>\n')                        
        images = [self.cutout.efluxsma_png,self.cutout.emagsma_png,self.cutout.sbfluxsma_png,self.cutout.sbmagsma_png]
        images = [os.path.basename(i) for i in images]
        labels = ['Enclosed Flux','Enclosed Magnitude','Surface Brightness','Surface Brightness']
        write_table(self.html,images=images,labels=labels)

    def write_mag_table(self):
        self.html.write('<h2>R-band AB Magnitudes</h2>\n')        
        labels=['mag24','mag25', 'mag26','mag_petro','mag_galfit','mag_photutil'] 
        myrow = vfha[self.vfindex]
        colnames = ['M24','M25','M26','PETRO_MAG','GAL_MAG','ELLIP_SUM_MAG']
        data = ["{:.2f}".format(myrow[c][0]) for c in colnames]
        write_text_table(self.html,labels,data)        

        self.html.write('<h2>Halpha-band AB Magnitudes/SFR</h2>\n')        
        labels=['mag16','mag17','mag_petro','SSFR_IN','SSFR_OUT'] 
        myrow = vfha[self.vfindex]
        colnames = ['HM16','HM17','HPETRO_MAG','SSFR_IN','SSFR_OUT']
        data = ["{:.2f}".format(myrow[c][0]) for c in colnames]
        write_text_table(self.html,labels,data)        
        
    def write_morph_table(self):
        self.html.write('<h2>Morphology Parameters</h2>\n')        
        labels=['Band','Gini','Asym','C30','Petro Con'] 
        myrow = vfha[self.vfindex]
        colnames = ['ELLIP_GINI','ELLIP_ASYM','C30','PETRO_CON']
        colnames2 = ['ELLIP_GINI2','ELLIP_HASYM','HC30','HPETRO_CON']
        data = ["{:.2f}".format(myrow[c][0]) for c in colnames]
        data.insert(0,'r')
        data2 = ["{:.2f}".format(myrow[c][0]) for c in colnames2]
        data2.insert(0,'Halpha')
        write_text_table(self.html,labels,data,data2=data2)
    def close_html(self):
        self.html.close()
# wrap

if __name__ == '__main__':
    # work through coadd directory
    #global vmain
    vfmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')
    vfha = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_ha.fits')    
    #cutout_source_dir = '/home/rfinn/research/Virgo/gui-output-2019-june/cutouts/'
    #telescope='INT'
    #run='2019-Jun'
    

    cutout_source_dir = '/home/rfinn/research/Virgo/gui-output-2019/cutouts/'
    telescope='INT'
    run='2019-Feb'
    os.chdir(cutout_source_dir)
    hdiflag = False
    
    # setup for HDI cutouts!
    cutout_source_dir = '/home/rfinn/research/Virgo/gui-output-all-hdi-20210308/cutouts/'
    cutout_source_dir = '/home/rfinn/research/Virgo/all-cutouts-20210309/'    
    telescope='WIYN'
    run='HDI'
    os.chdir(cutout_source_dir)
    hdiflag = True


    # setup for NGC cutouts!
    #cutout_source_dir = '/home/rfinn/research/Virgo/gui-output-NGC5846/cutouts/'
    telescope='WIYN'
    run='MOSAIC'
    os.chdir(cutout_source_dir)
    
    
    outdir = homedir+'/research/Virgo/html-dev/cutouts/'
    #outdir = homedir+'/research/Virgo/html-dev/cutouts/ngc5846-cutouts/'    

    # this should contain a list of all the galaxy folders
    flist1 = os.listdir(cutout_source_dir)

    flist1.sort()
    galindex = np.arange(len(flist))
    ngal = 0
    startgal = None
    startindex = 0
    if startgal is not None:
        for i,f in enumerate(flist1):
            if f.startswith(startgal):
                startindex=i
                print('starting with galaxy ',startgal)
                break
    for i in galindex[startindex:]: # loop through list
        subdir = flist[i]
        print(subdir)
        #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
        if os.path.isdir(subdir):
            print('##########################################')
            print('##########################################')        
            print('WORKING ON DIRECTORY: ',subdir)
            print('##########################################')
            print('##########################################')

            try:
                # move to subdirectory
                # adding the telescope and run so that we don't write over
                # images if galaxy was observed more than once
                gal_outdir = os.path.join(outdir,subdir+"")
                print('out directory for this galaxy = ',gal_outdir)
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                if not os.path.exists(gal_outdir):
                    os.mkdir(gal_outdir)
                cutoutdir = os.path.join(cutout_source_dir,subdir,"")
                p = cutout_dir(cutoutdir=cutoutdir,outdir=gal_outdir)
                p.runall()
            
                # define previous gal for html links
                if i > 0:
                    previous = (flist1[i-1])
                    #print('previous = ',previous)
                else:
                    previous = None
                if i < len(flist1)-1:
                    next = flist1[i+1]
                    #print('next = ',next)
                else:
                    next = None
                h = build_html_cutout(p,gal_outdir,previous=previous,next=next,tel=telescope,run=run)
                h.build_html()
            except:
                print('WARNING: problem building webpage for ',subdir)
        # just running on one directory for testing purposes
        #if i > 1:
        #    break

