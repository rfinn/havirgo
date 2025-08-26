#!/usr/bin/env python

"""
The goal is to use the reprojected legacy images  to make a stellar mass image.

log(M/L) = a_r + b_r x (g-r)

we already have a g-r image in the legacy subdirectory

need to convert r-band image to absolute magnitude - this is just a ZP offset


m = 22.5 - 2.5*log10(flux)

m - M = 5 log10(d_pc) - 5

to get distance, we use the hubble law and recession velocity:



    def taylor_mstar(self):
        # calc stellar mass
        logMstarTaylor=1.15+0.70*(self.a100sdss['gmi_corr']) -0.4*(self.a100sdss['absMag_i_corr'])
        ###
        # NSA abs mags are for H0=100, need to correct for our assumed cosmology
        # so need to add 5*np.log10(h), where h = cosmo.H(0)/100.
        #logMstarTaylor = logMstarTaylor - 0.4*(5*np.log10(cosmo.H(0).value/100.))
        # not doing this b/c absMag_i_corr already has distance corrected for H0

        # -0.68 + .7*gmi_cor + (Mi-4.56)/-2.5
        # add taylor_mstar column to a100sdss
        # only set values for galaxies with photflag == 1

Absolute magnitude of the sun:
https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf

Filter   Abs_Vega  Abs_AB
-------------------------
SDSS_r    4.53     4.65
DES_r     4.45     4.61
PS_r      4.53     4.64


TODO : use the legacy r-band image instead of our r-band image

"""
import os
import glob
import sys
import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
from astropy import stats, convolution

#from matplotlib import pyplot as plt
######################################################################
###  CONSTANTS
######################################################################

# coefficients for (g-r) to M/L ratio conversion
# Roediger & Courteau 2015, with Marigo TPAGB
ar = -0.647
br = 1.497

######################################################################
###  FILTER DEFINITIONS
######################################################################
# Halpha filter width in angstrom
filter_width_AA = {'BOK':80.48,'HDI':80.48,'INT':95,'MOS':80.48,'INT6657':80}

# central wavelength in angstroms
filter_lambda_c_AA = {'BOK':6620.52,'HDI':6620.52,'INT':6568,'MOS':6620.52,'INT6657':6657}

######################################################################
###  FUNCTIONS
######################################################################

def get_params_from_cutout_directory_name(dirname):
    t = os.path.basename(dirname).split('-')
    pointing = t[-1]
    dateobs = t[-2]
    telescope = t[-3]
    vfid = t[0]
    vfid_ned = os.path.basename(dirname).split('-'+telescope)[0]
    #print(telescope,dateobs,pointing,vfid_ned)
    return telescope,dateobs,pointing,vfid_ned

######################################################################
###  GALAXY CLASS
######################################################################

class galaxy():

    def __init__(self,dirname,vr,vcosmic,ra=None,dec=None):
        self.dirname = dirname
        self.vfid = dirname.split('-')[0]
        self.ra = ra
        self.dec = dec
        self.d_vr = cosmo.luminosity_distance(vr/3.e5) # dist in Mpc from recession velocity in km/s
        self.d_vcosmic = cosmo.luminosity_distance(vcosmic/3.e5) # dist in Mpc from flow-corrected recession vel in km/s
        
        self.telescope,self.dateobs,pointing,self.vfid_ned = \
            get_params_from_cutout_directory_name(dirname)
    def construct_filenames(self):
        """ build filename of relevant images """
        
        # reprojected legacy r-band image
        legacy_images = glob.glob(os.path.join('legacy/*r-ha.fits'))
        self.legacy_r2ha = legacy_images[0]

        # reprojected legacy r-band image
        legacy_images = glob.glob(os.path.join('legacy/*r.fits'))
        self.legacy_r = legacy_images[0]
        
        # reprojected legacy g-r image
        legacy_images = glob.glob(os.path.join('legacy/*gr-ha-smooth.fits'))
        self.legacy_gr = legacy_images[0]
        
        # our r-band image
        self.R = self.dirname+'-R.fits'
        # our gr-CS image (halpha image using color-based subtraction)

        # updating to use Gautam's auto cont sub image
        self.ha = self.dirname+'-CS-gr-auto.fits'

        # define the mask file
        maskfile = self.R.replace('-R.fits','-R-mask.fits')
        if not os.path.exists(maskfile):
            print("WARNING: no mask found")
            self.mask = None
        else:
            mask = fits.getdata(maskfile)
            self.mask = mask > 0
        
        pass
    
    def get_mstar_image(self,makeplots=False):
        """ create mstar image from rband and g-r image and vr  """
        # smoothing was 20, setting it to 5 to see if the effective radius changes
        smoothing=15
        rhdu = fits.open(self.legacy_r2ha)
        rsmooth = convolution.convolve_fft(rhdu[0].data, convolution.Box2DKernel(smoothing), allow_huge=True, nan_treatment='interpolate')

        # g-r image is in magnitudes already
        grhdu = fits.open(self.legacy_gr)
        mag_gr = grhdu[0].data

        # convert to absolute magnitude
        Mr_vr = 22.5 -2.5*np.log10(rsmooth) - 5*np.log10(self.d_vr.to('pc').value) + 5
        Mr_vcosmic = 22.5 -2.5*np.log10(rsmooth) - 5*np.log10(self.d_vcosmic.to('pc').value) + 5        


 
        self.logMstar_vr = ar + br*mag_gr + (Mr_vr - 4.61)/(-2.5)
        self.logMstar_vcosmic = ar + br*mag_gr + (Mr_vcosmic - 4.61)/(-2.5)

        # save mstar images
        outimage = self.dirname+'-logmstar-vr.fits'
        outimage = self.dirname+'-mstar-vr.fits'        
        hdu = fits.PrimaryHDU(10.**(self.logMstar_vr-7), header=rhdu[0].header) # making in linear units of 1e7 Msun
        hdu.writeto(outimage, overwrite=True) #sky-subtracted r-band image - use this for photometry

        outimage = self.dirname+'-logmstar-vcosmic.fits'
        outimage = self.dirname+'-mstar-vcosmic.fits'        
        hdu = fits.PrimaryHDU(10.**(self.logMstar_vcosmic-7), header=rhdu[0].header)
        hdu.writeto(outimage, overwrite=True) #sky-subtracted r-band image - use this for photometry

        rhdu.close()
        grhdu.close()
        
        
        pass

    def get_sfr_image(self,makeplots=False,verbose=True):
        """ read in gr-CS image and convert to SFR"""

        # read in continuum-subtracted image
        hahdu = fits.open(self.ha)
        self.haheader = hahdu[0].header
        hZP = hahdu[0].header['PHOTZP']


        # see notes in the overleaf document
        # https://www.overleaf.com/project/5ede62d7c2853b0001731d73

        # combination of conversion from AB to cgs, and then to SFR
        # 60.710 = 48.6/2.5 (from fnu_cgs to mAB conversion) + 41.27 (kennicutt & evans conversion)
        a = 10.**(-0.4*hZP - 60.710)

        # conversion from fnu to flambda to flux
        dlambda = filter_width_AA[self.telescope]
        clambda = filter_lambda_c_AA[self.telescope]
        # c = 3e18 A/s
        # dlambda = filter width in A
        # clambda = center wavelength in A
        b = 3.e18*dlambda/clambda**2 # frequency in Hz

        # convert flux to lumininosity
        c_vr = 4*np.pi*self.d_vr.cgs.value[0]**2
        c_vcosmic = 4*np.pi*self.d_vcosmic.cgs.value[0]**2
        # from halphagui

        #print(a)
        #print(b)
        #print(c_vr)
        if verbose:
            print(f"scale factors for sfr image: a={a:3.2e}, b={b:3.2e}, c={c_vr:3.2e}, product={a*b*c_vr:3.2e}")

        product = a*b*c_vr
        self.sfr_vr = hahdu[0].data*product

        product = a*b*c_vcosmic        
        self.sfr_vcosmic = hahdu[0].data*product
        #plt.figure()
        #plt.subplot(1,2,1)
        #plt.imshow(hahdu[0].data)
        #plt.subplot(1,2,2)
        #plt.imshow(self.sfr_vr)

        # write out images
        outimage = self.dirname+'-sfr-vr.fits'
        # why am I multiplying by a factor of 1000?
        #hdu = fits.PrimaryHDU(self.sfr_vr*1.e3, header=hahdu[0].header)
        hdu = fits.PrimaryHDU(self.sfr_vr, header=hahdu[0].header)
        hdu.writeto(outimage, overwrite=True) 
        
        outimage = self.dirname+'-sfr-vcosmic.fits'
        #hdu = fits.PrimaryHDU(self.sfr_vcosmic*1.e3, header=hahdu[0].header)
        hdu = fits.PrimaryHDU(self.sfr_vcosmic, header=hahdu[0].header)
        hdu.writeto(outimage, overwrite=True) 
        hahdu.close()        
        pass

    def get_ssfr_image(self):
        """ divide SFR image by mstar image  """
        from astropy.stats import sigma_clipped_stats
        #import ccdproc
        try:
            from photutils import make_source_mask
            # create mask to cut low SNR pixels based on SNR in SFR image
            mask = make_source_mask(self.sfr_vr,nsigma=3,npixels=5,dilate_size=5)
        except ImportError:
            from photutils.segmentation import detect_threshold, detect_sources
            from photutils.utils import circular_footprint
            footprint = circular_footprint(radius=10)
            # create threshold
            threshold = detect_threshold(self.sfr_vr,nsigma=3)#,dilate_size=5)
            # detect sources
            segment_img = detect_sources(self.sfr_vr, threshold, npixels=10)
            # make source mask
            mask = segment_img.make_source_mask(footprint=footprint)
            
        masked_data = np.ma.array(self.sfr_vr,mask=mask)
        #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

        mean,median,std = sigma_clipped_stats(masked_data,sigma=3.0,cenfunc=np.ma.mean)
        print(f'STD cut in SF image = {std:.3e}')
        flag = self.sfr_vr < 3*std
        self.logssfr = np.log10(self.sfr_vr) - self.logMstar_vr
        self.logssfr[flag] = np.nan
        self.haheader['SFRSTD']=float(f"{std:.3e}")
        outimage = self.dirname+'-ssfr.fits'
        # changing to a linear scale on May 30, 2025
        hdu = fits.PrimaryHDU(10.**(self.logssfr+10), header=self.haheader) # convert to units of 1e-10
        hdu.writeto(outimage, overwrite=True) 

    def save_figure(self,zoom=True,zoomfactor=2):
        # save figure
        from matplotlib import pyplot as plt
        from scipy.stats import scoreatpercentile
        from astropy.visualization import simple_norm
        from astropy.wcs import WCS
        from PIL import Image
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        

        imwcs = WCS(self.haheader)
        
        
        fig = plt.figure(figsize=(12,4))
        plt.subplots_adjust(wspace=0.3,left=.05)


        # set up list of images to 
        images = [self.logMstar_vr,self.sfr_vr,self.logssfr]
        ximsize,yimsize = self.logMstar_vr.shape
        titles = ["logMstar",r"$H\alpha \ SFR$",'logsSFR']
        percentile1 = .5
        percentile2 = 99.5
        stretch=['linear','asinh','linear']
        for i,cs in enumerate(images):
            if self.mask is not None:
                cs = np.ma.masked_where(self.mask,cs)

            plt.subplot(1,4,i+2)#,projection=imwcs)
            try:
                norm = simple_norm(cs, stretch=stretch[i],max_percent=percentile2,min_percent=percentile1)
                plt.imshow(cs, norm=norm,origin='lower',interpolation='nearest')#,vmin=v1,vmax=v2)
            except:
                print("WARNING: problem calculating the simple_norm - check images ") 
            #if i == 2:
            #    plt.imshow(cs, norm=norm,origin='lower',vmin=-1.2e-5,vmax=1e-4)
            #if i == 3:
            #    plt.imshow(cs, norm=norm,origin='lower',vmin=-10.5,vmax=-9)
            #else:
            #    plt.imshow(cs, norm=norm,origin='lower')#,vmin=v1,vmax=v2)

            if zoom:
                # TODO fix to use headers to zoom
                # check for RA and DEC, b/c image might not be centered
                #print("zooming")                    
                xsize,ysize = cs.shape
                delta = xsize//(zoomfactor*2)                    
                if self.ra is not None:
                    galcoord = SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
                    xcenter,ycenter = imwcs.world_to_pixel(galcoord)
                    xcenter = xcenter[0]
                    ycenter = ycenter[0]                        
                else:
                    xcenter = xsize//2
                    ycenter = ysize // 2

                xmin = xcenter - delta
                xmax = xcenter + delta

                ymin = ycenter - delta
                ymax = ycenter + delta

                if xmin < 1:
                    xmin = 1
                if ymin < 1:
                    ymin = 1
                if xmax > xsize:
                    xmax=xsize-1
                if ymax > ysize:
                    xmax=xsize-1
                #print([ymin,ymax,xmin,xmax]) 
                plt.axis([xmin,xmax,ymin,ymax])                
                xcoords = np.array([xmin,xmax])
                ycoords = np.array([ymin,ymax])

            #plt.imshow(cs)#,vmin=-0.015,vmax=.1)#,cmap='gray_r')

            plt.colorbar(fraction=.045)
            plt.title(titles[i],fontsize=16)
            plt.xticks([],[])
            plt.yticks([],[])

        # this next block zooms into center half of the image

        # zoom jpeg image

        sky = imwcs.pixel_to_world(xcoords,ycoords)
        # for debugging
        #print("xmin,xmax = ",xmin,xmax,ymin,xmax)
        #print("xcoords = ",xcoords)
        #print("sky = ",sky)
        
        #print(legacyr)
        legwcs = WCS(fits.getheader(self.legacy_r))
        # read in the Legacy jpeg image
        legacy_jpg = glob.glob('legacy/*.jpg')[0]        
        jpeg_data = Image.open(legacy_jpg)        
        
        # plot jpg as projection of legacy r-band
        plt.subplot(1,4,1,projection=legwcs)
        plt.imshow(jpeg_data)
        plt.title("Legacy")
        axleg = plt.gca()
        #zoom=False
        if zoom:
            # convert ra and dec of zoomed image into pixel
            # coordinates on jpeg image using the legacy wcs
            # to translate between the world and pixel coords
            x,y = legwcs.world_to_pixel(sky)
            plt.axis([x[0],x[1],y[0],y[1]])


        plt.savefig(self.dirname+'-mstar-sfr-ssfr.png')        
        
        return fig,axleg
    
if __name__ == '__main__':
    dirname = sys.argv[1]
    topdir = os.getcwd()

    if len(sys.argv) > 2:
        zoomflag=False
    else:
        zoomflag = True

    os.chdir(dirname)
    # get VFID from dirname
    vfid = os.path.basename(dirname).split('-')[0]
    # get nedname from dirname
    # this is complicated b/c some ned names are crazy, like MCG+10-23-067

    homedir = os.getenv("HOME")
    tabledir = homedir+'/research/Virgo/tables-north/v2/'
    # read in vf main table
    vfmain = Table.read(tabledir+'vf_v2_main.fits')
    
    # get row index for this galaxy
    galindex = np.arange(len(vfmain))[vfmain['VFID'] == vfid]
    
    # get vr from vf_v2_main.fits
    vr = vfmain['vr'][galindex]
    
    # get vcosmic from vf_vf_env.fits
    vfenv = Table.read(tabledir+'vf_v2_environment.fits')
    vcosmic = vfenv['Vcosmic'][galindex]

    # get RA, DEC

    ra = vfmain['RA'][galindex]
    dec = vfmain['DEC'][galindex]
    # initiate instance of galaxy class
    g = galaxy(dirname,vr, vcosmic,ra=ra,dec=dec)
    g.construct_filenames()
    g.get_mstar_image()
    g.get_sfr_image()
    g.get_ssfr_image()
    g.save_figure(zoom=zoomflag,zoomfactor=2.)
    os.chdir(topdir)

