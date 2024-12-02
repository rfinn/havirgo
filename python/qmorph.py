#!/usr/bin/env python

"""
statmorph using a different segmentation image when running gini, which is not what we want

I want to run gini on r and halpha, using the same segmentation image for both (based on r-band image)

"""

#!/usr/bin/env python

'''
This is going to be the wrapper to do photometry on the detected galaxies.

Would like to build this totally on photutils.

Useful references:


https://photutils.readthedocs.io/en/stable/segmentation.html
- detecting sources
- creating a segmentation image
- getting source properties (including total flux, Gini coefficient!!!)
- defining elliptical apertures for sources

'''

from photutils import detect_threshold, detect_sources

# changing to remove deprecated function source_properties
#from photutils import source_properties
from photutils.segmentation import SourceCatalog

from photutils import make_source_mask
from photutils import Background2D, MedianBackground
from photutils import EllipticalAperture
from photutils.utils import calc_total_error

from photutils.isophote import EllipseGeometry, Ellipse
from photutils import aperture_photometry
from photutils.morphology import gini

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from astropy.stats import sigma_clip, SigmaClip,sigma_clipped_stats
from astropy.visualization import simple_norm


import scipy.ndimage as ndi
import statmorph
from statmorph.utils.image_diagnostics import make_figure


from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

import numpy as np
import sys

import time
start_time = time.time()

import matplotlib
#matplotlib.use('Qt5Agg')


# modules in halphagui
import imutils

## filter information
## from https://www.noao.edu/kpno/mosaic/filters/
central_wavelength = {'4':6620.52,'8':6654.19,'12':6698.53,'16':6730.72,'R':6513.5,'r':6292.28,'inthalpha':6568.,'intha6657':6657,'intr':6240} # angstrom
dwavelength = {'4':80.48,'8':81.33,'12':82.95,'16':81.1,'R':1511.3,'r':1475.17,'inthalpha':95.,'intha6657':80,'intr':1347} # angstrom

# define colors - need this for plotting line and fill_between in the same color
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def display_image(image,percent=99.9,lowrange=False,mask=None,sigclip=True):
    lowrange=False
    if sigclip:
        clipped_data = sigma_clip(image,sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    #v1,v2=scoreatpercentile(image,[.5,99.5])            
    #plt.imshow(image, cmap='gray_r',vmin=v1,vmax=v2,origin='lower')    


def get_M20(catalog,objectIndex):
    """ 
    calculate M20 according to Lotz+2004 for central object only
    https://iopscience.iop.org/article/10.1086/421849/fulltext/  

    PARAMS:
    * catalog - this is a photutils.SourceCatalog
    * objectIndex - the object number in the catalog

    RETURNS:
    * M20

    NOTES:
    cat.data_ma = A 2D MaskedArray cutout from the data using the minimal bounding box of the source.
    cat.segment = A 2D ndarray cutout of the segmentation image using the minimal bounding box of the source.
    cat.cutout_centroid = The (x, y) coordinate, relative to the cutout data, of the centroid within the isophotal source segment.

    """
    # total second-order moment Mtot is the flux in each pixel fi multiplied by distance^2 of pixel to center,
    # summed over all pixels assigned to the segmentation map

    objNumber = catalog.label[objectIndex]

    # data_ma returns a 2D array with the unmasked values
    # using the min bounding box that fits the galaxy    
    dat = catalog.data_ma[objectIndex]

    # create flag for pixels associated with object in segmentation map
    segflag = catalog.segment[objectIndex] == objNumber

    # get the center coordinates of the object
    xc,yc = catalog.cutout_centroid[objectIndex]

    # can't make sense of moments that photutils includes in the catalog, so recalculating here
    ymax,xmax = catalog.data_ma[objectIndex].shape

    # create a meshgrid to represent pixels in segmentation
    X,Y = np.meshgrid(np.arange(xmax),np.arange(ymax))

    # calculate distance of each point from center of galaxy
    # this is for 2nd order moment
    distsq = (X-xc)**2 + (Y-yc)**2

    # second Moment total
    # segflag ensures that we are just counting the pixels assoc with object
    second_moment_tot = np.sum(dat[segflag]*distsq[segflag])

    ##
    # getting the second moment of 20% highest pixels
    ##

    # using https://github.com/vrodgom/statmorph/blob/master/statmorph/statmorph.py#L1430    
    # Calculate threshold pixel value
    # TODONE - update to fix M20
    sorted_pixelvals = np.sort(dat.flatten())
    flux_fraction = np.cumsum(sorted_pixelvals) / np.sum(sorted_pixelvals)
    sorted_pixelvals_20 = sorted_pixelvals[flux_fraction >= 0.8]

    threshold_brightest20 = sorted_pixelvals_20[0]

    #threshold_brightest20 = scoreatpercentile(dat[segflag].flatten(),80)

    # define flag for pixels that contain top 20% of total flux
    brightest20 = dat > threshold_brightest20

    # sum the second moment of brightest 20
    second_moment_20 = np.sum(dat[segflag & brightest20]*distsq[segflag & brightest20])

    # now calculate M20 as
    # M20 = log10(Sum_Mi/Mtot)

    M20 = np.log10(second_moment_20/second_moment_tot)

    return M20

def get_fraction_masked_pixels(catalog,objectIndex):
    """ 
    get area in the segmentation image, 
    masked area, and fraction of pixels masked
    """

    objNumber = catalog.label[objectIndex]
    dat = catalog.data[objectIndex]    
    masked_dat = catalog.data_ma[objectIndex]

    # create flag for pixels associated with object in segmentation map
    goodflag = catalog.segment[objectIndex] == objNumber

    # get number of pixels in the original segmentation image
    number_total = np.sum(goodflag)

    number_masked = number_total - np.sum(goodflag & masked_dat.mask)

    return number_total, number_masked, number_masked/number_total

# read in image and mask

# identify source for photometry

# run detect to detect source

# estimate ellipse parameters from source properties

# run 

class qmorph():
    """
    input:
    image - image to use for measuring morphology
    segmap - segmentation map, with one object identified
    """
    def __init__(self,image,segmap):

        # make sure image and segmap have same shape

        assert image.shape == segmap.shape



class ellipse():
    '''
    class to run photometry routines on image

    INPUT
    * image         - primary image (this is usually rband)
    * image2        - image2 is designed to be the Halpha image, 
                      but it can be any second image whereby you define 
                      the ellipse geometry using image 1, and
                      measure the photometry on image 1 and image 2
    * mask          - mask to apply when measuring photometry.  
                      this is usually created from the r-band image
    * image_frame   - image frame for plotting inside a gui; like if this is called from halphamain
    * use_mpl       - use mpl for testing purposes, before integrating with pyqt gui
    * napertures    - number of apertures for measuring photmetry in. default is 20.
    * apertures     - list of apetures to use, instead of generating them automatically
    * image2_filter - this is used to calculate the flux levels in the second filter.  
                      this should be one of standard halpha filters in our dictionary 
                      ('4','inthalpha','intha6657').
    * filter_ratio  - ratio of flux in image2/image1
    * psf           - used by statmorph (not sure if this is the image or the image name)
    * psf_ha        - used by statmorph



    '''
    def __init__(self, image, image2 = None, mask = None, image_frame=None, use_mpl=False, napertures=20,apertures=None, image2_filter=None, filter_ratio=None,psf=None,psf_ha=None,objra=None,objdec=None,fixcenter=False):
        '''  inputs described above '''

        self.image, self.header = fits.getdata(image, header=True)
        self.image_name = image

        # get image dimensions - will use this to determine the max sma to measure
        self.yimage_max, self.ximage_max = self.image.shape

        self.objra = objra
        self.objdec = objdec
        self.fixcenter = fixcenter
        
        self.pixel_scale = imutils.get_pixel_scale(self.header)        
        # check to see if obj position is passed in - need to do this for off-center objects
        if (objra is not None): # unmask central elliptical region around object
            # get wcs from mask image
            wcs = WCS(self.header)
            
            # get x and y coord of galaxy from (RA,DEC) using mask wcs
            #print(f"\nobject RA={self.objra:.4f}, DEC={self.objdec:.4f}\n")
            self.xcenter,self.ycenter = wcs.wcs_world2pix(self.objra,self.objdec,0)
            self.xcenter_ra = self.xcenter
            self.ycenter_dec = self.ycenter            
            # convert sma to pixels using pixel scale from mask wcs
            self.pixel_scale = wcs.pixel_scale_matrix[1][1]
            #self.objsma_pixels = self.objsma/(self.pixel_scale*3600)
        
        # image 2 is designed to be the Halpha image, but it can be any second
        # image whereby you define the ellipse geometry using image 1, and
        # measure the photometry on image 1 and image 2
        #
        # self.image2_flag is True is image2 is provided
        if image2 is not None:
            self.image2_name = image2
            self.image2,self.header2 = fits.getdata(image2, header=True)
            self.image2_flag = True
        else:
            self.image2_flag = False
            self.image2 = None
            self.image2_name = None
            self.header2 = None
        self.image2_filter = image2_filter
        self.filter_ratio = filter_ratio
        # will use the gain to calculate the noise in the image
        try:
            self.gain = self.header['GAIN']
        except KeyError:
            print("WARNING: no GAIN keyword in header. Setting gain=1")
            self.gain = 1.
        self.psf = psf
        self.psf_ha = psf_ha

        if self.psf is not None:
            self.psf_data = fits.getdata(self.psf)
        if self.psf_ha is not None:
            self.hpsf_data = fits.getdata(self.psf)            
        # the mask should identify all pixels in the cutout image that are not
        # associated with the target galaxy
        # these will be ignored when defining the shape of the ellipse and when measuring the photometry
        #
        # self.mask_flag is True if a mask is provided
        if mask is not None:
            self.mask_image, self.mask_header = fits.getdata(mask,header=True)
            self.mask_flag = True
            # convert to boolean array with bad pixels = True
            self.boolmask = np.array(self.mask_image,'bool')
            self.masked_image = np.ma.array(self.image, mask = self.boolmask)
            if self.image2_flag:
                self.masked_image2 = np.ma.array(self.image2, mask = self.boolmask)
        else:
            print('not using a mask')
            self.mask_flag = False
            self.masked_image = self.image
            if self.image2_flag:
                self.masked_image2 = self.image2
        # image frame for plotting inside a gui
        # like if this is called from halphamain.py
        self.image_frame = image_frame

        # alternatively, for plotting with matplotlib
        # use this if running this code as the main program
        self.use_mpl = use_mpl
        self.napertures = napertures
        # assuming a typical fwhm 
        self.fwhm = 3.5
    def get_noise_in_aper(self, flux, area):
        ''' calculate the noise in an area '''
        if self.sky_noise is not None:
            noise_e = np.sqrt(flux*self.gain + area*self.sky_noise*self.gain)
            noise_adu = noise_e/self.gain
        else:
            noise_adu = np.nan
        return noise_adu

    def run_two_image_phot(self,write1=False):
        ''' 
        batch all of the functions that we run for the gui, including:

        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_sb()
        self.convert_units()
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()
        '''

        #print("detect objects")
        self.detect_objects()
        #print("find central")        
        self.find_central_object() 
        #print("find ellipse guess")               
        self.get_ellipse_guess()
        #print("measure phot")                
        self.measure_phot()
        #print("get M20")                
        #self.get_all_M20()
        #print("get frac masked pixels")                
        self.get_all_frac_masked_pixels()
        #print("calc sb")         
        self.calc_sb()
        #print("convert units")                
        #self.convert_units()
        #print("get asym")        
        ##self.get_image2_gini()
        #try:
        #    self.get_asymmetry()
        #except:
        #    print("WARNING: problem measuring asymmetry")
        #    self.asym = -99
        #    self.asym_err = -99
        #    self.asym2 = -99
        #    self.asym2_err = -99
        
        #print("running statmorph - please be patient...")
        #print()
        ##self.run_statmorph()
        ##self.statmorph_flag = True
        #try:    
        #    self.run_statmorph()
        #    self.statmorph_flag = True
        #except:
        #    self.statmorph_flag = False            
        #    print("WARNING: problem running statmorph")
        #print("writing phot fits tables")
        #self.write_phot_tables()
        self.write_phot_fits_table2_simple()
        #self.get_sky_noise()

        #print()
        #print("finished with photutils")
        #print()
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()
    
    def run_for_gui(self,runStatmorphFlag=True):
        ''' 
        batch all of the functions that we run for the gui, including:

        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_sb()
        self.convert_units()
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()
        '''

        print("detect objects")
        self.detect_objects()
        print("find central")        
        self.find_central_object() 
        print("find ellipse guess")               
        self.get_ellipse_guess()
        print("measure phot")                
        self.measure_phot()
        print("get M20")                
        self.get_all_M20()
        print("get frac masked pixels")                
        self.get_all_frac_masked_pixels()
        print("calc sb")         
        self.calc_sb()
        print("convert units")                
        self.convert_units()
        print("get asym")        
        #self.get_image2_gini()
        try:
            self.get_asymmetry()
        except:
            print("WARNING: problem measuring asymmetry")
            self.asym = -99
            self.asym_err = -99
            self.asym2 = -99
            self.asym2_err = -99
        #self.run_statmorph()
        #self.statmorph_flag = True
        if runStatmorphFlag:
            #print("running statmorph - please be patient...")
            #self.run_statmorph()
            #self.statmorph_flag = True
            
            try:
                print("running statmorph - please be patient...")
                print()
                self.run_statmorph()
                self.statmorph_flag = True
            except:
                self.statmorph_flag = False            
                print("WARNING: problem running statmorph")
            if self.image2 is not None:
                #self.run_statmorph_image2()
                #self.statmorph_flag2 = True                
                try:
                    print("running statmorph on image 2- please be patient...")
                    print()
                    self.run_statmorph_image2()
                    self.statmorph_flag2 = True
                except:
                    self.statmorph_flag2 = False            
                    print("WARNING: problem running statmorph on image 2")

            
        print("writing tables")
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()

        print()
        print("finished with photutils")
        print()
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()
    def run_with_galfit_ellipse(self, xc,yc,BA=1,THETA=0):
        '''
        replicating run_for_gui(), but taking input ellipse geometry from galfit

        '''
        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()

        ### RESET ELLIPSE GEOMETRY USING GALFIT VALUES
        self.xcenter = float(xc)
        self.ycenter = float(yc)
        self.position = (self.xcenter, self.ycenter)
        # leave sma as it is defined in get_ellipse_guess
        # so that we measure photometry over the same approximate region
        # in practice this could be different areas if the ellipticity is very different
        # self.sma = r

        # need to reset b to be consistent with galfit ellipticity
        BA = float(BA)
        THETA = float(THETA)
        #print('THETA inside phot wrapper',THETA, BA)
        self.b = BA*self.sma
        self.eps = 1 - BA
        #print(self.b,self.eps,self.sma,BA)
        t = THETA
        if t < 0:
            self.theta = np.radians(180. + t)
        else:
            self.theta = np.radians(t) # orientation in radians
        # EllipticalAperture gives rotation angle in radians from +x axis, CCW
        self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)

        ### AND NOW BACK TO OUR REGULAR PROGRAMMING
        #print('measuring phot')
        self.measure_phot()
        #print('measuring phot')
        self.calc_sb()
        #print('measuring converting units')
        self.convert_units()
        #print('writing table')
        #self.get_image2_gini()
        try:
            self.get_asymmetry()
        except:
            print("\nWARNING: problem calculating asymmetry, probably b/c image is rectangular...")
            self.asym = -99
            self.asym_err = -99
            self.asym2 = -99
            self.asym2_err = -99
            print()
        self.write_phot_fits_tables(prefix='GAL')
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()
    def detect_objects(self, snrcut=1.5,npixels=10):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''
        # this is not right, because the mask does not include the galaxy
        # updating based on photutils documentation
        # https://photutils.readthedocs.io/en/stable/background.html
        
        # get a rough background estimate

        # I already compute sky sigma and store it in header
        # should look for that and use that as a threshold if it's available

        try:
            
            skystd = self.header['SKYSTD']
            self.sky_noise = skystd
            self.sky = self.header['SKYMED']
        except KeyError:
            print("WARNING: SKYSTD not found in ",self.image_name)
            self.sky_noise = np.nan

        # get the value for halpha
        try:
            if self.header2 is not None:
                self.sky_noise2 = self.header2['SKYSTD']
                self.sky2 = self.header['SKYMED']
            else:
                print("WARNING: SKYSTD not found in ",self.image2_name)
                self.sky_noise2 = np.nan
                self.sky2 = np.nan
        except KeyError:
            print("WARNING: SKYSTD not found in ",self.image2_name)
            self.sky_noise2 = np.nan
            self.sky2 = np.nan
        
        if self.mask_flag:
            if self.sky_noise is not np.nan:
                self.threshold = self.sky_noise
            else:
                self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
            self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        else:
            if self.sky_noise is not np.nan:
                self.threshold = self.sky_noise
            else:
            
                self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            #self.cat = source_properties(self.image, self.segmentation)
            self.cat = SourceCatalog(self.image, self.segmentation)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation)
            


    def detect_objectsv2(self, snrcut=1.5,npixels=10):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''
        # this is not right, because the mask does not include the galaxy
        # updating based on photutils documentation
        # https://photutils.readthedocs.io/en/stable/background.html
        
        # get a rough background estimate
        sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        if self.mask_flag:
            threshold = detect_threshold(self.image, nsigma=snrcut,sigclip_sigma=3.0, mask=self.boolmask)
        else:
            threshold = detect_threshold(self.image, nsigma=snrcut,sigclip_sigma=3.0)

        segmentation = detect_sources(self.image, threshold, npixels=npixels)


        # make an object mask, expanding the area using a circular footprint
        #mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)        
        mask = make_source_mask(self.image,1.5,10,dilate_size=11)
        #plt.figure()
        #plt.imshow(mask,origin="lower")
        mean, median, std = sigma_clipped_stats(self.image, sigma=3.0, mask=mask)
        self.sky = mean
        self.sky_noise = std
        self.image -= self.sky

        if self.image2 is not None:
            # measure halpha properties using same segmentation image
            mean, median, std = sigma_clipped_stats(self.image2, sigma=3.0, mask=mask)            
            #self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
            self.sky2 = mean
            self.sky_noise2 = std

            # subtract sky
            self.image2 -= self.sky2
        
        
        # now make a new segmentation image based on the new noise estimate
        # default is 1.5*std
        #self.segmentation = detect_sources(self.image, snrcut*std, npixels=npixels)


        # subtract the sky, again...
        #print("\nsky value = ",self.sky)
        #self.image -= self.sky

        ##
        # old code
        ##


        if self.mask_flag:
            self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
            self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        else:
            self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            #self.cat = source_properties(self.image, self.segmentation)
            self.cat = SourceCatalog(self.image, self.segmentation)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        
        #self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
        #self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
        #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
        #self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
        
    def detect_objects_old(self, snrcut=1.5,npixels=11):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''

        if self.mask_flag:
            self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
            self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        else:
            self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            #self.cat = source_properties(self.image, self.segmentation)
            self.cat = SourceCatalog(self.image, self.segmentation)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
            
        # get average sky noise per pixel
        # threshold is the sky noise at the snrcut level, so need to divide by this
        self.sky_noise = np.mean(self.threshold)/snrcut


    def get_all_M20(self):
        # as a kludge, I am going to set all objects' M20 equal to this value
        # in the end, I will only keep the value for the central object...
        M20 = get_M20(self.cat,self.objectIndex)
        allM20 = M20*np.ones(len(self.cat))
        self.cat.add_extra_property('M20',allM20)
        self.M20_1 = M20

        # repeat for image2 if it's included
        if self.image2 is not None:
            M20 = get_M20(self.cat2,self.objectIndex)
            allM20 = M20*np.ones(len(self.cat))
            self.cat2.add_extra_property('M20',allM20)
            self.M20_2 = M20
            
    def get_all_frac_masked_pixels(self):
        # as a kludge, I am going to set all objects' masked fraction equal to this value
        # in the end, I will only keep the value for the central object...
        ntotal,nmasked,frac_masked = get_fraction_masked_pixels(self.cat,self.objectIndex)
        allfmasked = frac_masked*np.ones(len(self.cat))
        self.cat.add_extra_property('MASKEDFRAC',allfmasked)
        self.masked_fraction = allfmasked
        self.pixel_area = ntotal
        self.masked_pixel_area = ntotal - nmasked
    def get_sky_noise(self):
        '''
        * get the noise in image1 and image2 
        * noise is stored as SKYERR in image header
          - units of sky noise are erg/s/cm^2/arcsec^2
        '''

        # I already calculate this in detect_objects
        # get sky noise for image 1
        #if self.mask_flag:
        #    threshold = detect_threshold(self.image, nsigma=1,mask=self.boolmask)
        #else:
        #    threshold = detect_threshold(self.image, nsigma=snrcut)

        # add sky noise to image 1 header
        
        sky_noise_erg = self.sky_noise*self.uconversion1/self.pixel_scale**2

        print('r sky noise = ',sky_noise_erg)
        try:
            self.header.set('PHOT_SKY','{:.2f}'.format(self.sky),'sky in ADU')
        except AttributeError:
            print("Warning, self.sky not found, setting to zero")
            self.sky = 0
        self.header.set('SKYNOISE','{:.2f}'.format(self.sky_noise),'sky noise in ADU')        
        self.header.set('SKYERR','{:.2e}'.format(sky_noise_erg),'sky noise in erg/s/cm^2/arcsec^2')
        # save files
        fits.writeto(self.image_name,self.image,header=self.header,overwrite=True)
        self.im1_skynoise = sky_noise_erg
        # get sky noise for image 2
        if self.image2 is not None:
            sky_noise_erg2 = self.sky_noise2*self.uconversion2/self.pixel_scale**2
            self.header2.set('PHOT_SKY','{:.2f}'.format(self.sky2),'sky in ADU')
            self.header2.set('SKYNOISE','{:.2f}'.format(self.sky_noise2),'sky noise in ADU')        
            self.header2.set('SKYERR','{:.2e}'.format(sky_noise_erg2),'sky noise in erg/s/cm^2/arcsec^2')
            fits.writeto(self.image2_name,self.image2,header=self.header2,overwrite=True)
            try:
                sky_noise_erg2 = self.sky_noise2*self.uconversion2/self.pixel_scale**2
                self.header2.set('PHOT_SKY','{:.2f}'.format(self.sky2),'sky in ADU')
                self.header2.set('SKYNOISE','{:.2f}'.format(self.sky_noise2),'sky noise in ADU')        
                self.header2.set('SKYERR','{:.2e}'.format(sky_noise_erg2),'sky noise in erg/s/cm^2/arcsec^2')
                fits.writeto(self.image2_name,self.image2,header=self.header2,overwrite=True)
            
            except AttributeError:
                print("Warning, self.sky not found, setting to zero")
                self.sky2 = 0
                self.sky_noise2 = 0
                sky_noise_erg2 = 0
            

            self.im2_skynoise = sky_noise_erg2
            print('ha sky noise = ',sky_noise_erg2)

    def find_central_object(self):
        ''' 
        find the central object in the image and get its objid in segmentation image.
        object is stored as self.objectIndex
        '''

        # TODONE - need to be able to handle objects that are not at the center - should have option to pass in RA/DEC and then do like in maskwrapper
        if self.objra is not None:
            #print()
            print("getting object position from RA and DEC")
            #print()
            xc = self.xcenter_ra
            yc = self.ycenter_dec
        else:
            ydim,xdim = self.image.shape
            xc = xdim/2
            yc = ydim/2            
        distance = np.sqrt((np.ma.array(self.cat.xcentroid) - xc)**2 + (np.ma.array(self.cat.ycentroid) - yc)**2)        
        #distance = np.sqrt((self.cat.xcentroid.value - xdim/2.)**2 + (self.cat.ycentroid.value - ydim/2.)**2)
        # save object ID as the row in table with source that is closest to center

        # check to see if len(distance) is > 1

        if len(distance) > 1:
            try:
                self.objectIndex = np.arange(len(distance))[(distance == min(distance))][0]
            except IndexError:
                print("another $#@$# version change???",np.arange(len(distance))[(distance == min(distance))],len(distance))
                print('x vars: ',self.cat.xcentroid, xc)
                print('y vars: ', self.cat.ycentroid, yc)                
                print(self.cat)
                sys.exit()
        else:
            self.objectIndex = 0
            print("WARNING: only one object in the SourceCatalog!",distance)
        #print(self.objectIndex)
        if self.image2 is not None:
            # the object index in cat 2 is not necessarily the same
            # not sure if this is something I changed or if this has always been the case...
            
            distance = np.sqrt((np.ma.array(self.cat2.xcentroid) - xc)**2 + (np.ma.array(self.cat2.ycentroid) - yc)**2)        
            # save object ID as the row in table with source that is closest to center
            self.objectIndex2 = np.arange(len(distance))[(distance == min(distance))][0]
            
        if self.objra is not None:
            # check that distance of this object is not far from the original position
            xcat = self.cat.xcentroid[self.objectIndex]
            ycat = self.cat.ycentroid[self.objectIndex]

            offset = np.sqrt((xcat-self.xcenter_ra)**2 + (ycat-self.ycenter_dec)**2)
            if offset > 100:
                print()
                print("Hold the horses - something is not right!!!")
            

            #print("")            
            #print(f"comparing xcenter {xcat:.1f} and from ra {self.xcenter_ra:.1f}")
            #print(f"comparing ycenter {ycat:.1f} and from dec {self.ycenter_dec:.1f}")
            #print()

    def get_mask_from_segmentation(self):
        # create a copy of the segmentation image
        # replace the object index values with zeros        
        segmap = self.segmentation.data == self.cat.label[self.objectIndex]

        # subtract this from segmentation

        mask_data = self.segmentation.data - segmap*self.cat.label[self.objectIndex]
        # smooth 
        segmap_float = ndi.uniform_filter(np.float64(mask_data), size=10)
        mask = segmap_float > 0.5

        self.mask_image = mask
        self.boolmask = mask
        self.mask_flag = True


        # turn the segmentation image into a boolean mask

        
    def run_statmorph(self):
        '''
        run statmorph on image1 and image2 (if provided).

        results are stored as self.morph and self.morph2

        summary figures are save as XX statmorph-r.pdf and statmore-ha.pdf
        '''
        # show original
        #plt.figure()
        #plt.imshow(self.segmentation.data)
        # need segmentation map of object only
        segmap = self.segmentation.data == self.cat.label[self.objectIndex]
        segmap_float = ndi.uniform_filter(np.float64(segmap), size=10)
        segmap = segmap_float > 0.5
        self.segmap = segmap

        if self.mask_image is not None:
            mask = self.mask_image > 0
        else:
            mask = None
        #plt.figure()
        #plt.imshow(segmap, origin='lower', cmap='gray')

        # run statmorph on r-band image
        if self.psf is None:
            source_morphs = statmorph.source_morphology(self.image, segmap, gain=self.gain,mask=mask)
        else:
            source_morphs = statmorph.source_morphology(self.image, segmap, gain=self.gain, psf=self.psf_data,mask=mask)
        self.source_morphs = source_morphs
        self.morph = source_morphs[0]
        fig = make_figure(self.morph)
        figname = self.image_name.split('.fits')[0]
        fig.savefig(figname+'statmorph-r.pdf')
        self.segmap = segmap
    def run_statmorph_image2(self):

        if self.mask_image is not None:
            mask = self.mask_image > 0
        else:
            mask = None
        
        if self.psf_ha is None:
            source_morphs2 = statmorph.source_morphology(self.image2, self.segmap, gain=self.gain)
        else:
            source_morphs2 = statmorph.source_morphology(self.image2, self.segmap, gain=self.gain, psf=self.hpsf_data,mask=mask)
        self.source_morphs2 = source_morphs2            
        self.morph2 = source_morphs2[0]
        fig2 = make_figure(self.morph2)
        figname = self.image_name.split('.fits')[0]        
        fig2.savefig(figname+'statmorph-ha.pdf')
    def get_image2_gini(self, snrcut=1.5):
        ''' 
        calculate gini coefficient for image2 using pixels that are associated with r-band object ID

        this also calculates the sum and mag of the pixels associated with the central galaxy 
        (not sure why this is done together...)
        
        '''
        if self.mask_flag:
            self.threshold2 = detect_threshold(self.image2, nsigma=snrcut, mask=self.boolmask)
            self.segmentation2 = detect_sources(self.image2, self.threshold2, npixels=10,mask=self.boolmask)
            #self.cat2 = source_properties(self.image2, self.segmentation2, mask=self.boolmask)
            cat2 = SourceCatalog(self.image2, self.segmentation2, mask=self.boolmask)            
        else:
            self.threshold2 = detect_threshold(self.image2, nsigma=snrcut)
            self.segmentation2 = detect_sources(self.image2, self.threshold2, npixels=10)
            #self.cat2 = source_properties(self.image2, self.segmentation2)
            cat2 = SourceCatalog(self.image2, self.segmentation2)            

        '''
        select pixels associated with rband image in the segmentation
        AND
        pixels that are above the SNR cut in the Halpha image (image2)
        '''
        self.gini_pixels = (self.segmentation.data == self.cat.label[self.objectIndex]) & (self.segmentation2.data > 0.)

        #self.tbl = self.cat.to_table()
        self.gini2 = gini(self.image2[self.gini_pixels])
        #self.source_sum2 = np.sum(self.image2[self.gini_pixels])
        #self.source_sum2_erg = self.uconversion1*self.source_sum2
        #self.source_sum2_mag = self.magzp2 - 2.5*np.log10(self.source_sum2)
    def get_asymmetry(self):
        '''
        * goal is to measure the assymetry of the galaxy about its center
        * going to measure asymmetry from pixels in the segmentation image only, so

        '''
        # TODO - need to be able to handle images that are not square
        
        # for pixels in segmentation image of central object
        # (can't figure out a way to do this without looping
        # calculate delta_x and delta_y from centroid

        self.object_pixels = self.segmentation.data == self.cat.label[self.objectIndex]

        #xc = self.cat.xcentroid[self.objectIndex].value
        #yc = self.cat.ycentroid[self.objectIndex].value
        xc = self.cat.xcentroid[self.objectIndex]
        yc = self.cat.ycentroid[self.objectIndex]
        row,col = np.where(self.object_pixels)

        grid_size = 3
        sum_diff = np.zeros((grid_size,grid_size),'f')
        source_sum = np.zeros((grid_size,grid_size),'f')
        for dxc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
            for dyc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                drow = np.array((row-(yc+dyc)),'i')
                dcol = np.array((col-(xc+dxc)),'i')
                row2 = np.array(((yc+dyc) -1*drow),'i')
                col2 = np.array(((xc+dxc) -1*dcol),'i')
                sum_diff[dyc,dxc] = np.sum(np.abs(self.masked_image[row,col] - self.masked_image[row2,col2]))
                # divide by the sum of the original pixel values for object
                source_sum[dyc,dxc] = np.sum(self.image[self.object_pixels])
        asym = sum_diff/source_sum
        #print(asym)
        self.asym = np.min(asym)
        self.asym_err = np.std(asym)
        r,c = np.where(asym == np.min(asym))
        self.asym_center = np.array([r+yc,c+xc])

        print('asymmetry = {:.3f}+/-{:.3f}'.format(self.asym,self.asym_err))
        
        if self.image2_flag:
            print("getting asym for image2")
            # using the same segmentation image at for r-band
            # is this the correct thing to do?  does segmentation2.data need to be > 0?
            self.object_pixels2 = (self.segmentation.data == self.cat.label[self.objectIndex]) #& (self.segmentation2.data > 0.)

            #xc = self.cat.xcentroid[self.objectIndex].value
            #yc = self.cat.ycentroid[self.objectIndex].value

            # using the r-band centroid
            xc = self.cat.xcentroid[self.objectIndex]
            yc = self.cat.ycentroid[self.objectIndex]
            row,col = np.where(self.object_pixels2)
            sum_diff = np.zeros((grid_size,grid_size),'f')
            source_sum = np.zeros((grid_size,grid_size),'f')


            # looks like I am measuring asymmetry myself here?
            # doesn't statmorph do this? - looks like I am not using it
            for dxc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                for dyc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                    drow = np.array((row-(yc+dyc)),'i')
                    dcol = np.array((col-(xc+dxc)),'i')
                    row2 = np.array(((yc+dyc) -1*drow),'i')
                    col2 = np.array(((xc+dxc) -1*dcol),'i')
                    sum_diff[dyc,dxc] = np.sum(np.abs(self.masked_image2[row,col] - self.masked_image2[row2,col2]))
                    # divide by the sum of the original pixel values for object
                    source_sum[dyc,dxc] = np.sum(self.image2[self.object_pixels2])
            asym2 = sum_diff/source_sum
            
            #print('asym2 = ',asym2,r,c)
            # measure halpha asymmetry at pixel where R-band asymmetry is minimum
            try:
                self.asym2 = asym2[r,c][0]
                #print(self.asym2)
            except IndexError:
                try:
                    r,c = np.where(asym == np.min(asym2))
                
                    self.asym2 = asym2[r,c][0]
                except IndexError:
                    self.asym2 = np.nan
                    self.asym2_err = np.nan
                    self.asym2_center = np.nan
                    return
            self.asym2_err = np.std(asym2)
            r,c = np.where(asym == np.min(asym2))
            self.asym2_center = np.array([r+yc,c+xc])
            #print('asymmetry2 = ',self.asym2)
            print('asymmetry = {:.3f}+/-{:.3f}'.format(self.asym2,self.asym2_err))
            '''
            # use all the same images as for r-band measurement
            self.object_pixels2 = (self.segmentation.data == self.cat.id[self.objectIndex])# & (self.segmentation2.data > 0.)

            xc = self.cat.xcentroid[self.objectIndex].value
            yc = self.cat.ycentroid[self.objectIndex].value
            row,col = np.where(self.object_pixels2)

            drow = np.array((row-yc),'i')
            dcol = np.array((col-xc),'i')
            row2 = np.array((yc -1*drow),'i')
            col2 = np.array((xc -1*dcol),'i')
            sum_diff = np.sum(np.abs(self.masked_image2[row,col] - self.masked_image2[row2,col2]))
            # divide by the sum of the original pixel values for object
            source_sum = np.sum(self.image2[self.object_pixels2])
        
        
            self.asym2b = sum_diff/source_sum
            print('asymmetry2 = ',self.asym2b)
            '''
        
    def get_ellipse_guess(self, r=2.5):
        '''
        this gets the guess for the ellipse geometry from the detection catalog 
        '''
        obj = self.cat[self.objectIndex]
        #self.xcenter = obj.xcentroid.value
        #self.ycenter = obj.ycentroid.value

        if not self.fixcenter:
            self.xcenter = obj.xcentroid
            self.ycenter = obj.ycentroid


        #if self.objra is not None:
        #    print("")            
        #    print(f"comparing xcenter {self.xcenter:.1f} and from ra {self.xcenter_ra:.1f}")
        #    print(f"comparing ycenter {self.ycenter:.1f} and from dec {self.ycenter_dec:.1f}")
        #    print()
        self.position = (self.xcenter, self.ycenter)
        #print(self.position,self.xcenter,obj.xcentroid,self.ycenter,obj.ycentroid)
        self.sma = obj.semimajor_sigma.value * r
        self.start_size = self.sma
        self.b = obj.semiminor_sigma.value * r
        self.eps = 1 - self.b/self.sma
        self.gini = obj.gini
        self.source_sum = self.cat[self.objectIndex].segment_flux
        self.sky_centroid = obj.sky_centroid
        # orientation is angle in radians, CCW relative to +x axis
        t = obj.orientation.value
        #print('inside get_ellipse_guess, orientation = ',obj.orientation)
        if t < 0: # convert to positive angle wrt +x axis
            self.theta = np.pi+obj.orientation.to(u.rad).value
        else:
            self.theta = obj.orientation.to(u.rad).value # orientation in radians
        # EllipticalAperture gives rotation angle in radians from +x axis, CCW
        try:
            self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        except ValueError:
            print("\nTrouble in paradise...")
            print(self.position,self.sma,self.b,self.theta)
            sys.exit()
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)
    def draw_guess_ellipse(self):
        ''' DRAW INITIAL ELLIPSE ON R-BAND CUTOUT '''
        #
        markcolor='magenta'
        markwidth=1
        obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,self.sma, self.sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
        self.markhltag = self.image_frame.canvas.add(obj)
        self.image_frame.fitsimage.redraw()

    def draw_guess_ellipse_mpl(self):
        ''' DRAW INITIAL ELLIPSE ON R-BAND CUTOUT '''
        #
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        #plt.imshow(self.masked_image, cmap='Greys', norm=norm , origin='lower')
        display_image(self.masked_image)#, cmap='Greys', norm=norm , origin='lower')        
        plt.colorbar()
        self.aperture.plot(color='k', lw=1.)
        plt.show()

    def fit_ellipse(self):
        ''' FIT ELLIPSE '''
        #
        # create instance of photutils.Ellipse
        # https://photutils.readthedocs.io/en/stable/isophote.html
        self.ellipse = Ellipse(self.masked_image, self.guess)
        self.isolist = self.ellipse.fit_image()#sfix_pa = True, step=.5)#, fix_eps=True, fix_center=True)
        self.table = self.isolist.to_table()
        
    def draw_fit_results(self):
        ''' DRAW RESULTING FIT ON R-BAND CUTOUT '''
        markcolor='cyan'
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma), np.max(self.isolist.sma), 3)
            objlist = []
            for sma in smas:
                iso = self.isolist.get_closest(sma)
                obj = self.image_frame.dc.Ellipse(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps), rot_deg = np.degrees(iso.pa), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            self.markhltag = self.image_frame.canvas.add(self.coadd.dc.CompoundObject(*objlist))
            self.image_frame.fitsimage.redraw()
        else:
            print('problem fitting ellipse')
    def draw_fit_results_mpl(self):
        ''' draw fit results in matplotlib figure '''
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm, origin='lower')
        apertures = []
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma)+2, np.max(self.isolist.sma), 12)
            objlist = []
            for sma in smas:
                iso = self.isolist.get_closest(sma)
                #print(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps),  np.degrees(iso.pa))
                apertures.append(EllipticalAperture((iso.x0,iso.y0),iso.sma, iso.sma*(1-iso.eps), theta = np.degrees(iso.pa)))
            for aperture in apertures:
                aperture.plot(color='white',lw=1.5)
        plt.show()
        #plt.close()
    def show_seg_aperture(self,plotname=None):
        ''' matplotlib plotting to show apertures   '''
        tbl1 = self.cat.to_table()
        cat = self.cat
        r=3.
        apertures = []
        for obj in cat:
            position = np.transpose((obj.xcentroid, obj.ycentroid))
            try:
                a = obj.semimajor_axis_sigma.value * r
                b = obj.semiminor_axis_sigma.value * r
            except AttributeError:
                a = obj.semimajor_sigma.value * r
                b = obj.semiminor_sigma.value * r
                
            theta = obj.orientation.to(u.rad).value
            #print(theta)
            apertures.append(EllipticalAperture(position, a, b, theta=theta))
    
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
        ax1.imshow(self.masked_image, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        #cmap = segm_deblend.make_cmap(random_state=12345)
        ax2.imshow(self.segmentation.data, origin='lower')
        ax2.set_title('Segmentation Image')
        for aperture in apertures:
            aperture.plot(axes=ax1, color='white', lw=1.5)
            aperture.plot(axes=ax2, color='white', lw=1.5)
        #plt.show()
        if plotname is not None:
            plt.savefig(plotname)
        #plt.close()
    def measure_phot(self):
        '''
        # alternative is to use ellipse from detect
        # then create apertures and measure flux

        # rmax is max radius to measure ellipse
        # could cut this off based on SNR
        # or could cut this off based on enclosed flux?
        # or could cut off based on image dimension, and do the cutting afterward
        
        #rmax = 2.5*self.sma
        '''
        
        '''
        this is how becky set the apertures
        a = [0]
        for i in range(1,500):
        a.append(a[i-1] + hwhm + (hwhm*i*.1))
        
        '''
        # rmax is set according to the image dimensions
        # look for where the semi-major axis hits the edge of the image
        # could by on side (limited by x range) or on top/bottom (limited by y range)
        # 
        #print('xcenter, ycenter, theta = ',self.xcenter, self.ycenter,self.theta)
        rmax = np.min([(self.ximage_max - self.xcenter)/abs(np.cos(self.theta)),\
                       (self.yimage_max - self.ycenter)/abs(np.sin(self.theta))])
        #print('print rmax, ximage_max, image.shape = ',rmax,self.ximage_max,self.image.shape)
        '''
        this is how becky set the apertures
        a = [0]
        for i in range(1,500):
        a.append(a[i-1] + hwhm + (hwhm*i*.1))
        
        '''

        # TODO - update apertures to make use of input apertures
        index = np.arange(80)
        apertures = (index+1)*.5*self.fwhm*(1+(index+1)*.1)
        #apertures = (index+1)*self.fwhm*(1+(index+1)*.1)
        # cut off apertures at edge of image
        self.apertures_a = apertures[apertures < rmax]
        print(f"\nNumber of apertures = {len(index)}\n")
        #print('number of apertures = ',len(self.apertures_a))
        #self.apertures_a = np.linspace(3,rmax,40)
        self.apertures_b = (1.-self.eps)*self.apertures_a
        self.area = np.pi*self.apertures_a*self.apertures_b # area of each ellipse


        self.flux1 = np.zeros(len(self.apertures_a),'f')
        self.flux1_err = np.zeros(len(self.apertures_a),'f')
        if self.image2_flag:
            self.flux2 = np.zeros(len(self.apertures_a),'f')
            self.flux2_err = np.zeros(len(self.apertures_a),'f')
        self.allellipses = []
        for i in range(len(self.apertures_a)):
            #print('measure phot: ',self.xcenter, self.ycenter,self.apertures_a[i],self.apertures_b[i],self.theta)
            #,ai,bi,theta) for ai,bi in zip(a,b)]
            # EllipticalAperture takes rotation angle in radians, CCW from +x axis
            ap = EllipticalAperture((self.xcenter, self.ycenter),self.apertures_a[i],self.apertures_b[i],self.theta)#,ai,bi,theta) for ai,bi in zip(a,b)]
            self.allellipses.append(ap)

            if self.mask_flag:
                # check for nans, and add them to the mask
                nan_mask = self.image == np.nan
                combined_mask =  self.boolmask | nan_mask
                self.phot_table1 = aperture_photometry(self.image, ap, mask=combined_mask)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, mask=combined_mask)
            else:
                # subpixel is the method used by Source Extractor
                self.phot_table1 = aperture_photometry(self.image, ap, method = 'subpixel', subpixels=5)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, method = 'subpixel', subpixels=5)
            self.flux1[i] = self.phot_table1['aperture_sum'][0]
            
            # calculate noise
            self.flux1_err[i] = self.get_noise_in_aper(self.flux1[i], self.area[i])
            if self.image2_flag:
                self.flux2[i] = self.phot_table2['aperture_sum'][0]
                self.flux2_err[i] = self.get_noise_in_aper(self.flux2[i], self.area[i])
    def draw_phot_apertures(self,plotname=None):
        ''' matplotlib plotting to show apertures; provide a plot name to save the output figure   '''
        tbl1 = self.cat.to_table()
        cat = self.cat
        r=3.
        apertures = []
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6))
        clipped_data = sigma_clip(self.image,sigma_lower=5,sigma_upper=5)#,grow=10)
        norm = simple_norm(clipped_data, stretch='asinh',percent=99)
        
        #display_image(self.image)
        ax1.imshow(self.image, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        #cmap = segm_deblend.make_cmap(random_state=12345)
        ax2.imshow(self.segmentation.data, origin='lower')
        ax2.set_title('Segmentation Image')
        # plot a subset of apertures
        nap = len(self.allellipses)
        plotaps = np.array([0,nap//2,-1],'i')
        for i in plotaps:
            aperture = self.allellipses[i]
            aperture.plot(axes=ax1, color='c', lw=1.5)
            aperture.plot(axes=ax2, color='white', lw=1.5)
        if plotname is not None:
            plt.savefig(plotname)

    def calc_sb(self):
        # calculate surface brightness in each aperture

        # first aperture is calculated differently
        self.sb1 = np.zeros(len(self.apertures_a),'f')
        self.sb1_err = np.zeros(len(self.apertures_a),'f')

        self.sb1[0] = self.flux1[0]/self.area[0]
        self.sb1_err[0] = self.get_noise_in_aper(self.flux1[0], self.area[0])/self.area[0]
        # outer apertures need flux from inner aperture subtracted
        for i in range(1,len(self.area)):
            self.sb1[i] = (self.flux1[i] - self.flux1[i-1])/(self.area[i]-self.area[i-1])
            self.sb1_err[i] = self.get_noise_in_aper((self.flux1[i] - self.flux1[i-1]),(self.area[i]-self.area[i-1]))/(self.area[i]-self.area[i-1])

        # calculate SNR to follow Becky's method of cutting off analysis where SNR = 2
        self.sb1_snr = np.abs(self.sb1/self.sb1_err)
        # repeat for image 2 if it is provided
        if self.image2_flag:
            self.sb2 = np.zeros(len(self.apertures_a),'f')
            self.sb2_err = np.zeros(len(self.apertures_a),'f')
            self.sb2[0] = self.flux2[0]/self.area[0]
            self.sb2_err[0] = self.get_noise_in_aper(self.flux2[0], self.area[0])/self.area[0]
            for i in range(1,len(self.area)):
                self.sb2[i] = (self.flux2[i] - self.flux2[i-1])/(self.area[i]-self.area[i-1])
                self.sb2_err[i] = self.get_noise_in_aper((self.flux2[i] - self.flux2[i-1]),(self.area[i]-self.area[i-1]))/(self.area[i]-self.area[i-1])
            self.sb2_snr = np.abs(self.sb2/self.sb2_err)

    def convert_units(self):
        '''
        ###########################################################
        ### SET UP INITIAL PARAMETERS TO CALCULATE CONVERSION
        ### FROM ADU/S TO PHYSICAL UNITS
        ###########################################################
        '''

        self.pixel_scale = imutils.get_pixel_scale(self.header)
        try:
            self.magzp = float(self.header['PHOTZP'])

        except:
            print("WARNING: no PHOTZP keyword in image header. \nAssuming ZP=22.5")
            self.magzp = 22.5
        #print('mag zp = ',self.magzp)
        filter = self.header
        # multiply by bandwidth of filter to convert from Jy to erg/s/cm^2
        bandwidth1 = 3.e8*dwavelength['R']*1.e-10/(central_wavelength['R']*1.e-10)**2
        # need to figure out how to adjust automatically
        bandwidth1 = 3.e8*dwavelength['r']*1.e-10/(central_wavelength['r']*1.e-10)**2        
        self.uconversion1 = 3631.*10**(self.magzp/-2.5)*1.e-23*bandwidth1
        if self.image2_filter:
            bandwidth2 = 3.e8*dwavelength[self.image2_filter]*1.e-10/(central_wavelength[self.image2_filter]*1.e-10)**2
            try:
                self.magzp2 = float(self.header2['PHOTZP'])
                self.uconversion2 = 3631.*10**(self.magzp2/-2.5)*1.e-23*bandwidth2
            except:
                # use 25 as default ZP if none is provided in header
                self.uconversion2 = 3631.*10**(25/-2.5)*1.e-23*bandwidth2
        if self.filter_ratio is not None:
            if self.image2_flag:
                self.uconversion2b = self.filter_ratio*self.uconversion1
        else:
            self.uconversion2b = None
            
        ###########################################################
        ### CONVERT UNITS TO
        ### FLUX -> ERG/S/CM^2
        ### FLUX -> MAG
        ### SURFACE BRIGHTNESS -> ERG/S/CM^2/ARCSEC^2
        ### SURFACE BRIGHTNESS -> MAG/ARCSEC^2
        ###########################################################
        self.sky_noise_erg = self.sky_noise*self.uconversion1/self.pixel_scale**2
        self.flux1_erg = self.uconversion1*self.flux1
        self.flux1_err_erg = self.uconversion1*self.flux1_err
        self.source_sum_erg = self.uconversion1*self.source_sum
        self.source_sum_mag = self.magzp - 2.5*np.log10(self.source_sum)
        self.mag1 = self.magzp - 2.5*np.log10(self.flux1)
        self.mag1_err = self.mag1 - (self.magzp - 2.5*np.log10(self.flux1 + self.flux1_err))
        self.sb1_erg_sqarcsec = self.uconversion1*self.sb1/self.pixel_scale**2
        self.sb1_erg_sqarcsec_err = self.uconversion1*self.sb1_err/self.pixel_scale**2
        self.sb1_mag_sqarcsec = self.magzp - 2.5*np.log10(self.sb1/self.pixel_scale**2)
        self.sb1_mag_sqarcsec_err = self.sb1_mag_sqarcsec - (self.magzp - 2.5*np.log10((self.sb1 + self.sb1_err)/self.pixel_scale**2))
        if self.image2_flag:
            self.flux2_erg = self.uconversion2*self.flux2
            self.flux2_err_erg = self.uconversion2*self.flux2_err
            self.source_sum2 = self.cat2.segment_flux[self.objectIndex]
            self.source_sum2_erg = self.uconversion2*self.cat2.segment_flux[self.objectIndex]
            self.source_sum2_mag = self.magzp2 - 2.5*np.log10(self.source_sum)
            

            self.mag2 = self.magzp2 - 2.5*np.log10(self.flux2)
            self.mag2_err = self.mag2 - (self.magzp2 - 2.5*np.log10(self.flux2 + self.flux2_err))
            self.sb2_erg_sqarcsec = self.uconversion2*self.sb2/self.pixel_scale**2
            self.sb2_erg_sqarcsec_err = self.uconversion2*self.sb2_err/self.pixel_scale**2
            self.sb2_mag_sqarcsec = self.magzp2 - 2.5*np.log10(self.sb2/self.pixel_scale**2)
            # add error to flux and calculate magnitude, then take difference with original mag
            self.sb2_mag_sqarcsec_err = self.sb2_mag_sqarcsec - (self.magzp2 - 2.5*np.log10((self.sb2+self.sb2_err)/self.pixel_scale**2))
            # this next set uses the filter ratio and the filter 1 flux conversion to
            # convert narrow-band flux (filter 2) to physical units.
            if self.uconversion2b:
                conversion = self.uconversion2b
                self.flux2b_erg = conversion*self.flux2
                self.flux2b_err_erg = conversion*self.flux2_err
                self.sb2b_erg_sqarcsec = conversion*self.sb2/self.pixel_scale**2
                self.sb2b_erg_sqarcsec_err = conversion*self.sb2_err/self.pixel_scale**2
                self.sb2b_mag_sqarcsec = self.magzp2 - 2.5*np.log10(conversion*self.sb2/self.pixel_scale**2)
                self.sb2b_mag_sqarcsec_err = self.sb2b_mag_sqarcsec - (self.magzp2 - 2.5*np.log10(conversion*(self.sb2+self.sb2_err)/self.pixel_scale**2))
                
    def write_phot_tables(self):
        '''
        write out photometry for image and image2 in ascii format
        '''
        
        # radius enclosed flux
        outfile = open(self.image_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

        #outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
        #outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
        outfile.write('# radius flux flux_err sb sb_err sb_snr flux_erg flux_erg_err mag mag_err sb_ergsqarc sb_err_ergsqarc sb_magsqarc sb_err_magsqarc \n')
        for i in range(len(self.apertures_a)):
            s='%.2f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e '% \
                          (self.apertures_a[i],self.flux1[i],self.flux1_err[i],\
                           self.sb1[i], self.sb1_err[i], \
                          self.sb1_snr[i], \
                          self.flux1_erg[i], self.flux1_err_erg[i],\
                          self.mag1[i], self.mag1_err[i], \
                          self.sb1_erg_sqarcsec[i],self.sb1_erg_sqarcsec[i], \
                          self.sb1_mag_sqarcsec[i],self.sb1_mag_sqarcsec[i])
            s=s+'\n'
            outfile.write(s)
        outfile.close()

        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            outfile = open(self.image2_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles
    
            #outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
            #outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
            s = '# radius flux flux_err sb sb_err sb_snr flux_erg flux_erg_err mag mag_err sb_ergsqarc sb_err_ergsqarc sb_magsqarc sb_err_magsqarc'
            if self.uconversion2b:
                s = s +' fluxb_erg fluxb_erg_err sbb_ergsqarc sbb_err_ergsqarc sbb_magsqarc sbb_err_magsqarc'
            s = s + '\n'
            outfile.write(s)
            for i in range(len(self.apertures_a)):
                s = '%.2f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e '% \
                              (self.apertures_a[i],self.flux2[i],self.flux2_err[i],\
                               self.sb2[i], self.sb2_err[i],self.sb2_snr[i],\
                              self.flux2_erg[i], self.flux2_err_erg[i],\
                              self.mag2[i], self.mag2_err[i], \
                              self.sb2_erg_sqarcsec[i],self.sb2_erg_sqarcsec[i], \
                              self.sb2_mag_sqarcsec[i],self.sb2_mag_sqarcsec[i])
                if self.uconversion2b:
                    s=s+' %.3e %.3e %.3e %.3e %.3e %.3e'% \
                      (self.flux2b_erg[i], self.flux2b_err_erg[i],\
                      self.sb2b_erg_sqarcsec[i],self.sb2b_erg_sqarcsec[i], \
                      self.sb2b_mag_sqarcsec[i],self.sb2b_mag_sqarcsec[i])
                s = s+'\n'
                outfile.write(s)

            outfile.close()
    def write_phot_fits_tables(self, prefix=None):
        ''' write out photometry for image and image2 in fits format '''

        if prefix is None:
             outfile = self.image_name.split('.fits')[0]+'_phot.fits'
        else:
             outfile = self.image_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
        print('photometry outfile = ',outfile)

        data = [self.apertures_a*self.pixel_scale,self.apertures_a, \
             self.flux1,self.flux1_err,\
             self.sb1, self.sb1_err, \
             self.sb1_snr, \
             self.flux1_erg, self.flux1_err_erg,\
             self.mag1, self.mag1_err, \
             self.sb1_erg_sqarcsec,self.sb1_erg_sqarcsec_err, \
             self.sb1_mag_sqarcsec,self.sb1_mag_sqarcsec_err]

        names = ['sma_arcsec','sma_pix','flux','flux_err',\
                 'sb', 'sb_err', \
                 'sb_snr', \
                 'flux_erg', 'flux_erg_err',\
                 'mag', 'mag_err', \
                 'sb_erg_sqarcsec','sb_erg_sqarcsec_err', \
                 'sb_mag_sqarcsec','sb_mag_sqarcsec_err']

        units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '',\
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,\
                 u.mag,u.mag,\
                 u.erg/u.s/u.cm**2/u.arcsec**2,u.erg/u.s/u.cm**2/u.arcsec**2,\
                 u.mag/u.arcsec**2,u.mag/u.arcsec**2]


        #self.sky_noise,self.sky_noise_erg]
        #'sky_noise_ADU_sqpix','sky_noise_erg_sqarcsec']
        #u.adu/u.s/u.pixel**2,u.erg/u.s/u.cm**2/u.arcsec**2]        
        columns = []
        for i in range(len(data)):
            columns.append(Column(data[i],name=names[i],unit=units[i]))
        
        t = Table(columns)
        t.write(outfile, format='fits', overwrite=True)

        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            if prefix is None:
                outfile = self.image2_name.split('.fits')[0]+'_phot.fits'
            else:
                outfile = self.image2_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
    
            data = [self.apertures_a*self.pixel_scale,self.apertures_a, \
                self.flux2,self.flux2_err,\
                self.sb2, self.sb2_err, \
                self.sb2_snr, \
                self.flux2_erg, self.flux2_err_erg,\
                self.mag2, self.mag2_err, \
                self.sb2_erg_sqarcsec,self.sb2_erg_sqarcsec_err, \
                self.sb2_mag_sqarcsec,self.sb2_mag_sqarcsec_err]
            names = ['sma_arcsec','sma_pix','flux','flux_err',\
                'sb', 'sb_err', \
                'sb_snr', \
                'flux_erg', 'flux_erg_err',\
                'mag', 'mag_err', \
                'sb_erg_sqarcsec','sb_erg_sqarcsec_err', \
                'sb_mag_sqarcsec','sb_mag_sqarcsec_err']
            units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '',\
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,\
                 u.mag,u.mag,\
                 u.erg/u.s/u.cm**2/u.arcsec**2,u.erg/u.s/u.cm**2/u.arcsec**2,\
                 u.mag/u.arcsec**2,u.mag/u.arcsec**2]
            columns = []
            for i in range(len(data)):
                columns.append(Column(data[i],name=names[i],unit=units[i]))
            if self.uconversion2b:
                data = [self.flux2_erg, self.flux2_err_erg,\
                      self.sb2_erg_sqarcsec,self.sb2_erg_sqarcsec_err, \
                      self.sb2_mag_sqarcsec,self.sb2_mag_sqarcsec_err]
                names = ['flux2_erg', 'flux2_err_erg',\
                         'sb2_erg_sqarcsec','sb2_erg_sqarcsec_err', \
                         'sb2_mag_sqarcsec','sb2_mag_sqarcsec_err']
                units = [u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,\
                         u.erg/u.s/u.cm**2/u.arcsec**2,u.erg/u.s/u.cm**2/u.arcsec**2,\
                         u.mag/u.arcsec**2,u.mag/u.arcsec**2]
                for i in range(len(data)):
                    columns.append(Column(data[i],name=names[i],unit=units[i]))
            t = Table(columns)
            t.write(outfile, format='fits', overwrite=True)


    def write_phot_fits_table1_simple(self, prefix=None):
        ''' write out photometry for image and image2 in fits format '''

        if prefix is None:
             outfile = self.image_name.split('.fits')[0]+'_phot.fits'
        else:
             outfile = self.image_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
        print('photometry outfile = ',outfile)

        data = [self.apertures_a*self.pixel_scale,self.apertures_a, \
             self.flux1,self.flux1_err,\
             self.sb1, self.sb1_err, \
             self.sb1_snr]
             #self.flux1_erg, self.flux1_err_erg,\
             #self.mag1, self.mag1_err, \
             #self.sb1_erg_sqarcsec,self.sb1_erg_sqarcsec_err, \
             #self.sb1_mag_sqarcsec,self.sb1_mag_sqarcsec_err]

        names = ['sma_arcsec','sma_pix','flux','flux_err',\
                 'sb', 'sb_err', \
                 'sb_snr', ]

        units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '']


        #self.sky_noise,self.sky_noise_erg]
        #'sky_noise_ADU_sqpix','sky_noise_erg_sqarcsec']
        #u.adu/u.s/u.pixel**2,u.erg/u.s/u.cm**2/u.arcsec**2]        
        columns = []
        for i in range(len(data)):
            columns.append(Column(data[i],name=names[i],unit=units[i]))
        
        t = Table(columns)
        t.write(outfile, format='fits', overwrite=True)
    def write_phot_fits_table2_simple(self, prefix=None):
        """ write out phot for second image only - don't want to overwrite R phot """
        if prefix is None:
             outfile = self.image_name.split('.fits')[0]+'_phot.fits'
        else:
             outfile = self.image_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
        
        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            if prefix is None:
                outfile = self.image2_name.split('.fits')[0]+'_phot.fits'
            else:
                outfile = self.image2_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
    
            data = [self.apertures_a*self.pixel_scale*3600,self.apertures_a, \
                self.flux2,self.flux2_err,\
                self.sb2, self.sb2_err, \
                self.sb2_snr]
            names = ['sma_arcsec','sma_pix','flux','flux_err',\
                'sb', 'sb_err', \
                'sb_snr']
            units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '']
            columns = []
            for i in range(len(data)):
                columns.append(Column(data[i],name=names[i],unit=units[i]))
            self.tab2_simple = Table(columns)
            self.tab2_simple.write(outfile, format='fits', overwrite=True)


    def draw_phot_results(self):
        ''' DRAW RESULTING FIT ON R-BAND CUTOUT, for gui '''
        markcolor='cyan'
        objlist=[]
        markwidth=1.5
        for sma in self.apertures_a:
            obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,sma, sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
            objlist.append(obj)
            #print(self.xcenter,self.ycenter,sma, sma*(1-self.eps), self.theta, np.degrees(self.theta))
        self.markhltag = self.image_frame.canvas.add(self.image_frame.dc.CompoundObject(*objlist))
        self.image_frame.fitsimage.redraw()
        
    def draw_phot_results_mpl(self):
        ''' draw results in matplotlib figure '''
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm, origin='lower')

        apertures = []
        for sma in self.apertures_a:
            apertures.append(EllipticalAperture((self.xcenter,self.ycenter),sma, sma*(1-self.eps), theta = self.theta))
            
        for aperture in apertures:
            aperture.plot(color='white',lw=1.5)
        plt.show()
    def plot_profiles(self):
        ''' enclosed flux and surface brightness profiles, save figure '''
        plt.close("all")        
        plt.figure(figsize=(10,4))
        plt.subplots_adjust(wspace=.3)
        plt.subplot(2,2,1)
        #plt.plot(self.apertures_a,self.flux1,'bo')
        plt.errorbar(self.apertures_a,self.flux1,self.flux1_err,fmt='b.')
        plt.title('R-band')
        #plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Enclosed flux')
        plt.gca().set_yscale('log')
        if self.image2_flag:
            plt.subplot(2,2,2)
            plt.errorbar(self.apertures_a,self.flux2,self.flux2_err,fmt='b.')
            #plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Enclosed flux')
            plt.title('H-alpha')
            plt.gca().set_yscale('log')
        # plot surface brightness vs radius
        plt.subplot(2,2,3)
        #plt.plot(self.apertures_a,self.flux1,'bo')
        plt.errorbar(self.apertures_a,self.sb1,self.sb1_err,fmt='b.')
        plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Surface Brightess')
        plt.gca().set_yscale('log')
        if self.image2_flag:
            plt.subplot(2,2,4)
            plt.errorbar(self.apertures_a,self.sb2,self.sb2_err,fmt='b.')
            plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Surface Brightness')
            plt.gca().set_yscale('log')
        #plt.show()
        plt.savefig(self.image_name.split('.fits')[0]+'-enclosed-flux.png')
    def plot_fancy_profiles(self):
        # plot enclosed flux        
        fig = plt.figure(figsize=(10,4))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95,wspace=.3)

        labels = ['R','Halphax100']
        alphas = [1,.4,.6,.4]
        x = self.apertures_a*self.pixel_scale
        fluxes = [self.flux1_erg,self.flux2_erg]
        flux_errs = [self.flux1_err_erg,self.flux2_err_erg]
        plt.subplot(1,2,1)
        for i,t in enumerate(fluxes):
            y0 = fluxes[i]            
            y1 = y0+flux_errs[i]
            y2 = y0-flux_errs[i]

            if (i == 1) + (i == 3):
                y0=y0*100
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(x,y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(x,y0,'-',lw=2,color=mycolors[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend(loc='lower right')

        plt.subplot(1,2,2)

        fluxes = [self.sb1_erg_sqarcsec,self.sb2_erg_sqarcsec]
        flux_errs = [self.sb1_erg_sqarcsec_err,self.sb2_erg_sqarcsec_err]
        for i,t in enumerate(fluxes):
            y0 = fluxes[i]            
            y1 = y0+flux_errs[i]
            y2 = y0-flux_errs[i]

            if (i == 1) + (i == 3):
                y0=y0*100
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(x,y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(x,y0,'-',lw=2,color=mycolors[i])
        
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Surface Brightness',fontsize=16)

        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')        
        plt.legend(loc='upper right')
            
        plt.savefig(self.image_name.split('.fits')[0]+'_enclosed_flux_fancy.png')        
        #plt.close(fig)
        
        
if __name__ == '__main__':
    image = 'MKW8-18216-R.fits'
    mask = 'MKW8-18216-R-mask.fits'
    image2 = 'MKW8-18216-CS.fits'
    nsaid='18045'
    prefix = 'MKW8-'
    nsaid='110430'
    nsaid='157146'
    prefix = 'NRGs27-'
    image = prefix+nsaid+'-R.fits'
    mask = prefix+nsaid+'-R-mask.fits'
    image2 = prefix+nsaid+'-CS.fits'
    # testing on 2017 pointing 1
    # second galaxy has clear halpha but profile is not fit
    # want to make sure we record some size
    image = 'v17p01-N119230-A742747-R.fits'
    rphot_table = 'v17p01-N119230-A742747-R_phot.fits'
    image2 = 'v17p01-N119230-A742747-CS.fits'
    haphot_table = 'v17p01-N119230-A742747-CS_phot.fits'
    mask = 'v17p01-N119230-A742747-R-mask.fits'
    image = 'v17p01-N118647-A8219-R.fits'
    rphot_table = 'v17p01-N118647-A8219-R_phot.fits'
    image2 = 'v17p01-N118647-A8219-CS.fits'
    haphot_table = 'v17p01-N118647-A8219-CS_phot.fits'
    mask = 'v17p01-N118647-A8219-R-mask.fits'
    myfilter = '4'
    myratio = .0406
    
    # testing on 2019 pointing 1
    # second galaxy has clear halpha but profile is not fit
    # want to make sure we record some size
    image = 'VFID3623-CGCG118-019-v19p001-R.fits'
    
    rphot_table = 'VFID3623-CGCG118-019-v19p001-R-phot.fits'
    image2 = 'VFID3623-CGCG118-019-v19p001-CS.fits'
    haphot_table = 'VFID3623-CGCG118-019-v19p001-CS-phot.fits'
    mask = 'VFID3623-CGCG118-019-v19p001-R-mask.fits'
    
    myfilter = 'inthalpha'
    myratio = .0356
    #image = 'MKW8-18037-R.fits'
    #mask = 'MKW8-18037-R-mask.fits'
    #image = 'r-18045-R.fits'
    #mask = 'r-18045-R-mask.fits'

    prefix = 'VFID0501-UGC09556-BOK-20210315-VFID0501'
    prefix = 'VFID2772-NGC2964-HDI-20180313-p019'    
    image = prefix+'-R.fits'
    rphot_table = prefix+'-R-phot.fits'
    image2 = prefix+'-CS.fits'
    haphot_table = prefix+'-CS-phot.fits'
    mask = prefix+'-R-mask.fits'
    myfilter = '4'
    myratio = .0497427

    try:
        e = ellipse(image,mask=mask, image2=image2, use_mpl=True,image2_filter=myfilter, filter_ratio=myratio)
    except FileNotFoundError:
        print("so sorry, but no images were loaded")
        print("try e = ellipse(imagename) to start")
    ## print('detect objects')
    ## e.detect_objects()
    ## print('find central')
    ## e.find_central_object()
    ## print('get guess')
    ## e.get_ellipse_guess()
    ## print('draw guess')
    ## e.draw_guess_ellipse_mpl()
    ## print('fit ellipse')
    ## e.fit_ellipse()
    ## print('plot results')
    ## e.draw_fit_results_mpl()

    print("--- %s seconds ---" % (time.time() - start_time))


