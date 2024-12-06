#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from astropy.modeling.models import Sersic2D
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.utils import lazyproperty
from photutils.segmentation import detect_threshold, detect_sources
from photutils.morphology import gini
import time
import statmorph
import os
from astropy.io import fits
from statmorph.utils.image_diagnostics import make_figure

import myStatmorph

def get_median_sky(image1):
    threshold = detect_threshold(image1, 1.)
    #print(threshold)

    npixels = 15  # minimum number of connected pixels

    # create a segmentation map
    # this will contain all objects with npixels above threshold
    convolved_image = convolve(image1, psfsky)
    segmap_sky = detect_sources(convolved_image, threshold, npixels)
    
    skyregion = np.ma.array(image1,mask=segmap_sky.data > 0)
    mean,median,std = sigma_clipped_stats(skyregion,sigma=3.0,cenfunc=np.ma.median)
    #print(f"median sky = {median}")
    return median

class galaxy():

    def __init__(self,image1,image2,mask=None):
        # read in the images and headers
        self.image1_filename = image1
        self.image2_filename = image2        
        self.image1, self.header1 = fits.getdata(image1,header=True)
        self.image2, self.header2 = fits.getdata(image2,header=True)

        if mask is not None:
            mask = fits.getdata(mask)
            self.mask_bool = mask <1
        else:
            self.mask_bool = np.zeros_like(self.image1,'bool')

        self.image1 = np.ma.array(self.image1,mask=self.mask_bool)
        self.image2 = np.ma.array(self.image2,mask=self.mask_bool)        
    def subtract_background(self):
        median = get_median_sky(self.image1)
        image1_minus_bkg = self.image1 - median
        pass

    @lazyproperty
    def segmap(self):
        # psf for making segmentation image for galaxy
        kernel = Gaussian2DKernel(9)
        kernel.normalize()  # make sure kernel adds up to 1
        psf = kernel.array  # we only need the numpy array
        #plt.imshow(psf, origin='lower', cmap='gray')


        threshold = detect_threshold(self.image1, 2)
        #print(threshold)

        npixels = 15  # minimum number of connected pixels

        # create a segmentation map
        # this will contain all objects with npixels above threshold
        convolved_image = convolve(self.image1, psf)
        segmap = detect_sources(convolved_image, threshold, npixels)

        # estimate the position of the galaxy as the center of the x and y dimensions
        # this assumes that the object is centered in the image
        ymax, xmax = self.image1.shape
        #print(xmax,ymax)
        xc = xmax //2
        yc = ymax //2
        #print(f"xcenter = {xc}, ycenter={yc}")

        # get value of segmentation image at center
        print(f"segmentation value at center of image = {segmap.data[yc,xc]}")
        objid = segmap.data[yc,xc]


        # create a new image with ones where segmap == objid and zeros otherwise
        newsegmap = np.array(segmap.data == objid,'i')
        plt.figure(figsize=(8,6))
        plt.imshow(newsegmap,origin="lower")
        plt.title("segmentation map")
        plt.colorbar()
        return newsegmap

    def run_statmorph(self):
        self.morph1 = myStatmorph.myStatmorph(self.image1,self.segmap,1,gain=1,cutout_extent=1.5)
        self.morph2 = myStatmorph.myStatmorph(self.image2,self.segmap,1,gain=1,cutout_extent=1.5)        

    def plot_statmorph(self):
        fig = make_figure(self.morph1)
        im2name = os.path.basename(self.image1_filename)
        outfile = im2name.replace('.fits','-statmorph.pdf')
        plt.savefig(outfile)


        fig = make_figure(self.morph2)
        im2name = os.path.basename(self.image2_filename)
        outfile = im2name.replace('.fits','-statmorph.pdf')
        plt.savefig(outfile)
        

    

    
