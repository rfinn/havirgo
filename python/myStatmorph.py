#!/usr/bin/env python



import numpy as np
#import matplotlib.pyplot as plt
#from astropy.visualization import simple_norm
#from astropy.modeling.models import Sersic2D
#from astropy.convolution import convolve, Gaussian2DKernel
from astropy.utils import lazyproperty
#from photutils.segmentation import detect_threshold, detect_sources
#from photutils.morphology import gini
#import time

import statmorph
#import os
#from astropy.io import fits

class myStatmorph(statmorph.SourceMorphology):

    """
    add on to statmorph 

    * changed to use the same segmentation map for the gini coefficient calculation 
      as it does for the other morph calculations

    """
    
    @lazyproperty
    def _segmap_gini(self):
        '''overwriting function so that it uses the reg segmap'''
        #self._image[self._slice_stamp]        
        segmap = np.array(self._segmap.data== 1,'i')
        return segmap[self._slice_stamp]
    

    def print_quantities(self):
        print("Measured Quantities")
        for q in statmorph.statmorph._quantity_names:
            print(f"{q}")
