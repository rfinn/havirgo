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

class myStatmorph(statmorph.SourceMorphology):

    def testprop(self):
        pass
    
    @lazyproperty
    def _segmap_gini(self):
        '''overwriting function so that it uses the reg segmap'''

        #self._image[self._slice_stamp]        
        segmap = np.array(self._segmap.data== 1,'i')
        return segmap[self._slice_stamp]
    

