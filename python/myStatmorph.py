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

    @lazyproperty
    def gini(self):
        '''overwriting function so that it uses the reg segmap'''
        #print("gini before new calc = ",self.gini)
        print("calculating gini2\n")

        self._segmap_gini = np.array(self._segmap.data == 1)
        image = self._cutout_stamp_maskzeroed.flatten()
        segmap = self._segmap_gini.flatten()

        sorted_pixelvals = np.sort(np.abs(image[segmap]))
        n = len(sorted_pixelvals)
        if n <= 1 or np.sum(sorted_pixelvals) == 0:
            warnings.warn('[gini] Not enough data for Gini calculation.',
                          AstropyUserWarning)
            self.flag = 2
            return -99.0  # invalid

        indices = np.arange(1, n+1)  # start at i=1
        gini = (np.sum((2*indices-n-1) * sorted_pixelvals) /
                (float(n-1) * np.sum(sorted_pixelvals)))
        #print(gini)
        return gini
