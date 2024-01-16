#!/usr/bin/env python

import os
import sys

homedir = os.getenv("HOME")

# add virgo respository
sys.path.append(os.path.join(homedir),'github/Virgo/programs/')

from readtablesv2 import vtables


class haplots(vtables):


    def compare_sfrs(self):
        """get halpha-corrected SFRs vs SED sfrs """

        pass

    
