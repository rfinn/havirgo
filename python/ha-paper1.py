#!/usr/bin/env python

from astropy.table import Table
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
import os

import sys
homedir = os.getenv("HOME")
import virgoCommon

sys.path.append(os.path.join(homedir,'/github/Virgo/programs/'))

from readtablesv2 import vtables

class haplots(vtables):
    def get_detect_flag(self):
        self.detect_flag = self.halpha['ELLIP_HSUM'] > 0

        pass

    def plot_qc_galfit(self):
        """compare galfit parameters with NSA """

        # index, Re, PA, BA
        nsa_params = [self.nsav0[i] for i in ['SERSIC_N','SERSIC_TH50','SERSIC_BA','SERSIC_PHI']]
        labels = ['GAL_N','GAL_RE','GAL_BA','GAL_PA']
        galfit_params = [self.halpha[i] for i in ['GAL_N','GAL_RE','GAL_BA','GAL_PA']]
        flag = self.detect_flag & self.main['NSAV0flag']
        plt.figure(figsize=(8,6))
        plt.subplots_adjust(wspace=.3,hspace=.3)
        for i in range(len(galfit_params)):
            plt.subplot(2,2,i+1)
            plt.scatter(nsa_params[i][flag],galfit_params[i][flag],label=labels[i],alpha=.5)
            if i < 2:
                plt.gca().set_xscale('log')
                plt.gca().set_yscale('log')                
            plt.xlabel("NSA "+labels[i].replace("GAL_",""),fontsize=14)
            plt.ylabel("GALFIT "+labels[i],fontsize=14)
            plt.legend()
            
    def plot_Ttype_logmstar(self):
        """  test"""
        pass

if __name__ == "__main__":
    # read in v2 tables
    h = haplots(tabledir='/home/rfinn/research/Virgo/tables-north/v2/',tableprefix='vf_v2_')
    h.read_all()
    h.get_detect_flag()
