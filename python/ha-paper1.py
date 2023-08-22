#!/usr/bin/env python

from astropy.table import Table

import seaborn as sns
from pandas import DataFrame as pdf

from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
import os

import sys
homedir = os.getenv("HOME")
sys.path.append(os.path.join(homedir,'/github/Virgo/programs/'))
import virgoCommon

from scipy.stats import spearmanr

from readtablesv2 import vtables

plotdir = homedir+"/research/Virgo/plots/halpha/"
############################
## FUNCTIONS
############################

def scatter_hist(x, y, ax, ax_histx, ax_histy,nbins=20,alpha=1):
    """https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html """
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, alpha=0.5)

    # now determine nice limits by hand:
    binwidth = (np.max(x)-np.min(x))/nbins
    xbins = np.linspace(np.min(x), np.max(x) + binwidth,nbins)
    ax_histx.hist(x, bins=xbins)
    binwidth = (np.max(y)-np.min(y))/nbins    
    ybins = np.linspace(np.min(y), np.max(y) + binwidth,nbins)    
    ax_histy.hist(y, bins=ybins, orientation='horizontal')


class haplots(vtables):
    def get_sample_flag(self):
        # SNR in halpha
        snr_flag = np.abs(self.halpha['HF_TOT']/self.halpha['HF_TOT_ERR']) > 3
        snr_flag = (self.halpha['HM16'] > 10)
        # stellar mass cut
        mass_flag = self.magphys['logMstar'] > 8.

        # morph flag?


        # inclination flag?

        ba_flag = self.halpha['GAL_BA'] < 0.85

        # fraction of unmasked pixels within segmentation > 50%

        
        self.sampleflag = snr_flag  & mass_flag & ba_flag

    def get_environment(self):
        """  goal is to get environment variable - cluster, filament, field """
        env = []
        for i in range(len(self.env)):
            if self.env['cluster_member'][i]:
                env.append("cluster")
            elif self.env['pure_field'][i]:
                env.append("field")
            elif self.env["filament_member"][i]:
                env.append("filament")
            #elif self.env['rich_group_memb'][i]:
            #    env.append("rich_group")
            #elif self.env['poor_group_memb'][i]:
            #    env.append("poor_group")
            else:
                env.append("other")
        self.environment = env
        
        
    def correct_halpha_flux(self):
        self.hacorr = self.halpha['FILT_COR']
    def get_detect_flag(self):
        self.detect_flag = self.halpha['ELLIP_HSUM'] > 0

        pass
    def plot_fwhm_r_h(self):
        x = self.halpha['R_FWHM']
        y = self.halpha['H_FWHM'] - x

        std = np.std(y)
        flag = self.main['HAobsflag'] & (x > 0)

        x = x[flag]
        y = y[flag]        
        # FROM https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html
        # Create a Figure, which doesn't have to be square.
        fig = plt.figure()#layout='constrained')
        # Create the main axes, leaving 25% of the figure space at the top and on the
        # right to position marginals.
        ax = fig.add_gridspec(left=.2,top=0.75, right=0.75).subplots()
        # The main axes' aspect can be fixed.
        ax.set(aspect=1)
        # Create marginal axes, which have 25% of the size of the main axes.  Note that
        # the inset axes are positioned *outside* (on the right and the top) of the
        # main axes, by specifying axes coordinates greater than 1.  Axes coordinates
        # less than 0 would likewise specify positions on the left and the bottom of
        # the main axes.
        ax_histx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
        ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
        # Draw the scatter plot and marginals.
        scatter_hist(x, y, ax, ax_histx, ax_histy)
        ax.set_xlabel("R-band FWHM (arcsec)",fontsize=16)
        ax.set_ylabel(r"$\rm H\alpha - R \ FWHM \ (arcsec)$ ",fontsize=16)


        
        #plt.plot(x[flag],y[flag],'bo',alpha=.7)

        # plot reference lines
        plt.sca(ax)
        plt.axhline(y=0,ls='-',color='k')

        plt.axhline(y=std,ls='--',color='k')
        plt.axhline(y=-1*std,ls='--',color='k')        
        print(f"H FWHM - R FWHM = {np.mean(y):.2f} +/- {std:.2f}")
        plt.savefig(plotdir+"R-H-FWHM.pdf")
        plt.savefig(plotdir+"R-H-FWHM.png")        
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

    def plot_qc_qmorph(self):
        """
        plot r-band quantitative morphology parameters

        C30 vs Gini

        Gini vs M20 - color code by asym

        C30 vs GAL_N

        PETRO_CON vs GAL_N
        """
        
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(wspace=.3,hspace=.3)

        xvars = ['ELLIP_GINI','ELLIP_M20','GAL_N','ELLIP_GINI']
        yvars = ['C30','ELLIP_GINI','C30','ELLIP_ASYM']
        cvars = [self.halpha['ELLIP_ASYM'],self.halpha['ELLIP_ASYM'],self.magphys['logMstar'],self.magphys['logMstar']]
        labels = ['ELLIP_ASYM','ELLIP_ASYM','logMstar','logMstar']        
        flag = self.detect_flag
        vmins = [.1,.1,8.5,8.5]
        vmaxs = [.8,.8,10.5,10.5]        
        for i in range(len(xvars)):
            plt.subplot(2,2,i+1)
            plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=cvars[i][flag],alpha=.5,vmin=vmins[i],vmax=vmaxs[i])
            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')                
            plt.xlabel(xvars[i],fontsize=14)
            plt.ylabel(yvars[i],fontsize=14)
            plt.colorbar(label=labels[i])
            plt.legend()
    def plot_qc_qmorph_ha_r(self):
        """
        plot r-band quantitative morphology parameters

        C30 vs Gini

        Gini vs M20 - color code by asym

        C30 vs GAL_N

        PETRO_CON vs GAL_N
        """
        
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(wspace=.3,hspace=.3)

        xvars = ['ELLIP_GINI','ELLIP_M20','C30','ELLIP_ASYM']
        yvars = ['ELLIP_GINI2','ELLIP_HM20','HC30','ELLIP_HASYM']        

        cvars = [self.halpha['ELLIP_ASYM'],self.halpha['ELLIP_ASYM'],self.magphys['logMstar'],self.magphys['logMstar']]
        labels = ['ELLIP_ASYM','ELLIP_ASYM','logMstar','logMstar']        
        flag = self.sampleflag
        vmins = [.1,.1,8.5,8.5]
        vmaxs = [.8,.8,10.5,10.5]
        allax = []
        for i in range(len(xvars)):
            plt.subplot(2,2,i+1)
            plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=self.magphys['logsSFR'][flag],alpha=.6,vmin=-12,vmax=-9)
            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')                
            plt.xlabel(xvars[i],fontsize=14)
            plt.ylabel(yvars[i],fontsize=14)
            x1,x2 = plt.xlim()
            xline = np.linspace(x1,x2,100)
            plt.plot(xline,xline,'k--')
            if i == 2:
                plt.ylim(0,1.2)

            plt.legend()
            allax.append(plt.gca())
        plt.colorbar(label='magphys log sSFR',ax=allax,fraction=.08)

    def plot_delta_gini_ssfr(self):
        plt.figure(figsize=(8,6))

        y = self.halpha['ELLIP_GINI2']-self.halpha['ELLIP_GINI']
        x = self.magphys['logsSFR']
        c = self.magphys['logSFR']
        c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.magphys['logMstar'] > 0) & ~y.mask #& (self.magphys['logsSFR'] > -11.5)
        plt.scatter(x[flag],y[flag],c=c[flag],vmin=-1,vmax=1.5)
        cb = plt.colorbar()
        cb.set_label(label="$\Sigma_5$",size=16)
        rho,p = spearmanr(x[flag],y[flag])

        plt.xlabel(r"$sSFR \ (M_\odot/yr/M_\odot)$",fontsize=16)
        plt.ylabel(r"$\Delta \ Gini \ (H\alpha - R)$",fontsize=16)
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")
        #return x[flag],y[flag]
    def plot_delta_asym_sigma5(self):
        plt.figure(figsize=(8,6))

        y = self.halpha['ELLIP_HASYM']-self.halpha['ELLIP_ASYM']
        y = self.halpha['ELLIP_GINI2']-self.halpha['ELLIP_GINI']
        x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR']
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.magphys['logMstar'] > 0) & ~y.mask & (self.magphys['logsSFR'] > -11)
        plt.scatter(x[flag],y[flag],c=c[flag])
        cb = plt.colorbar()
        cb.set_label(label="$sSFR$",size=16)
        rho,p = spearmanr(x[flag],y[flag])

        plt.xlabel(r"$\Sigma_5$",fontsize=16)
        plt.ylabel(r"$\Delta \ Asym \ (H\alpha - R)$",fontsize=16)
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")
        #return x[flag],y[flag]    

    def plot_delta_asym_gini(self):
        plt.figure(figsize=(8,6))

        y = self.halpha['ELLIP_HASYM']-self.halpha['ELLIP_ASYM']
        x = self.halpha['ELLIP_GINI2']-self.halpha['ELLIP_GINI']
        #x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR']
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.magphys['logMstar'] > 0) & ~y.mask #& (self.magphys['logsSFR'] > -11)
        plt.scatter(x[flag],y[flag],c=c[flag])
        cb = plt.colorbar()
        cb.set_label(label="$sSFR$",size=16)
        rho,p = spearmanr(x[flag],y[flag])

        plt.ylabel(r"$\Delta \ Asym \ (H\alpha - R)$",fontsize=16)
        plt.xlabel(r"$\Delta \ Gini \ (H\alpha - R)$",fontsize=16)        
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")

    def plot_delta_c30_gini(self):
        plt.figure(figsize=(8,6))

        y = self.halpha['HC30']-self.halpha['C30']
        x = self.halpha['ELLIP_GINI2']-self.halpha['ELLIP_GINI']
        #x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR']
        c = self.halpha['ELLIP_HASYM']-self.halpha['ELLIP_ASYM']        
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.magphys['logMstar'] > 0) & ~x.mask #& (self.magphys['logsSFR'] > -11)
        plt.scatter(x[flag],y[flag],c=c[flag],vmin=-.2,vmax=.7)
        cb = plt.colorbar()
        #cb.set_label(label="$sSFR$",size=16)
        cb.set_label(label="$\Delta \ Asym $",size=16)        
        rho,p = spearmanr(x[flag],y[flag])

        plt.ylabel(r"$\Delta \ C30 \ (H\alpha - R)$",fontsize=16)
        plt.xlabel(r"$\Delta \ Gini \ (H\alpha - R)$",fontsize=16)        
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")    
        #return x[flag],y[flag]
    def plot_ha_sizes(self):
        """ compare various size measurements for halpha """

        tab = Table(self.halpha['HR16','HR17','HR_F25','HR_F50','HR_F75','HPETRO_R','HPETRO_R50','HPETRO_R90'])
        tab.add_column(self.environment,name='environment')

        # make table into a pandas dataframe

        df = pdf(data=tab)

        # use seaborn to make pair plot!

if __name__ == "__main__":
    # read in v2 tables
    h = haplots(tabledir='/home/rfinn/research/Virgo/tables-north/v2/',tableprefix='vf_v2_')
    h.read_all()
    h.get_detect_flag()
    h.get_sample_flag()
    h.get_environment()
