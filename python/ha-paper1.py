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
#print(homedir,os.path.join(homedir,'/github/Virgo/programs/'))
sys.path.append(homedir+'/github/Virgo/programs/')
import virgoCommon

from scipy.stats import spearmanr
from scipy.stats import median_absolute_deviation as MAD

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
        # get mass weighted center
        badMassFlag = self.magphys['logMstar_med'] < 2

        combinedMass = self.magphys['logMstar_med'] * ~badMassFlag + self.magphys['logMstar_best'] * badMassFlag

        self.combinedMass = combinedMass
        mass_flag = combinedMass > 7.8

        # morph flag?


        # inclination flag?

        self.ba_flag = self.halpha['GAL_BA'] < 0.85
        # fraction of unmasked pixels within segmentation > 50%

        self.ephot_flag = (self.halpha['M26'] > 10)& (self.halpha['M24'] <20)
        self.ephot_flag =  (self.halpha['M24'] > 0)          
        self.rad_flag = self.halpha['R24'] > 0
        self.galfit_flag = (self.halpha['GAL_RE'] < 50) & (self.halpha['GAL_MAG'] < 20)  & (self.halpha['GAL_N'] < 8) 
        self.morph_flag = (self.halpha['ELLIP_ASYM'] > -90)& (self.halpha['C30'] < 1) 
        self.smorph_flag = (self.halpha['SMORPH_XCENTROID'] > 0) & (self.halpha['SMORPH_S'] > -.1)


        self.ephot_hflag = (self.halpha['HM16'] > 0)& (self.halpha['HM17'] >0) 
        self.rad_hflag = (self.halpha['HR16'] > 0) & (self.halpha['HR17'] < 200)
        self.morph_hflag = (self.halpha['HC30'] < 1) &(self.halpha['ELLIP_HASYM'] > -10) & \
            (self.halpha['ELLIP_HASYM'] < 5) &  (self.halpha['ELLIP_HGINI'] > 0)  &\
            (self.halpha['ELLIP_HGINI']< 1) 
        self.smorph_hflag = (self.halpha['SMORPH_HXCENTROID'] > 0) & (self.halpha['SMORPH_S'] > -.1) &\
            (self.halpha['SMORPH_HGINI'] > -90) &(self.halpha['SMORPH_HM20'] > -90) &\
            (self.halpha['SMORPH_HS'] > -.5)

        
        self.sampleflag = snr_flag  & mass_flag #& ba_flag

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
        self.rdetect_flag = self.halpha['M24'] > 0        

        pass
    def plot_fwhm_r_h(self):
        x = self.halpha['R_FWHM']
        y = self.halpha['H_FWHM'] - x


        flag = self.main['HAobsflag'] & (x > 0)

        std = np.std(y[flag])
        ave = np.mean(y[flag])
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
        plt.axhline(y=ave,ls='-',color='k')

        plt.axhline(y=ave+std,ls='--',color='k')
        plt.axhline(y=ave-std,ls='--',color='k')        
        print(f"H FWHM - R FWHM = {np.mean(y):.2f} +/- {std:.2f}")
        plt.savefig(plotdir+"R-H-FWHM.pdf")
        plt.savefig(plotdir+"R-H-FWHM.png")

    def compare_m24_ephot(self):
        """compare r-band mag 24 with JM's ephot  """
        #plt.figure(figsize=(8,10))
        plt.close("all")
        fig,axs = plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(6,8),squeeze=True)

        axfontsize=16

        plt.subplots_adjust(hspace=.05)
        sampleflag = self.rdetect_flag & (self.ephot['FLUX_SB24_R'] > 0)
        # color by the SMA_SB24
        color = np.log10(self.ephot['SMA_SB24'][sampleflag])
        v1 = 2
        v2 = 100
        gr = 2.5*np.log10(self.ephot['FLUX_SB24_G']) - 2.5*np.log10(self.ephot['FLUX_SB24_R'])
        #color = gr[self.rdetect_flag]
        
        # define magnitude
        nsamag = 22.5 - 2.5*np.log10(self.ephot['FLUX_SB24_R'])
        dy = self.halpha['M24'][sampleflag] - nsamag[sampleflag]

        self.mag_offset_jm = dy
        # plot 1:1
        axs[0].scatter(nsamag[sampleflag],self.halpha['M24'][sampleflag],c=color,s=20)#,vmin=v1,vmax=v2)#,vmin=2,vmax=100)
        xl = np.linspace(10,18,100)
        axs[0].plot(xl,xl,'r--')

        axs[0].set_ylabel('photutils mag 24',fontsize=axfontsize)
        #ax[0].set_xlabel('NSA v0 r-band mag',fontsize=20)

        #axs[0].set_axis([])

        #################################
        # plot residuals
        #################################


        sp = axs[1].scatter(nsamag[sampleflag],dy,c=color,s=20)#,vmin=2,vmax=100)

        allax=[]
        for a in axs:
            plt.sca(a)
            allax.append(plt.gca())

        cb = plt.colorbar(sp,fraction=.08,ax=allax)
        #cb.set_label(label='g-r JM',fontsize=14)
        cb.set_label(label='log10(SMA_24)',fontsize=14)
        #axs[1].axhline(y=0,ls='-',c='k')
        axs[1].set_xlabel('JM mag 24',fontsize=axfontsize)
        axs[1].set_ylabel('photutils - JM',fontsize=axfontsize)


        ###################################
        # get std and mean offset
        ###################################
        statflag = (dy < 99) & (dy > -99)
        magmean = np.nanmean(dy[statflag])
        magstd = np.nanstd(dy[statflag])
        magmedian = np.nanmedian(dy[statflag])
        magmad = MAD(dy[statflag])
        plt.axhline(y=magmean,ls='-',c='r',label='mean offset')
        plt.axhline(y=magmean+magstd,ls='--',c='r',label='mean+/-std')
        plt.axhline(y=magmean-magstd,ls='--',c='r')

        plt.axhline(y=magmedian,ls='-',c='c',label='median offset')        
        plt.axhline(y=magmedian+magmad,ls=':',c='c',label='med+/-mad')
        plt.axhline(y=magmedian-magmad,ls=':',c='c')
        plt.text(0.05,0.9,f"mean/std={magmean:.2f}+/-{magstd:.2f}",fontsize=12,transform=plt.gca().transAxes,horizontalalignment='left')
        plt.text(0.05,0.8,f"median/MAD={magmedian:.2f}+/-{magmad:.2f}",fontsize=12,transform=plt.gca().transAxes,horizontalalignment='left')
        #plt.ylim(-1,1)
        plt.legend()

        # save outliers

    def compare_ephot_jm_geometry(self):
        """compare ephot ellipse parameters with JM parameters """

        # index, Re, PA, BA
        jm_labels = ['SMA_SB24','BA_MOMENT','PA_MOMENT']
        jm_params = [self.ephot[i] for i in jm_labels]
        labels = ['R24','ELLIP_EPS','ELLIP_THETA']
        ephot_params = [self.halpha[i] for i in ['R24','ELLIP_EPS','ELLIP_THETA']]

        color = 22.5 - 2.5*np.log10(self.ephot['FLUX_SB24_R'])
        v1 = 10
        v2 = 21
        clabel = 'JM M24'

        color = self.ephot['BA_MOMENT']
        v1 = 0
        v2 = 1
        clabel = 'JM B/A'

        
        flag = self.rdetect_flag 
        plt.figure(figsize=(10,8))
        plt.subplots_adjust(wspace=.5,hspace=.4)
        allax = []
        for i in range(len(ephot_params)):
            plt.subplot(2,2,i+1)
            if i == 0:
                # make y axis a fractional difference
                dy = ephot_params[i][flag] - jm_params[i][flag]
                dy = dy/jm_params[i][flag]
            elif i == 1: # convert ephot EPS to B/A, EPS = 1-B/A

                dy = (1-ephot_params[i][flag]) - jm_params[i][flag]
                plt.scatter(jm_params[i][flag],dy,label=labels[i],alpha=.5,c=color[flag],vmin=v1,vmax=v2)

                # save output to inspect outliers

                self.ba_offset_jm = dy
            elif i == 2: # make definition of  PA consistent
                # photutils and JM use ref axis that differs by 90 deg
                
                x = jm_params[i][flag]
                above90 = x > 90
                x[above90] = x[above90]-180
                dy = ephot_params[i][flag] - x - 90
                #dy = ephot_params[i][flag] #- jm_params[i][flag]                
                #plt.scatter(x,dy,label=labels[i],alpha=.5,c=color[flag],vmin=v1,vmax=v2)
            else:
                dy = ephot_params[i][flag] - jm_params[i][flag]
            plt.scatter(jm_params[i][flag],dy,label=labels[i],alpha=.7,c=color[flag],vmin=v1,vmax=v2)

            ###################################
            # get std and mean offset
            ###################################
            statflag = (dy < 99) & (dy > -99)
            magmean = np.nanmean(dy[statflag])
            magstd = np.nanstd(dy[statflag])
            magmedian = np.nanmedian(dy[statflag])
            magmad = MAD(dy[statflag])
            plt.axhline(y=magmean,ls='-',c='r',label='mean offset')
            plt.axhline(y=magmean+magstd,ls='--',c='r',label='mean+/-std')
            plt.axhline(y=magmean-magstd,ls='--',c='r')

            plt.axhline(y=magmedian,ls='-',c='c',label='median offset')        
            plt.axhline(y=magmedian+magmad,ls=':',c='c',label='med+/-mad')
            plt.axhline(y=magmedian-magmad,ls=':',c='c')
            plt.text(0.05,0.9,f"mean/std={magmean:.2f}+/-{magstd:.2f}",fontsize=12,transform=plt.gca().transAxes,horizontalalignment='left')
            plt.text(0.05,0.8,f"median/MAD={magmedian:.2f}+/-{magmad:.2f}",fontsize=12,transform=plt.gca().transAxes,horizontalalignment='left')

            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')
            #if i == 0:
            #    plt.ylim(-5,5)
            #elif i == 1:
            #    plt.axis([0,100,-100,100])
            #plt.axhline(ls='--',c='k')
            if i == 0:
                plt.ylabel("(Photutils - JM)/JM",fontsize=14)
            else:
                plt.ylabel("Photutils - JM",fontsize=14)                
            plt.xlabel("JM "+jm_labels[i].replace("GAL_",""),fontsize=14)
            if i ==0:
                plt.gca().set_xscale('log')
            
            allax.append(plt.gca())
            if i == 2:
                plt.legend(bbox_to_anchor=(1.5,1))

        cb = plt.colorbar(ax=allax,fraction=0.08,label='JM M24')
        cb.set_label(label=clabel,fontsize=14)
        
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
            
    def plot_qc_galfit_residuals(self):
        """compare galfit parameters with NSA """

        # index, Re, PA, BA
        nsa_params = [self.nsav0[i] for i in ['SERSIC_N','SERSIC_TH50','SERSIC_BA','SERSIC_PHI']]
        labels = ['GAL_N','GAL_RE','GAL_BA','GAL_PA']
        galfit_params = [self.halpha[i] for i in ['GAL_N','GAL_RE','GAL_BA','GAL_PA']]
        flag = self.detect_flag & self.main['NSAV0flag']
        plt.figure(figsize=(8,6))
        plt.subplots_adjust(wspace=.2,hspace=.4)
        for i in range(len(galfit_params)):
            plt.subplot(2,2,i+1)
            if i == 3: # make definition of galfit and NSA PA consistent
                # NSA PA
                x = nsa_params[i][flag]
                above90 = x > 90
                x[above90] = x[above90]-180
                dy = galfit_params[i][flag] - x
                plt.scatter(x,dy,label=labels[i],alpha=.5)
            else:
                dy = galfit_params[i][flag] - nsa_params[i][flag]
                plt.scatter(nsa_params[i][flag],dy,label=labels[i],alpha=.5)

            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')
            if i == 0:
                plt.ylim(-5,5)
            elif i == 1:
                plt.axis([0,100,-100,100])
            plt.axhline(ls='--',c='k')
            if i%2 == 0:
                plt.ylabel("GALFIT - NSA",fontsize=14)
            plt.xlabel("NSA "+labels[i].replace("GAL_",""),fontsize=14)
            
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
        yvars = ['ELLIP_HGINI','ELLIP_HM20','HC30','ELLIP_HASYM']        

        cvars = [self.halpha['ELLIP_ASYM'],self.halpha['ELLIP_ASYM'],self.magphys['logMstar'],self.magphys['logMstar']]
        labels = ['ELLIP_ASYM','ELLIP_ASYM','logMstar','logMstar']        
        flag = self.sampleflag & self.morph_flag & self.morph_hflag
        vmins = [.1,.1,8.5,8.5]
        vmaxs = [.8,.8,10.5,10.5]
        allax = []
        for i in range(len(xvars)):
            plt.subplot(2,2,i+1)
            #plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=self.magphys['logsSFR'][flag],alpha=.6,vmin=-12,vmax=-9)
            plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=self.magphys['logMstar'][flag],alpha=.6,vmin=8,vmax=10.5)            
            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')                
            plt.xlabel(xvars[i],fontsize=14)
            plt.ylabel(yvars[i],fontsize=14)
            x1,x2 = plt.xlim()
            xline = np.linspace(x1,x2,100)
            plt.plot(xline,xline,'k--')
            #if i == 2:
            #    plt.ylim(0,1.2)

            #plt.legend()
            allax.append(plt.gca())
        #plt.colorbar(label='magphys log sSFR',ax=allax,fraction=.08)
        plt.colorbar(label='magphys log Mstar',ax=allax,fraction=.08)
    def get_asym_gals(self,hasym=2.,rasym=.45):
        # print out ids for galaxies with high asymmetry in halpha but med in r
        flag = (self.halpha['ELLIP_HASYM'] > hasym) & (self.halpha['ELLIP_ASYM'] > 1)

        flag = flag & self.sampleflag & self.morph_flag & self.morph_hflag
        flag = flag & (self.magphys['logMstar'] > 8.) & (self.magphys['logsSFR'] > -10.5)
        print("VFIDs for gals with high Halpha asym and regular r-band asym")
        print(self.main['VFID'][flag])
    def plot_qc_qmorph_smorph_r(self):
        """
        plot r-band quantitative morphology parameters

        C30 vs Gini

        Gini vs M20 - color code by asym

        C30 vs GAL_N

        PETRO_CON vs GAL_N
        """
        sflag = self.halpha['SMORPH_XCENTROID'] > 0
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(wspace=.3,hspace=.3)

        xvars = ['ELLIP_GINI','ELLIP_M20','C30','ELLIP_ASYM']
        yvars = ['SMORPH_GINI','SMORPH_M20','SMORPH_C','SMORPH_A']        

        cvars = [self.halpha['ELLIP_ASYM'],self.halpha['ELLIP_ASYM'],self.magphys['logMstar'],self.magphys['logMstar']]
        labels = ['ELLIP_ASYM','ELLIP_ASYM','logMstar','logMstar']        
        flag = self.sampleflag & sflag
        vmins = [.1,.1,8.5,8.5]
        vmaxs = [.8,.8,10.5,10.5]
        allax = []
        color = 22.5 -2.5*np.log10(self.ephot['FLUX_SB24_R'])
        color = self.magphys['logsSFR']        
        for i in range(len(xvars)):
            plt.subplot(2,2,i+1)
            plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=color[flag],alpha=.6,vmin=-12,vmax=-9)
            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')                
            plt.xlabel(xvars[i],fontsize=14)
            plt.ylabel(yvars[i],fontsize=14)
            x1,x2 = plt.xlim()
            xline = np.linspace(x1,x2,100)
            
            if i == 2:
                plt.plot(xline,5*xline,'k--')
            else:
                plt.plot(xline,xline,'k--')
            #plt.legend()
            allax.append(plt.gca())
        plt.colorbar(label='magphys log sSFR',ax=allax,fraction=.08)

    def plot_qc_qmorph_smorph_ha(self):
        """
        plot r-band quantitative morphology parameters

        C30 vs Gini

        Gini vs M20 - color code by asym

        C30 vs GAL_N

        PETRO_CON vs GAL_N
        """
        sflag = (self.halpha['SMORPH_XCENTROID'] > 0) & (self.halpha['SMORPH_HXCENTROID'] > 0) \
            & (self.halpha['SMORPH_HM20'] > -2)& (self.halpha['ELLIP_HASYM'] > -2)& (self.halpha['ELLIP_HASYM'] < 2)
        print("Number with good statmorph flag = ",np.sum(sflag))
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(wspace=.3,hspace=.3)

        xvars = ['ELLIP_HGINI','ELLIP_HM20','HC30','ELLIP_HASYM']
        yvars = ['SMORPH_HGINI','SMORPH_HM20','SMORPH_HC','SMORPH_HA']        

        cvars = [self.halpha['ELLIP_ASYM'],self.halpha['ELLIP_ASYM'],self.magphys['logMstar'],self.magphys['logMstar']]
        labels = ['ELLIP_ASYM','ELLIP_ASYM','logMstar','logMstar']        
        flag = self.sampleflag & sflag
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
            
            if i == 2:
                plt.plot(xline,5*xline,'k--')
            else:
                plt.plot(xline,xline,'k--')
            plt.legend()
            allax.append(plt.gca())
        plt.colorbar(label='magphys log sSFR',ax=allax,fraction=.08)

    def plot_smorph_ha_r(self):
        """
        plot r-band quantitative morphology parameters

        C30 vs Gini

        Gini vs M20 - color code by asym

        C30 vs GAL_N

        PETRO_CON vs GAL_N
        """
        
        plt.figure(figsize=(10,6))
        plt.subplots_adjust(wspace=.5,hspace=.5)

        xvars = ['SMORPH_GINI','SMORPH_M20','SMORPH_C','SMORPH_A','SMORPH_S']
        yvars = ['SMORPH_HGINI','SMORPH_HM20','SMORPH_HC','SMORPH_HA','SMORPH_HS']        

        cvars = [self.halpha['ELLIP_ASYM'],self.halpha['ELLIP_ASYM'],self.magphys['logMstar_med'],self.magphys['logMstar_med']]
        labels = ['ELLIP_ASYM','ELLIP_ASYM','logMstar','logMstar']        
        flag = self.sampleflag & self.smorph_flag & self.smorph_hflag
        vmins = [.1,.1,8.5,8.5]
        vmaxs = [.8,.8,10.5,10.5]
        allax = []
        for i in range(len(xvars)):
            plt.subplot(2,3,i+1)
            #plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=self.magphys['logsSFR'][flag],alpha=.6,vmin=-12,vmax=-9)
            plt.scatter(self.halpha[xvars[i]][flag],self.halpha[yvars[i]][flag],c=self.magphys['logMstar_med'][flag],alpha=.6,vmin=8,vmax=10.5)            
            #if i < 2:
            #    plt.gca().set_xscale('log')
            #    plt.gca().set_yscale('log')                
            plt.xlabel(xvars[i],fontsize=14)
            plt.ylabel(yvars[i],fontsize=14)
            x1,x2 = plt.xlim()
            xline = np.linspace(x1,x2,100)
            plt.plot(xline,xline,'k--')
            #if i == 2:
            #    plt.ylim(0,1.2)

            #plt.legend()
            allax.append(plt.gca())
        #plt.colorbar(label='magphys log sSFR',ax=allax,fraction=.08)
        plt.colorbar(label='magphys log Mstar',ax=allax,fraction=.08)

    def plot_delta_gini_ssfr(self):
        plt.figure(figsize=(8,6))

        y = self.halpha['ELLIP_HGINI']-self.halpha['ELLIP_GINI']
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
        y = self.halpha['ELLIP_HGINI']-self.halpha['ELLIP_GINI']
        x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR']
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag &self.morph_flag & self.morph_hflag  & (self.magphys['logMstar'] > 0) #&  (self.magphys['logsSFR'] > -11)
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
        x = self.halpha['ELLIP_HGINI']-self.halpha['ELLIP_GINI']
        #x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR']
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.magphys['logMstar'] > 0) & self.morph_flag & self.morph_hflag #& (self.magphys['logsSFR'] > -11)
        plt.scatter(x[flag],y[flag],c=c[flag])
        cb = plt.colorbar()
        cb.set_label(label="$sSFR$",size=16)
        rho,p = spearmanr(x[flag],y[flag])

        plt.ylabel(r"$\Delta \ Asym \ (H\alpha - R)$",fontsize=16)
        plt.xlabel(r"$\Delta \ Gini \ (H\alpha - R)$",fontsize=16)        
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")

    def plot_delta_asym_m20(self,usephot=True,printids=False,asymcut=0.5,asymmax=10):
        plt.figure(figsize=(8,6))
        if usephot:
            y = self.halpha['ELLIP_HASYM']-self.halpha['ELLIP_ASYM']
            x = self.halpha['ELLIP_HGINI']-self.halpha['ELLIP_GINI']
            morph_flag = self.morph_flag & self.morph_hflag
        else:
            y = self.halpha['SMORPH_HA']-self.halpha['SMORPH_A']
            x = self.halpha['SMORPH_HGINI']-self.halpha['SMORPH_GINI']
            morph_flag = self.smorph_flag & self.smorph_hflag
            
        #x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR']
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.magphys['logMstar'] > 0) & morph_flag #& (self.magphys['logsSFR'] > -11)
        plt.scatter(x[flag],y[flag],c=c[flag])
        cb = plt.colorbar()
        cb.set_label(label="$sSFR$",size=16)
        rho,p = spearmanr(x[flag],y[flag])

        plt.ylabel(r"$\Delta \ Asym \ (H\alpha - R)$",fontsize=16)
        plt.xlabel(r"$\Delta \ M_{20}\ (H\alpha - R)$",fontsize=16)        
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")
        if printids:
            print(self.main['VFID'][(y> asymcut) & (y < asymmax) & flag])
    def plot_delta_asym_gini(self):
        plt.figure(figsize=(8,6))

        y = self.halpha['ELLIP_HASYM']-self.halpha['ELLIP_ASYM']
        x = self.halpha['ELLIP_HGINI']-self.halpha['ELLIP_GINI']
        #x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR_med']
        #c = np.log10(self.env['n5th_2D'])
        flag = self.sampleflag   & (self.combinedMass > 0) & self.morph_flag & self.morph_hflag #& (self.magphys['logsSFR'] > -11)
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
        x = self.halpha['ELLIP_HGINI']-self.halpha['ELLIP_GINI']
        #x = np.log10(self.env['n5th_2D'])
        c = self.magphys['logsSFR_med']
        c = self.halpha['ELLIP_HASYM']-self.halpha['ELLIP_ASYM']        
        #c = np.log10(self.env['n5th_2D'])
        #flag = self.sampleflag   & (self.magphys['logMstar_med'] > 0) & ~x.mask #& (self.magphys['logsSFR'] > -11)
        flag = self.sampleflag   & (self.combinedMass > 0)  & self.morph_flag & self.morph_hflag # & ~x.mask #& (self.magphys['logsSFR'] > -11)
        plt.scatter(x[flag],y[flag],c=c[flag])#,vmin=-.2,vmax=3,alpha=1)
        cb = plt.colorbar()
        #cb.set_label(label="$sSFR$",size=16)
        cb.set_label(label="$\Delta \ Asym $",size=16)        
        rho,p = spearmanr(x[flag],y[flag])

        plt.ylabel(r"$\Delta \ C30 \ (H\alpha - R)$",fontsize=16)
        plt.xlabel(r"$\Delta \ Gini \ (H\alpha - R)$",fontsize=16)        
        s = f"Spearman rank: \nrho={rho:.2f}, pvalue={p:.1e}"
        plt.text(0.05,0.95,s,transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment="top")
        #plt.xlim(-1,1)
        #plt.ylim(-.3,.75)
        #return x[flag],y[flag]
    def plot_ha_sizes(self):
        """ compare various size measurements for halpha """

        tab = Table(self.halpha['HR16','HR17','HR_F25','HR_F50','HR_F75','HPETRO_R','HPETRO_R50','HPETRO_R90'])
        tab.add_column(self.environment,name='environment')

        # make table into a pandas dataframe

        df = pdf(data=np.array(tab))

        sns.pairplot(df[self.sampleflag],hue="environment",corner=True)#, diag_kind="hist")        

        # use seaborn to make pair plot!

if __name__ == "__main__":
    # read in v2 tables
    h = haplots(tabledir='/home/rfinn/research/Virgo/tables-north/v2/',tableprefix='vf_v2_')
    h.read_all()
    h.get_detect_flag()
    h.get_sample_flag()
    h.get_environment()
