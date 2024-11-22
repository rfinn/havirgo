#!/usr/bin/env

"""
* following on from duplicates.py
* adding cross match with vf_v2 catalogs


"""

from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.table import Table, Column
import os
import numpy as np
from scipy.stats import median_abs_deviation as MAD
from scipy.stats import ks_2samp, spearmanr
import seaborn as sb
import sys
import collections

homedir = os.getenv("HOME")
sys.path.append(os.path.join(homedir,'github/Virgo/programs/'))
from readtablesv2 import vtables
sys.path.append(os.path.join(homedir,'github/havirgo/python/'))
from duplicates import duplicates

mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

class hanalysis(duplicates):
    def test(self):
        print('hi')

    def get_haflag(self):
        haflag = np.abs(v.halpha['HF_TOT']/v.halpha['HF_TOT_ERR']) > 5
    def cut_vf_tables(self):
        """cut the vf tables to match the halpha table """
        pass
    def plot_sfr_mstar(self,hsmorphmaxflag=3,rsmorphmaxflag=2):

        plt.figure(figsize=(8,6))        
        # plot add vf galaxies

        self.get_smorph_flag(hsmorphmaxflag=hsmorphmaxflag,rsmorphmaxflag=rsmorphmaxflag)
        haflag = d.smorphflag
        x = self.v.magphys['logMstar_med']
        y = self.v.magphys['logSFR_med']
        flag1 = (x > 7)
        hflag1 = self.v.magphys['logMstar_med'][d.htab['VFINDEX']] > 5
                                                    
        flag1b = flag1 &  ~self.v.halpha['HAflag'] & ~self.v.main['COflag']# valid magphys fit        
        plt.plot(x[flag1b],y[flag1b],'k.',alpha=.1,label='All VF')        

        # halpha sample
        
        flag3 = flag1 & self.v.main['COflag']
        flag4 = hflag1 & self.smorphflag

        


        flag2 = hflag1 & ~self.smorphflag
        #flag2 = flag1 & self.v.halpha['HAflag'] & ~haflag

        plt.plot(x[self.htab['VFINDEX'][flag2]],y[self.htab['VFINDEX'][flag2]],'bv',c=mycolors[3],mec='k',alpha=.7,label=f'Halpha ~STAMORPH ({np.sum(flag2)})',markersize=9)        
        plt.plot(x[self.htab['VFINDEX'][flag4]],y[self.htab['VFINDEX'][flag4]],'bs',c=mycolors[0],mec='k',alpha=.8,label=f'Halpha STATMORPH ({np.sum(flag4)})',markersize=9)

        #plt.plot(x[flag3],y[flag3],'bo',c=mycolors[1],alpha=.8,label='Primary CO Sample',markersize=5)

        plt.xlabel('$\log(M_\star/M_\odot)$',fontsize=22)
        plt.ylabel('$\log(SFR/M_\odot/yr)$',fontsize=22)
        plt.legend()
        #plot_BV_MS(plt.gca())
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(6.9,11.1)

        # add sSFR cut
        xline = np.linspace(7,10.75,100)
        yline = -12 + xline
        plt.plot(xline,yline,'k--')
    def plot_mstar_hist(self,hsmorphmaxflag=3,rsmorphmaxflag=2):

        plt.figure(figsize=(8,6))        
        # plot add vf galaxies

        self.get_smorph_flag(hsmorphmaxflag=hsmorphmaxflag,rsmorphmaxflag=rsmorphmaxflag)
        haflag = d.smorphflag
        x = self.v.magphys['logMstar_med']
        y = self.v.magphys['logSFR_med']
        flag1 = (x > 7)
        hflag1 = self.v.magphys['logMstar_med'][d.htab['VFINDEX']] > 5
                                                    
        flag1b = flag1 &  ~self.v.halpha['HAflag'] & ~self.v.main['COflag']# valid magphys fit        
        #plt.plot(x[flag1b],y[flag1b],'k.',alpha=.1,label='All VF')        

        # halpha sample
        
        flag3 = flag1 & self.v.main['COflag']
        flag4 = hflag1 & self.smorphflag

        


        flag2 = hflag1 & ~self.smorphflag
        #flag2 = flag1 & self.v.halpha['HAflag'] & ~haflag
        mybins = np.array([ 9.1999979,  9.5999557,  9.9999135, 10.3998713, 10.7998291])
        

        
        #plt.plot(x[flag3],y[flag3],'bo',c=mycolors[1],alpha=.8,label='Primary CO Sample',markersize=5)

        plt.xlabel('$\log(M_\star/M_\odot)$',fontsize=22)
 
        plt.legend()
        #plot_BV_MS(plt.gca())
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        #plt.xlim(6.9,11.1)

        # add sSFR cut
        xline = np.linspace(7,10.75,100)
        yline = -12 + xline
        #plt.plot(xline,yline,'k--')

    def plot_sizeratio_environment(self):
        plt.figure(figsize=(8,6))        
        # plot add vf galaxies
        envs = ['pure_field','poor_group_memb','filament_member','rich_group_memb']#,'cluster_member']
        haflag = self.smorphflag & (self.v.magphys['logsSFR_med'][self.htab['VFINDEX']] > -12)
        sizeratio = self.htab['SMORPH_HR80']/self.htab['SMORPH_R80']
        mybins = np.linspace(0,3,20)
        for i in range(len(envs)):
            eflag = self.v.env[envs[i]][self.htab['VFINDEX']]
            print(envs[i], np.sum(eflag & self.smorphflag))
                                            
            xvar = sizeratio[eflag & haflag]                            
            plt.hist(xvar,alpha=.9,lw=2,bins=mybins,color=mycolors[i],label=envs[i]+f" ({len(xvar)})",histtype='step')
            plt.axvline(np.nanmedian(xvar),color=mycolors[i],ls='--',alpha=.5)
        plt.xlabel("R80(Ha)/R80(r)",fontsize=16)
        plt.ylabel("Number",fontsize=16)        
    
        plt.legend()
        flag0 = self.v.env[envs[0]][self.htab['VFINDEX']]
        flag2 = self.v.env[envs[2]][self.htab['VFINDEX']]        
        self.sizeratio = sizeratio
        res = ks_2samp(sizeratio[flag0 & haflag],sizeratio[flag2 & haflag])
        print(res)

    def plot_sizeratio_conc_environment(self):
        plt.figure(figsize=(8,6))        
        # plot add vf galaxies
        envs = ['pure_field','poor_group_memb','filament_member','rich_group_memb']#,'cluster_member']
        haflag = self.smorphflag & (self.v.magphys['logsSFR_med'][self.htab['VFINDEX']] > -12)
        sizeratio = self.htab['SMORPH_HR80']/self.htab['SMORPH_R80']
        xvariable = self.htab['SMORPH_R20']/self.htab['SMORPH_R80']
        mybins = np.linspace(0,3,20)
        for i in range(len(envs)):
            eflag = self.v.env[envs[i]][self.htab['VFINDEX']]
            print(envs[i], np.sum(eflag & self.smorphflag))
                                            
            yvar = sizeratio[eflag & haflag]
            xvar = self.htab['SMORPH_C'][eflag & haflag]
            #xvar = xvariable[eflag & haflag]            
            plt.subplot(4,1,i+1)
            plt.scatter(xvar,yvar,alpha=.5,c=mycolors[i],label=envs[i])
            plt.legend()
            plt.axis([1,7,0,2])
            #plt.axvline(np.nanmedian(xvar),color=mycolors[i],ls='--',alpha=.5)
            if i == 2:
                plt.ylabel("R80(Ha)/R80(r)",fontsize=16)
            if i < 3:
                plt.xticks([],[])
            print(spearmanr(xvar,yvar))
            xline = np.linspace(1,6,100)
            yline = -0.3*(xline-2)+.9
            plt.plot(xline,yline,'k--',alpha=.4)
        plt.xlabel("SMORPH C",fontsize=16)
        yvar = sizeratio[haflag]
        xvar = self.htab['SMORPH_C'][haflag]
        
        print("all environments: ",spearmanr(xvar,yvar))
    def sbplot_sizeratio_conc_environment(self):
        plt.figure(figsize=(8,6))        
        # plot add vf galaxies
        envs = ['pure_field','poor_group_memb','filament_member','rich_group_memb']#,'cluster_member']
        haflag = self.smorphflag & (self.v.magphys['logsSFR_med'][self.htab['VFINDEX']] > -12)
        sizeratio = self.htab['SMORPH_HR80']/self.htab['SMORPH_R80']
        mybins = np.linspace(0,3,20)
        for i in range(len(envs)):
            eflag = self.v.env[envs[i]][self.htab['VFINDEX']]
            print(envs[i], np.sum(eflag & self.smorphflag))
                                            
            yvar = sizeratio[eflag & haflag]
            xvar = self.htab['SMORPH_C'][eflag & haflag]
            plt.subplot(4,1,i+1)
            plt.scatter(xvar,yvar,alpha=.9,c=mycolors[i],label=envs[i])
            plt.legend()
            plt.axis([1,7,0,2])
            #plt.axvline(np.nanmedian(xvar),color=mycolors[i],ls='--',alpha=.5)
            if i == 2:
                plt.ylabel("R80(Ha)/R80(r)",fontsize=16)
            if i < 3:
                plt.xticks([],[])
            print(spearmanr(xvar,yvar))
        plt.xlabel("SMORPH C",fontsize=16)        
        


        
    def plot_giniratio_environment(self):
        plt.figure(figsize=(8,6))        
        # plot add vf galaxies
        envs = ['pure_field','poor_group_memb','filament_member','rich_group_memb']#,'cluster_member']
        haflag = self.smorphflag & (self.v.magphys['logsSFR_med'][self.htab['VFINDEX']] > -12)
        sizeratio = self.htab['SMORPH_HGINI'] - self.htab['SMORPH_GINI']
        mybins = np.linspace(-.3,.6,20)        
        #sizeratio = self.htab['SMORPH_HC'] - self.htab['SMORPH_C']        
        #mybins = np.linspace(-2,2,20)
        for i in range(len(envs)):
            eflag = self.v.env[envs[i]][self.htab['VFINDEX']]
            print(envs[i], np.sum(eflag & self.smorphflag))
                                            
            xvar = sizeratio[eflag & haflag]                            
            plt.hist(xvar,alpha=.9,lw=2,bins=mybins,color=mycolors[i],label=envs[i]+f" ({len(xvar)})",histtype='step')
            plt.axvline(np.nanmedian(xvar),color=mycolors[i],ls='--',alpha=.5)
        plt.xlabel("Gini(Ha) - Gini(r)",fontsize=16)
        plt.ylabel("Number",fontsize=16)        
    
        plt.legend()
        flag0 = self.v.env[envs[0]][self.htab['VFINDEX']]
        flag2 = self.v.env[envs[2]][self.htab['VFINDEX']]        

        res = ks_2samp(sizeratio[flag0 & haflag],sizeratio[flag2 & haflag])
        print(res)
if __name__ == '__main__':

    ####################################################
    # read in the vf tables
    ####################################################
    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')    
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v2/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_v2_', help = 'prefix for tables; default is vf_v2')                               
    args = parser.parse_args()

    # not sure what this is for...
    if args.tabledir.startswith('/home/rfinn/'):
        homedir = os.getenv("HOME")
        args.tabledir = args.tabledir.replace('/home/rfinn',homedir)

    # read in tables
    v = vtables(args.tabledir,args.tableprefix) 
    v.read_all()


    tabledir = homedir+"/research/Virgo/halpha-tables/"
    #hafilename = os.path.join(tabledir,"hgui_csgrphot_combined.fits")
    hafilename = os.path.join(tabledir,"hgui_csgrphot_combined_2024-Oct-18.fits")    
    googlesheet = tabledir+'hagalaxies-including-duplicates.csv'
    d = hanalysis(hafilename,googlesheet,v)
    d.get_sample()
    d.get_duplicates()
    d.get_smorph_flag()

    
