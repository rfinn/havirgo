#!/usr/bin/env python

from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.table import Table, Column
import os
import numpy as np
from scipy.stats import median_abs_deviation as MAD
import sys
import collections


homedir = os.getenv("HOME")
sys.path.append(os.path.join(homedir,'github/Virgo/programs/'))
from readtablesv2 import vtables

def pairplot_linear(tab,cols,dupindex1,dupindex2,colorcolumn='M24',remove_string=None):
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(hspace=.5,wspace=.35)

    nplot = 1
    allax = []

    npair = len(dupindex1)
    mycolor = tab[colorcolumn]
    for c in cols:
        plt.subplot(4,4,nplot)
        x1 = tab[c][dupindex1]
        x2 = tab[c][dupindex2]
        dx = x2 - x1
        colorcolumn = colorcolumn
        mycolor = 0.5*(tab[colorcolumn][dupindex1] +tab[colorcolumn][dupindex2])
        plt.scatter(x1,x2,c=mycolor,s=20,alpha=.7)#,vmin=1,vmax=1.3)

        xmin,xmax = plt.xlim()
        xline = np.linspace(xmin,xmax,100)
        plt.plot(xline,xline,'k--')
        #plt.title(c)
        if remove_string is not None:
            ctitle = c.replace(remove_string,'')
        else:
            ctitle = c
        plt.title(f"{ctitle} ({npair})")
        allax.append(plt.gca())
        med = np.nanmedian(dx)
        mad =MAD(dx)
        #mad = np.nanstd(dx)
        if nplot == 5:
            plt.text(0.05,0.95,f"MED/MAD =\n {med:.2e},{mad:.2e}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
        else:
            plt.text(0.05,0.95,f"MED/MAD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
        #print(f"number of pairs = {npair}")
        plt.plot(xline,xline+mad,'k:')
        plt.plot(xline,xline-mad,'k:')
        nplot += 1
    cb = plt.colorbar(ax=allax, fraction=0.08)
    cb.set_label(colorcolumn,fontsize=16)        
    #plt.show()


def pairplot_residuals(tab,cols,dupindex1,dupindex2,colorcolumn='M24',remove_string=None):
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(hspace=.5,wspace=.35)

    nplot = 1
    allax = []

    npair = len(dupindex1)        
    for c in cols:
        plt.subplot(4,4,nplot)
        x1 = tab[c][dupindex1]
        x2 = tab[c][dupindex2]
        dx = x2 - x1
        colorcolumn = colorcolumn
        mycolor = tab[colorcolumn][dupindex1]
        plt.scatter(x1,dx,c=mycolor,s=20,alpha=.7)#,vmin=1,vmax=1.3)

        plt.axhline(y=0,color='k',ls='--')        
        #plt.title(c)
        if remove_string is not None:
            ctitle = c.replace(remove_string,'')
        else:
            ctitle = c
        plt.title(f"{ctitle} ({npair})")
        
        allax.append(plt.gca())
        med = np.nanmedian(dx)
        mad =MAD(dx)
        #mad = np.std(dx)
        mad = np.nanstd(dx)
        if nplot == 5:
            #plt.text(0.05,0.95,f"MED/MAD =\n {med:.2e},{mad:.2e}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
            plt.text(0.05,0.95,f"MED/STD =\n {med:.2e},{mad:.2e}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
        else:
            #plt.text(0.05,0.95,f"MED/MAD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
            plt.text(0.05,0.95,f"MED/STD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)            
        #print(f"number of pairs = {npair}")
        plt.axhline(y=mad,color='k',ls=':')
        plt.axhline(y=-mad,color='k',ls=':')        

        nplot += 1
    cb = plt.colorbar(ax=allax, fraction=0.08)
    cb.set_label(colorcolumn,fontsize=16)        
    #plt.show()
    
class duplicates():
    def __init__(self,hatab_filename,googlesheet,vf):
        """
        ARGS:
        * hatab = stacked output from running photwrapper on CS gr images
        * googlesheet = spreadsheet with ratios and list of galaxies to exclude

        * vf = main vf tables
        """
        self.fulltab = Table.read(hatab_filename)

        fullvfid = np.zeros(len(self.fulltab),'i')
        for i in range(len(self.fulltab)):
            fullvfid[i] = int(self.fulltab['VFID'][i].replace('VFID',''))

        # this is the corresonding row in the full halpha tables
        self.vfindex = fullvfid

        # read in spreadsheet with ratios and list of galaxies to exclude
        self.gspread = Table.read(googlesheet)

        tel_int = 1*(self.fulltab['TEL'] == 'BOK') + \
          2*(self.fulltab['TEL'] == 'HDI') + \
          3*(self.fulltab['TEL'] == 'INT') + \
          4*(self.fulltab['TEL'] == 'MOS') 
        self.fulltab.add_column(Column(tel_int,name='TELNUM'))
    def get_sample(self):
        print("number of galaxies observed in halpha = ",len(self.fulltab))
        self.htab = self.fulltab[~self.fulltab['badflag']]
        self.hvfindex = self.vfindex[~self.fulltab['badflag']]
        print("number after removing bad flag = ",len(self.htab))

        self.rflag =  (self.htab['M26'] > 10)& (self.htab['M24'] <20) \
            & (self.htab['C30'] < 1) \
            & (self.htab['GAL_RE'] < 50) & (self.htab['GAL_MAG'] < 20)  \
            & (self.htab['GAL_N'] < 8) & (self.htab['GAL_MAG'] > 0) \
            & (self.htab['ELLIP_ASYM'] > -90) & (self.htab['PETRO_MAG'] > 0)

        

        self.hflag = (self.htab['HM16'] > 0)& (self.htab['HM17'] >0) \
            & (self.htab['HC30'] < 1) &(self.htab['ELLIP_HASYM'] > -10)\
            & (self.htab['ELLIP_HASYM'] < 5) \
            & (self.htab['ELLIP_HGINI'] > 0)  &  (self.htab['ELLIP_HGINI']< 1) \
            & (self.htab['HR17'] < 200)

        
    def get_duplicates(self):
        '''get indices of duplicates in the sample '''
        self.dupindex1 = []
        self.dupindex2 = []

        # this will contain the vfid for any galaxy observed more than once
        duplist2 = ([item for item, count in collections.Counter(self.htab['VFID']).items() if (count > 1)])

        # collect indices in htab for all the duplicates
        for vfid in duplist2:
            matchflag = self.htab['VFID'] == vfid
            indices = np.arange(len(self.htab))[matchflag]
            
            if len(indices) > 1:
                self.dupindex1.append(indices[0])
                self.dupindex2.append(indices[1])
                
            if len(indices) > 2:
                self.dupindex1.append(indices[0])
                self.dupindex2.append(indices[2])
                
                self.dupindex1.append(indices[1])
                self.dupindex2.append(indices[2])                

        # convert to integer arrays
        self.dupindex1 = np.array(self.dupindex1,'i')
        self.dupindex2 = np.array(self.dupindex2,'i')        
        print(f"number of duplicate observations = {len(self.dupindex1)}")

        
    def plot_rparams(self):
        cols = ['ELLIP_GINI','ELLIP_M20','C30','ELLIP_ASYM',\
                'ELLIP_SUM','PETRO_R','R24','R25',\
                'GAL_MAG','GAL_RE','GAL_N','GAL_BA',\
                'M24','M25','M26','PETRO_MAG']
        keepflag = self.rflag[self.dupindex1] & self.rflag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=.5,wspace=.35)

        pairplot_linear(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24')
        plt.show()


    def plot_rparams_residuals(self):
        cols = ['ELLIP_GINI','ELLIP_M20','C30','ELLIP_ASYM',\
                'ELLIP_SUM','PETRO_R','R24','R25',\
                'GAL_MAG','GAL_RE','GAL_N','GAL_BA',\
                'M24','M25','M26','PETRO_MAG']


        keepflag = self.rflag[self.dupindex1] & self.rflag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]

        pairplot_residuals(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24')
        plt.show()

    def plot_hparams(self):
        cols = ['ELLIP_HGINI','ELLIP_HASYM','ELLIP_HM20','HC30',\
        'ELLIP_HSUM','HPETRO_R','HR16','HR17',\
        'HM16','HM17','HF_TOT','HPETRO_MAG']

        keepflag = self.hflag[self.dupindex1] & self.hflag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=.5,wspace=.35)

        pairplot_linear(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24')
        plt.show()


    def plot_hparams_residuals(self):
        cols = ['ELLIP_HGINI','ELLIP_HASYM','ELLIP_HM20','HC30',\
        'ELLIP_HSUM','HPETRO_R','HR16','HR17',\
        'HM16','HM17','HF_TOT','HPETRO_MAG']


        keepflag = self.hflag[self.dupindex1] & self.hflag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]

        pairplot_residuals(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24')
        plt.show()

    def plot_rstatmorph(self):
        cols = ['SMORPH_XCENTROID','SMORPH_YCENTROID','SMORPH_RPETRO_CIRC','SMORPH_RPETRO_ELLIP',\
                    'SMORPH_RHALF_ELLIP','SMORPH_R20','SMORPH_R80','SMORPH_GINI',\
                    'SMORPH_M20','SMORPH_F_GM20','SMORPH_S_GM20','SMORPH_C',\
                    'SMORPH_A','SMORPH_S']
        flag =  (self.htab['SMORPH_XCENTROID'] > 0) \
          & (self.htab['SMORPH_S'] > -.5) & (self.htab['SMORPH_A'] > -.5) \
          & (self.htab['SMORPH_FLAG'] < 2)

        #flag =  (self.htab['SMORPH_FLAG'] < 1)

            
          
        
        keepflag = flag[self.dupindex1] & flag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=.5,wspace=.35)

        pairplot_linear(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24',remove_string='SMORPH_')
        pairplot_residuals(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24',remove_string='SMORPH_H')                
    def plot_hstatmorph(self):
        cols = ['SMORPH_HXCENTROID','SMORPH_HYCENTROID','SMORPH_HRPETRO_CIRC','SMORPH_HRPETRO_ELLIP',\
                    'SMORPH_HRHALF_ELLIP','SMORPH_HR20','SMORPH_HR50','SMORPH_HR80','SMORPH_HGINI',\
                    'SMORPH_HM20','SMORPH_HF_GM20','SMORPH_HS_GM20','SMORPH_HC',\
                    'SMORPH_HA','SMORPH_HS']
        flag = (self.htab['SMORPH_HXCENTROID'] > 0) & (self.htab['SMORPH_HFLAG'] < 3)\
          &  (self.htab['SMORPH_HM20'] > -50) & (self.htab['SMORPH_HS'] > -.5) & (self.htab['SMORPH_HA'] > -1)
          #& (self.htab['SMORPH_HFLUX_ELLIP'] < 100000)
          
          #&  ##& (self.htab['M24'] > 10)
        
        keepflag = flag[self.dupindex1] & flag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=.5,wspace=.35)

        pairplot_linear(self.htab,cols,dupindex1,dupindex2,colorcolumn='M24',remove_string='SMORPH_H')
        pairplot_residuals(self.htab,cols,dupindex1,dupindex2,colorcolumn='TELNUM',remove_string='SMORPH_H')        
        plt.show()        
        
        
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
    hafilename = os.path.join(tabledir,"hgui_csgrphot_combined_2024-Oct-17.fits")    
    googlesheet = tabledir+'hagalaxies-including-duplicates.csv'
    d = duplicates(hafilename,googlesheet,v)
    d.get_sample()
    d.get_duplicates()
    
