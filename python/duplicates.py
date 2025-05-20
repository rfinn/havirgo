#!/usr/bin/env python

from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.table import Table, Column
import os
import numpy as np
from scipy.stats import median_abs_deviation as MAD
from scipy import stats
import sys
import collections


homedir = os.getenv("HOME")
sys.path.append(os.path.join(homedir,'github/Virgo/programs/'))
from readtablesv2 import vtables

def pairplot_linear(tab,cols,dupindex1,dupindex2,colorcolumn='M24',\
                        remove_string=None,v1=None,v2=None):
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
        mycolor = np.abs(tab[colorcolumn][dupindex1] - tab[colorcolumn][dupindex2])
        if v1 is not None:
            plt.scatter(x1,x2,c=mycolor,s=20,alpha=.7,vmin=v1,vmax=v2)
        else:
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
    cb.set_label('$\Delta$'+colorcolumn,fontsize=16)        
    #plt.show()


def pairplot_residuals(tab,cols,dupindex1,dupindex2,colorcolumn='M24',remove_string=None,v1=None,v2=None):
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
        mycolor = np.abs(tab[colorcolumn][dupindex1] - tab[colorcolumn][dupindex2])
        if (v1 is not None) and (v2 is not None):
            plt.scatter(x1,dx,c=mycolor,s=20,alpha=.7,vmin=v1,vmax=v2)#,vmin=1,vmax=1.3)
        else:
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
        #mad = np.nanstd(dx)
        if nplot == 5:
            #plt.text(0.05,0.95,f"MED/MAD =\n {med:.2e},{mad:.2e}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
            plt.text(0.05,0.95,f"MED/MAD =\n {med:.2e},{mad:.2e}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
        else:
            #plt.text(0.05,0.95,f"MED/MAD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
            plt.text(0.05,0.95,f"MED/MAD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)            
        #print(f"number of pairs = {npair}")
        plt.axhline(y=mad,color='k',ls=':')
        plt.axhline(y=-mad,color='k',ls=':')        

        nplot += 1
    cb = plt.colorbar(ax=allax, fraction=0.08)
    cb.set_label('$\Delta$'+colorcolumn,fontsize=16)        
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
        self.v = vf
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

    def plot_rstatmorph(self,rsmorphmaxflag=2):
        cols = ['SMORPH_XCENTROID','SMORPH_YCENTROID','SMORPH_RPETRO_CIRC','SMORPH_RPETRO_ELLIP',\
                    'SMORPH_RHALF_ELLIP','SMORPH_R20','SMORPH_R80','SMORPH_GINI',\
                    'SMORPH_M20','SMORPH_F_GM20','SMORPH_S_GM20','SMORPH_C',\
                    'SMORPH_A','SMORPH_S']
        flag =  (self.htab['SMORPH_XCENTROID'] > 0) \
          & (self.htab['SMORPH_S'] > -.5) & (self.htab['SMORPH_A'] > -.5) \
          & (self.htab['SMORPH_FLAG'] < rsmorphmaxflag)

        #flag =  (self.htab['SMORPH_FLAG'] < 1)

            
          
        
        keepflag = flag[self.dupindex1] & flag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=.5,wspace=.35)
        pairplot_residuals(self.htab,cols,dupindex1,dupindex2,colorcolumn='R_FWHM',remove_string='SMORPH_') 
        pairplot_linear(self.htab,cols,dupindex1,dupindex2,colorcolumn='R_FWHM',remove_string='SMORPH_')

    def plot_hstatmorph(self,hsmorphmaxflag=2,remove_tel=None,keep_tel=None):
        cols = ['SMORPH_HXCENTROID','SMORPH_HYCENTROID','SMORPH_HRPETRO_CIRC','SMORPH_HRPETRO_ELLIP',\
                    'SMORPH_HRHALF_ELLIP','SMORPH_HR20','SMORPH_HR50','SMORPH_HR80','SMORPH_HGINI',\
                    'SMORPH_HM20','SMORPH_HF_GM20','SMORPH_HS_GM20','SMORPH_HC',\
                    'SMORPH_HA','SMORPH_HS']
        flag = (self.htab['SMORPH_HXCENTROID'] > 0) & (self.htab['SMORPH_HFLAG'] < hsmorphmaxflag)\
          &  (self.htab['SMORPH_HM20'] > -50) & (self.htab['SMORPH_HS'] > -.5) & (self.htab['SMORPH_HA'] > -1)
        if remove_tel is not None:
            flag = flag & (self.htab['TEL'] != remove_tel)

        if keep_tel is not None:
            flag = flag & (self.htab['TEL'] == keep_tel)
          #& (self.htab['SMORPH_HFLUX_ELLIP'] < 100000)
          
          #&  ##& (self.htab['M24'] > 10)
        
        keepflag = flag[self.dupindex1] & flag[self.dupindex2]
        dupindex1 = self.dupindex1[keepflag]
        dupindex2 = self.dupindex2[keepflag]
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=.5,wspace=.35)

        color = 'H_FWHM'
        color='SMORPH_HFLAG'
        pairplot_linear(self.htab,cols,dupindex1,dupindex2,colorcolumn=color,remove_string='SMORPH_',v1=0,v2=4)
        pairplot_residuals(self.htab,cols,dupindex1,dupindex2,colorcolumn=color,remove_string='SMORPH_',v1=0,v2=4)        
        plt.show()

    def get_smorph_flag(self,hsmorphmaxflag=3,rsmorphmaxflag=2):
        hflag = (self.htab['SMORPH_HXCENTROID'] > 0) \
          & (self.htab['SMORPH_HFLAG'] < hsmorphmaxflag)\
          &  (self.htab['SMORPH_HM20'] > -50) \
          & (self.htab['SMORPH_HS'] > -.5) & (self.htab['SMORPH_HA'] > -1)


        rflag =  (self.htab['SMORPH_XCENTROID'] > 0) \
          & (self.htab['SMORPH_S'] > -.5) & (self.htab['SMORPH_A'] > -.5) \
          & (self.htab['SMORPH_FLAG'] < rsmorphmaxflag)
          
        flag = hflag & rflag
        self.smorphflagh = hflag
        self.smorphflagr = rflag        
        self.smorphflag = flag

    def compare_ratios(self,hsmorphmaxflag=2,rsmorphmaxflag=2):
        """compare the r vs halpha parameters for galaxies that meet both r and halpha smorph flag """
        cols = ['R50','RHALF_ELLIP','RHALF_CIRC',\
                    'R80','RMAX_ELLIP','RMAX_CIRC',\
                    'FLUX_ELLIP','FLUX_CIRC',\
                    'GINI','M20','C','A','S']

        cols = ['R20','R50','R80','RMAX_ELLIP',\
                    'FLUX_ELLIP','FLUX_CIRC',\
                    'GINI','M20','C','A','S']

        self.get_smorph_flag(hsmorphmaxflag=hsmorphmaxflag,rsmorphmaxflag=rsmorphmaxflag)

        flag = self.smorphflag 


        #print(f"\nNumber that meet r and halpha SMORPH flags = {len(dupindex1)}")        
        
        tab = self.htab#[keepflag]

        
        plt.figure(figsize=(10,8))
        plt.subplots_adjust(hspace=.5,wspace=.35)

        nplot = 1
        allax = []
        colorcolumn = 'M24'

        mycolor = self.htab['SMORPH_HFLAG'] + self.htab['SMORPH_FLAG']
        mycolor = self.htab['M24'] 

        npair = np.sum(flag)
        print("number that meet halpha and r smorph flag = ",npair)

        v1 = 10
        v2=18
        remove_string='SMORPH_'
        #v1=None
        for c in cols:
            plt.subplot(3,4,nplot)
            hc = 'SMORPH_H'+c
            rc = 'SMORPH_'+c
            if remove_string is not None:
                ctitle = c.replace(remove_string,'')
            else:
                ctitle = c
                
            if c.startswith('FLUX'):
                x1 = np.log10(tab[rc][flag])
                x2 = np.log10(tab[hc][flag])
                
            else: # take difference
                x1 = tab[rc][flag]
                x2 = tab[hc][flag]
                
            mytitle = f"{ctitle} ({npair})"
            dx = x2 - x1
            #mycolor = 0.5*(tab[colorcolumn][dupindex1] + tab[colorcolumn][dupindex2])
            ccolor = mycolor[flag]
            if v1 is not None:
                plt.scatter(x1,x2,c=ccolor,s=20,alpha=.7,vmin=v1,vmax=v2)
            else:
                plt.scatter(x1,x2,c=ccolor,s=20,alpha=.7)#,vmin=1,vmax=1.3)

            xmin,xmax = plt.xlim()
            xline = np.linspace(xmin,xmax,100)
            plt.plot(xline,xline,'k--')


            plt.title(mytitle)
            allax.append(plt.gca())
            med = np.nanmedian(dx)
            mad =MAD(dx)
            mad = np.nanstd(dx)
            if nplot == 5:
                plt.text(0.05,0.95,f"MED/STD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)
            else:
                plt.text(0.05,0.95,f"MED/STD =\n {med:.2f},{mad:.2f}",transform=plt.gca().transAxes,horizontalalignment="left",verticalalignment='top',fontsize=10)

            plt.plot(xline,xline+mad,'k:')
            plt.plot(xline,xline-mad,'k:')

            #if c.startswith('R'):
            #   plt.gca().set_xscale('log')
            #   plt.gca().set_yscale('log')
            #if c.startswith('FLUX'):
            #   plt.gca().set_xscale('log')
            #   plt.gca().set_yscale('log')                   
            if nplot in [1,5,9]:
                plt.ylabel('Halpha',fontsize=16)
            if nplot > 8:
                plt.xlabel("R-band",fontsize=16)
            nplot += 1
        cb = plt.colorbar(ax=allax, fraction=0.08)
        cb.set_label(colorcolumn,fontsize=16)        

    def remove_duplicates(self):
        dtab = self.htab # has bad galaxies removed
        ftab = dtab#[(dtab['FILT_COR'] < 1.2)]
        fduplist = ([item for item, count in collections.Counter(ftab['VFID']).items() if (count > 1)])

        # create a flag for each subset of measurements
        ephot_flag = (dtab['M26'] > 10)& (dtab['M24'] <20)  
        rad_flag = dtab['R24'] > 0
        galfit_flag = (dtab['GAL_RE'] < 50) & (dtab['GAL_MAG'] < 20)  & (dtab['GAL_N'] < 8) 
        morph_flag = (dtab['ELLIP_ASYM'] > -90)& (dtab['C30'] < 1) 
        smorph_rflag = (dtab['SMORPH_XCENTROID'] > 0) & (dtab['SMORPH_S'] > -.1) & (dtab['SMORPH_FLAG'] < 2)


        ephot_hflag = (dtab['HM16'] > 0)& (dtab['HM17'] >0) 
        rad_hflag = (dtab['HR16'] > 0) & (dtab['HR17'] < 200)
        morph_hflag = (dtab['HC30'] < 1) &(dtab['ELLIP_HASYM'] > -10)&(dtab['ELLIP_HASYM'] < 5) &  (dtab['ELLIP_HGINI'] > 0)  &  (dtab['ELLIP_HGINI']< 1) 
        smorph_hflag = (dtab['SMORPH_HXCENTROID'] > 0)  & (dtab['SMORPH_HFLAG'] < 2)& (dtab['SMORPH_S'] > -.1) &(dtab['SMORPH_HGINI'] > -90) &(dtab['SMORPH_HM20'] > -90) &(dtab['SMORPH_HS'] > -.5)




        galfit_flag = np.array(galfit_flag,'i')
        ephot_flag = np.array(ephot_flag,'i')
        rad_flag = np.array(rad_flag,'i')
        morph_flag = np.array(morph_flag,'i')

        smorph_rflag = np.array(smorph_rflag,'i')
        #goodarea = dtab['ELLIP_UNMASKED_AREA']/dtab['ELLIP_TOTAL_AREA']
        ephot_hflag = np.array(ephot_hflag,'i')
        rad_hflag = np.array(rad_hflag,'i')
        morph_hflag = np.array(morph_hflag,'i')
        smorph_hflag = np.array(smorph_hflag,'i')
        # should compare two values
        filter_flag = (ftab['FILT_COR']< 1.1)        


        # allpointing = ftab['POINTING']
        tocut = []
        neasy = 0
        for vfid in fduplist:
            # find the rows in the dtab
            matchrows = ftab['VFID'] == vfid
            matchindex = np.arange(len(ftab))[matchrows]
            print()
            print(vfid, np.sum(matchrows))
            print()
            rflags = []
            hflags = []
            rsmorph = []
            hsmorph = []
            filtcor = []
            pointing = []
            rfwhm = []
            hfwhm = []
            rskynoise = []
            hskynoise = []
            for i in matchindex:
                rflags.append(np.sum(galfit_flag[i] + ephot_flag[i] + rad_flag[i] + morph_flag[i]))
                hflags.append(np.sum(ephot_hflag[i] + rad_hflag[i] + morph_hflag[i]))
                hsmorph.append(smorph_hflag[i])
                rsmorph.append(smorph_rflag[i])                
                filtcor.append(ftab['FILT_COR'][i])
                pointing.append(ftab['POINTING'][i])
                rfwhm.append(ftab['R_FWHM'][i])
                hfwhm.append(ftab['H_FWHM'][i])
                rskynoise.append(np.log10(ftab['R_SKYNOISE'][i]))
                hskynoise.append(np.log10(ftab['H_SKYNOISE'][i]))
            # pick the winner that has the most r/halpha measurements
            rwin = matchindex[np.where(rflags == np.max(rflags))] 
            hwin = matchindex[np.where(hflags == np.max(hflags))] 
            # pick the winner that has lowest statmorph flags
            rmorphwin = matchindex[np.where(rsmorph == np.min(rsmorph))] 
            hmorphwin = matchindex[np.where(hsmorph == np.min(hsmorph))] 
    
            # pick the winner that has the lowest filter correction
            fwin = matchindex[np.where(filtcor == np.min(filtcor))]
    
            # pick the winner that has the best image quality
            rim = matchindex[np.where(rfwhm == np.min(rfwhm))] 
            him = matchindex[np.where(hfwhm == np.min(hfwhm))] 
    
            # was going to add sky noise as well, but these are not as reliable
            rnoise = matchindex[np.where(rskynoise == np.min(rskynoise))] 
            hnoise = matchindex[np.where(hskynoise == np.min(hskynoise))] 
    
            # this checks to see if one image better filter correction 
            # (then len(fwin = 1), whereas if there is a tie len(fwin)=2)
            # check filter correction - if one is >5% smaller, keep this image
            if (len(fwin) == 1) & (len(filtcor) == 2):
                if np.abs(filtcor[1] - filtcor[0]) > 0.05:
                    # winner is image with smaller filter cor??
                    print("\t\t Winner by filter correction = ",ftab['POINTING'][fwin[0]])
                    for m in matchindex:
                        if m == fwin[0]: # fwin has only one entry
                            continue
                        else:
                            tocut.append(m)
                    neasy += 1
                    continue
    
            # if analysis flags are the same, pick image with best seeing
            if (len(rwin) > 1) & (len(hwin) > 1) & (len(him) == 1):
                print("\t\t Winner by image quality = ",ftab['POINTING'][him[0]])
                neasy += 1
                for m in matchindex:
                    if m == him[0]:
                        continue
                    else:
                        tocut.append(m)
                continue
        
            winners = [rwin,hwin,fwin,rim,him]
            listofwinners = []
            for w in winners:
                for k in w:
                    listofwinners.append(k)
            print("listofwinners = ",listofwinners)
    
            mode = stats.mode(np.array(listofwinners))
            if mode[1][0] > 1:
                print("\t\t Overall winner by mode ",ftab['POINTING'][mode[0][0]])
                neasy += 1
                for m in matchindex:
                    if m == mode[0][0]:
                        continue
                    else:
                        tocut.append(m)
                continue
                                                     
            for w in winners:
                if len(w) < 2:
                    print("\t",ftab['POINTING'][w[0]])
            print("\trflags:", rflags)
            print("\thflags:", hflags)
            print("\tfilt cor:", filtcor)
            print("\tpointing:",pointing)
            #print(f"\tR FWHM: {rfwhm:.2f}")
            #print(f"\tH FWHM: {hfwhm:.2f}")
            #print(f"\tR skynoise: {rskynoise:.2f}")
            #print(f"\tH skynoise: {hskynoise:.2f}")
            print(f"\tR FWHM:",rfwhm)
            print(f"\tH FWHM:" ,hfwhm)
            print(f"\tR skynoise:", rskynoise)
            print(f"\tH skynoise:",hskynoise)
            
            noinput = True
            print("indices = ",matchindex)

        print('rows to cut = ',tocut)
        print("neasy = ",neasy)

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
    #hafilename = os.path.join(tabledir,"hgui_csgrphot_combined_2024-Oct-18.fits")
    hafilename = os.path.join(tabledir,"hgui_csgrphot_combined_2025-May-20.fits")    
    googlesheet = tabledir+'hagalaxies-including-duplicates.csv'
    d = duplicates(hafilename,googlesheet,v)
    d.get_sample()
    d.get_duplicates()
    d.get_smorph_flag()
    
