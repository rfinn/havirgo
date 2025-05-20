#!/usr/bin/env python
"""
GOAL:
measure photometry of the continuum subtracted images that are created using the g-r color image


PROCEDURE:
* this will run in each cutout directory
* can be run in parallel

USAGE:
* this takes the directory name as input, like

python ~/github/havirgo/python/get_gr_cont_phot.py VFID6352-NGC5806-BOK-20210418-VFID6406


* to run in parallel

parallel --eta  python ~/github/havirgo/python/get_gr_cont_phot.py :::: virgo-cutouts.txt 


"""
import sys
import os
from datetime import date
import numpy as np

from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table, Column
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

homedir = os.getenv("HOME")
sys.path.append(homedir+"/github/halphagui/")

from halphamain import create_output_table

from photwrapper import ellipse
from fit_profile import profile, dualprofile, rprofile, haprofile, ratio_error

def get_params_from_name(image_name):
    t = os.path.basename(image_name).split('-')
    #print(t)
    if len(t) == 5:
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]
    elif len(t) == 6: # meant to catch negative declinations
        telescope = t[3]
        dateobs = t[4]
        pointing = t[5]
    else:
        print("ruh roh - trouble getting info from ",image_name, len(t))
        print(t)
        for s in t:
            pointing = t[-1]
            if s.startswith('20') & (len(s) == 8):
                dateobs = s
            elif s in ['BOK','INT','MOS','HDI']:
                telescope = s
        print(f"here is what I think: pointing={pointing}, dateobs={dateobs},telescope={telescope}")
                
            
    return telescope,dateobs,pointing


class output_table():
    def __init__(self,ngal,ids=None,prefix=None,pixelscale=None,rwcs=None):
        """
        ngal: not sure what this is or how I was planning to retrofit this to get the values for the new cs images
        from comments above, I was planning to run on each galaxy separately, so ngal should be one
        """
        self.ngalaxies = ngal
        # we are only writing output for one galaxy with this program,
        # so everything goes in row zero
        self.igal = 0
        if ids is None:
            ids = np.ones(1)
        c1 = Column(ids, name='VFID', description='VFID')
        self.table = Table([c1])

        user = os.getenv('USER')
        today = date.today()
        str_date_today = today.strftime('%Y-%b-%d')
        if prefix is None:
            self.output_table = 'data-'+user+'-'+str_date_today+'.fits'
        else:
            self.output_table = prefix+'-'+user+'-'+str_date_today+'.fits'
        self.pixelscale = pixelscale
        self.rwcs = rwcs
        self.prefix = prefix
    def add_ellipse(self):
        #####################################################################
        # ellipse output
        # xcentroid, ycentroid, eps, theta, gini, sky_centroid, area, background_mean, source_sum, source_sum_err
        #####################################################################
        #e0 = Column(np.zeros(self.ngalaxies,'bool'), name='BADGAL',description='bad galaxy flag - maybe partial coverage')
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_EPS',description='axis ratio from ellipse')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_THETA', unit=u.degree,description='position angle from ellipse')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_GINI',description='gini coeff from ellipse')
        e6 = Column(np.zeros(self.ngalaxies), name='ELLIP_HGINI',description='gini coeff method 2')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_M20',description='M20 for r image')
        e8 = Column(np.zeros(self.ngalaxies), name='ELLIP_HM20',description='M20 for Halpha image ')
        e9 = Column(np.zeros(self.ngalaxies,'e'), name='ELLIP_UNMASKED_AREA',description='unmasked source area from photutils')
        e9b = Column(np.zeros(self.ngalaxies,'e'), name='ELLIP_TOTAL_AREA',description='total source area from photutils')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM', unit = u.erg/u.s/u.cm**2,description='total flux from ellipse')
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM_MAG', unit = u.mag,description='mag from ellipse')
        e12 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM',description='asym from ellipse')
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM_ERR')
        e14 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM', unit=u.erg/u.s/u.cm**2,description='HA flux from ellipse')
        e15 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM_MAG', unit=u.mag,description='HA mag from ellipse')
        e16 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM',description='HA asymmetry from ellipse')
        e17 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM_ERR')
        e18 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKYNOISE',description='R skynoise in erg/s/cm^2/arcsec^2')
        e19 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKYNOISE',description='HA skynoise in erg/s/cm^2/arcsec^2')
        e20 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKY',description='R sky level in ADU')
        e21 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKY',description='HA sky level in ADU')

        # photutils radii
        e22 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R30',description='photutils R flux frac 30')
        e23 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R50',description='photutils R flux frac 50')
        e24 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R90',description='photutils R flux frac 90')
        
        e25 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR30',description='photutils Halpha flux frac 30')
        e26 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR50',description='photutils Halpha flux frac 50')
        e27 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR90',description='photutils Halpha flux frac 90')        

        
        self.table.add_columns([e1,e2,e3,e4,e5,e6,e7,e8, e9, e9b,e10, e11, e12, e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26,e27])

    def add_pointing_params(self):
        c13 = Column(np.zeros(self.ngalaxies,dtype='|S60'),name='POINTING', description='string specifying year and pointing')
        c14 = Column(np.zeros(self.ngalaxies,dtype='|S3'),name='TEL', description='telescope/instrument')
        c15 = Column(np.zeros(self.ngalaxies,dtype='i'),name='DATE-OBS', description='string specifying date of observation')                        
        self.table.add_columns([c13,c14,c15])
    def add_photutils(self):
        #####################################################################
        # profile fitting using photutils geometry
        #####################################################################
        #
        # r-band parameters
        #
        self.fields_r = ['R24','R25','R26','R_F25','R24V','R25V','R_F50','R_F75','M24','M25','M26', 'F_30R24','F_R24','C30',\
                    'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        self.units_r = [u.arcsec,u.arcsec,u.arcsec,u.arcsec,u.arcsec,\
                   u.arcsec,u.arcsec,u.arcsec,\
                   u.mag, u.mag, u.mag, \
                   u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,'',\
                   u.arcsec,u.erg/u.s/u.cm**2,u.arcsec, u.arcsec,'',u.mag
                   ]
        self.descriptions= ['isophotal radius at 24mag/sqarc AB',\
                            'isophotal radius at 25mag/sqarc AB',\
                            'isophotal radius at 26mag/sqarc AB',\
                            'radius that encloses 25% of total flux',\
                            'isophotal radius at 24mag/sqarc Vega',\
                            'isophotal radius at 24mag/sqarc Vega',\
                            'radius that encloses 50% of total flux',\
                            'radius that encloses 75% of total flux',\
                            'isophotal mag within R24',\
                            'isophotal mag within R25',\
                            'isophotal mag within R26',\
                            'flux within 30% of R24',\
                            'flux within R24',\
                            'C30 = flux w/in 0.3 r24 / flux w/in r24',\
                            'petrosian radius: where sb is 0.2 times mean sb',\
                            'flux enclosed within 2xpetro radius',\
                            'radius enclosing 50% of petrosian flux',\
                            'radius enclosing 90% of petrosian flux',\
                            '90% petro radius / 50% petro radius',\
                            'magnitude of petrosian flux']
        
        i=0
        for f,unit in zip(self.fields_r,self.units_r):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='ellipse '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=unit,description='ellipse '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        #
        # Halpha parameters
        #

        self.fields_ha = ['R16','R17',\
                  'R_F25','R_F50','R_F75',\
                  'M16','M17', \
                  'F_30R24','F_R24','C30',\
                  'R_F95R24','F_TOT',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG'
                  ]
        self.units_ha = [u.arcsec,u.arcsec,\
                 u.arcsec,u.arcsec, u.arcsec, \
                 u.mag, u.mag, \
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2, '',\
                 u.arcsec,u.erg/u.s/u.cm**2,\
                 u.arcsec,u.erg/u.s/u.cm**2,u.arcsec, u.arcsec,'',u.mag]
        self.descriptions_ha= ['HA isophotal radius at 16erg/s/cm^2',\
                            'HA isophotal radius at 17erg/s/cm^2',\
                            'HA radius that encloses 25% of total flux',\
                            'HA radius that encloses 50% of total flux',\
                            'HA radius that encloses 75% of total flux',\
                            'HA isophotal radius at 16erg/s/cm^s',\
                            'HA isophotal radius at 17erg/s/cm^2',\
                            'HA flux within 30% of R-band R24',\
                            'HA flux within R-band R24',\
                            'HA C30 = flux w/in 0.3 R-band r24 / flux w/in R-band r24',\
                            'HA flux within 30% of R-band R24',\
                            'HA total flux',\
                            'petrosian radius: where sb is 0.2 times mean sb',\
                            'flux enclosed within 2xpetro radius',\
                            'radius enclosing 50% of petrosian flux',\
                            'radius enclosing 90% of petrosian flux',\
                            '90% petro radius / 50% petro radius',\
                            'magnitude of petrosian flux']
        
        i=0
        for f,unit in zip(self.fields_ha,self.units_ha):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f,description='ellipse '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f, unit=unit,description='ellipse '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        f='LOG_SFR_HA'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=u.M_sun/u.yr,description='log10 of HA SFR in Msun/yr')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR',unit=u.M_sun/u.yr)
        c3 = Column(np.zeros(self.ngalaxies,'bool'),name=f+'_FLAG')
        print('testing: colname = ',f+'_FLAG')
        self.table.add_columns([c1,c2,c3])

        ######################################################################
        ### LAST TWO QUANTITIES, I SWEAR!
        ######################################################################        
        
        f='SSFR_IN'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_columns([c1,c2])

        f='SSFR_OUT'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(R) outside 0.3R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_columns([c1,c2])


        self.add_flags()
        
        self.table.add_column(Column(np.zeros(self.ngalaxies,dtype='U50'), name='COMMENT'))
        #print(self.table)

    def add_statmorph(self):
        #####################################################################
        # statmorph output
        #####################################################################

        # rband area
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RPETRO_CIRC', unit='arcsec',description='rpetro circ')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RPETRO_ELLIP', unit='arcsec',description='rpetro ellip')
        e6 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R20', unit='arcsec',description='R20')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R50', unit='arcsec',description='R50')
        e8 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R80', unit='arcsec',description='R80')
        
        e9 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_FLUX_CIRC',description='statmorph flux circ')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RHALF_CIRC',description='statmorph rhalf circ')        
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RMAX_CIRC',description='statmorph rmax circ')
        
        e12 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_FLUX_ELLIP', description='statmorph flux ellip')                
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RHALF_ELLIP',description='statmorph rhalf ellip')        
        e14 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RMAX_ELLIP',description='statmorph rmax ellip')
        
        e15 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_GINI',description='statmorph gini')
        e16 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_M20',description='statmorph M20')        
        e17 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_F_GM20',description='statmorph F(G,M20)')
        e18 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_S_GM20',description='statmorph S(G,M20)')
        e19 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_C',description='statmorph concentration')
        e20 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_A',description='statmorph asymmetry')
        e21 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_S',description='statmorph smoothness')
        e22 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_SERSIC_AMP',description='statmorph sersic amplitude')
        e23 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_RHALF',description='statmorph sersic rhalf')
        e24 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_N',description='statmorph sersic n')
        e25 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_XC',description='statmorph sersic xc')
        e26 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_YC',description='statmorph sersic yc')
        e27 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_ELLIP',description='statmorph sersic ellip')
        e28 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_THETA',description='statmorph sersic theta')
        e29 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SERSIC_CHISQ',description='statmorph sersic chisq dof')
        e30 = Column(np.zeros(self.ngalaxies,'i'), name='SMORPH_SERSIC_FLAG',description='statmorph sersic flag')        
        e31 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_SKY_MEAN',description='statmorph sky mean')
        e32 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_SKY_MED',description='statmorph sky med')
        e33 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_SKY_STD',description='statmorph sky sigma')
        e34 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_SNR_PIXEL',description='statmorph snr per pixes')                       
        e35 = Column(np.zeros(self.ngalaxies,'i'), name='SMORPH_FLAG',description='statmorph flag')                 

        ## Halpha parameters
        h1 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HXCENTROID', unit='pixel',description='xcentroid from ellipse')
        h2 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HYCENTROID', unit='pixel',description='ycentroid from ellipse')
        h3 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRPETRO_CIRC', unit='arcsec',description='rpetro circ')
        h4 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRPETRO_ELLIP', unit='arcsec',description='rpetro ellip')

        h6 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR20', unit='arcsec',description='R20')
        h7 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR50', unit='arcsec',description='R50')
        h8 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR80', unit='arcsec',description='R80')
    
        h9 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_HFLUX_CIRC',description='statmorph flux circ')
        h10 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRHALF_CIRC',description='statmorph rhalf circ')        
        h11 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRMAX_CIRC',description='statmorph rmax circ')
        
        h12 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_HFLUX_ELLIP', description='statmorph flux ellip')                
        h13 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRHALF_ELLIP',description='statmorph rhalf ellip')        
        h14 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRMAX_ELLIP',description='statmorph rmax ellip')
        
        h15 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HGINI',description='statmorph gini')
        h16 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HM20',description='statmorph M20')        
        h17 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HF_GM20',description='statmorph F(G,M20)')
        h18 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HS_GM20',description='statmorph S(G,M20)')
        h19 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HC',description='statmorph concentration')
        h20 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HA',description='statmorph asymmetry')
        h21 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HS',description='statmorph smoothness')
        h22 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_AMP',description='statmorph sersic amplitude')
        h23 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_RHALF',description='statmorph sersic rhalf')
        h24 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_N',description='statmorph sersic n')
        h25 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_XC',description='statmorph sersic xc')
        h26 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_YC',description='statmorph sersic yc')
        h27 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_ELLIP',description='statmorph sersic ellip')
        h28 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_THETA',description='statmorph sersic theta')
        h29 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSERSIC_CHISQ',description='statmorph sersic chisq dof')
        h30 = Column(np.zeros(self.ngalaxies,'i'), name='SMORPH_HSERSIC_FLAG',description='statmorph sersic flag')        
        h31 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_HSKY_MEAN',description='statmorph sky mean')
        h32 = Column(np.zeros(self.ngalaxies,'d'), name='SMORPH_HSKY_MED',description='statmorph sky med')
        h33 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSKY_STD',description='statmorph sky sigma')
        h34 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HSNR_PIXEL',description='statmorph snr per pixes')                      
        h35 = Column(np.zeros(self.ngalaxies,'i'), name='SMORPH_HFLAG',description='statmorph flag')                 
        
        
        
        self.table.add_columns([e1,e2,e3,e4,e6,e7,e8,e9,e10,\
                                e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,\
                                e21,e22,e23,e24,e25,e26,e27,e28,e29,e30,\
                                e31,e32,e33,e34,e35,\
                                h1,h2,h3,h4,h6,h7,h8,h9,h10,\
                                h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,\
                                h21,h22,h23,h24,h25,h26,h27,h28,h29,h30,\
                                h31,h32,h33,h34,h35])


    def add_flags(self):
        '''
        these are common comments that the user will be able to select
        '''
        names = ['CONTSUB_FLAG','MERGER_FLAG','SCATLIGHT_FLAG','ASYMR_FLAG','ASYMHA_FLAG','OVERSTAR_FLAG','OVERGAL_FLAG','PARTIAL_FLAG','EDGEON_FLAG','NUC_HA']
        descriptions =  ['Cont Sub Prob','merger/tidal','scattered light','asymmetric R-band', 'asymmetric Ha','foreground star', 'foreground gal','galaxy is edge-on','galaxy is only partially covered by mosaic','nuclear ha emission']
        for i,n in enumerate(names):
            #print(n)
            c = Column(np.zeros(self.ngalaxies,'bool'),name=n,description=descriptions[i])
            self.table.add_column(c)
        
        
    def write_fits_table(self):
        if self.prefix is not None:
            # this is not working when running gui - need to feed in the r-band image name
            # feed in current working directory name to get params
            gal_id = os.path.basename(os.getcwd())
            telescope,dateobs,p = get_params_from_name(gal_id)
            for i in range(len(self.table)):
                self.table['POINTING'][i] = gal_id
                self.table['TEL'][i] = telescope
                self.table['DATE-OBS'] = dateobs
        self.table.write(self.output_table, format='fits', overwrite=True)

    def write_ellipse_output(self,e):
        ### SAVE DATA TO TABLE
        self.e = e
        fields = ['XCENTROID','YCENTROID','EPS','THETA','GINI','HGINI',\
                  'M20','HM20',\
                  'UNMASKED_AREA','TOTAL_AREA',\
                  'SUM','SUM_MAG','ASYM','ASYM_ERR',\
                  'HSUM','HSUM_MAG','HASYM','HASYM_ERR']#,'SUM_ERR']
        values = [e.xcenter, e.ycenter,e.eps, np.degrees(e.theta), \
                  e.cat.gini[e.objectIndex],e.cat2.gini[e.objectIndex],\
                  e.M20_1,e.M20_2,\
                  e.cat[e.objectIndex].area.value*self.pixelscale*self.pixelscale,\
                  e.masked_pixel_area*self.pixelscale*self.pixelscale,\
                  e.source_sum_erg, e.source_sum_mag,e.asym, e.asym_err, \
                  e.source_sum2_erg,e.source_sum2_mag,e.asym2,e.asym2_err]
        for i,f in enumerate(fields):
            colname = 'ELLIP_'+f
            #print(colname)
            self.table[colname][self.igal]=float('%.2e'%(values[i]))

        # update sky noise
        fields = ['R_SKYNOISE','H_SKYNOISE']
        values = [e.im1_skynoise,e.im2_skynoise]
        for i,f in enumerate(fields):
            print(values[i])
            self.table[colname] = values[i]
        if self.rwcs is not None:
            e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_RA', unit=u.deg,description='R-band center RA from photutil centroid')
            e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_DEC', unit=u.deg,description='R-band center DEC from photutil centroid')
            self.table.add_columns([e1,e2])
            ra,dec = rwcs.wcs_pix2world(e.xcenter,e.ycenter,0)
            self.table['ELLIP_RA'][self.igal]=ra
            self.table['ELLIP_DEC'][self.igal]=dec

        # write out phot table
        colnames = ['area',
                    'background_mean',
                    'bbox_xmax',
                    'bbox_xmin',
                    'bbox_ymax',
                    'bbox_ymin',
                    'cxx',
                    'cxy',
                    'cyy',
                    'eccentricity',
                    'ellipticity',
                    'elongation',
                    'equivalent_radius',
                    'fwhm',
                    'gini',
                    'inertia_tensor',
                    'kron_flux',
                    'kron_fluxerr',
                    'kron_radius',
                    'local_background',
                    'moments', # multi-dimensional
                    'moments_central',
                    'orientation',
                    'perimeter',
                    'segment_flux',
                    'segment_fluxerr',
                    'semimajor_sigma',
                    'semiminor_sigma',
                    'xcentroid',
                    'ycentroid']



        # calculate fractional radii, but these are circular, and in pixels
        print("calculating fluxfrac in write_ellipse_output")
        r30 = e.cat.fluxfrac_radius(0.3)*pixelscale*u.arcsec/u.pixel
        r50 = e.cat.fluxfrac_radius(0.5)*pixelscale*u.arcsec/u.pixel
        r90 = e.cat.fluxfrac_radius(0.9)*pixelscale*u.arcsec/u.pixel

        e.cat.add_extra_property('PHOT_R30',r30)
        e.cat.add_extra_property('PHOT_R50',r50)
        e.cat.add_extra_property('PHOT_R90',r90)

        print("calculating fluxfrac for halpha in write_ellipse_output")
        # this 
        #r30 = e.cat2.fluxfrac_radius(0.3)*pixelscale*u.arcsec/u.pixel
        #r50 = e.cat2.fluxfrac_radius(0.5)*pixelscale*u.arcsec/u.pixel
        #r90 = e.cat2.fluxfrac_radius(0.9)*pixelscale*u.arcsec/u.pixel

        e.cat2.add_extra_property('PHOT_R30',r30)
        e.cat2.add_extra_property('PHOT_R50',r50)
        e.cat2.add_extra_property('PHOT_R90',r90)
        
        print("write fluxfrac to tables")        
        # write these out to the main table
        fields = ['R30','R50','R90']
                  #'HR30','HR50','HR90']
        
        values = [e.cat.PHOT_R30[e.objectIndex].value,\
                  e.cat.PHOT_R50[e.objectIndex].value,\
                  e.cat.PHOT_R90[e.objectIndex].value]
                  #e.cat2.PHOT_R30[e.objectIndex].value,\
                  #e.cat2.PHOT_R50[e.objectIndex].value,\
                  #e.cat2.PHOT_R90[e.objectIndex].value]
        for i,f in enumerate(fields):
            colname = 'ELLIP_'+f
            #print(colname,values[i])
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i].value))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except AttributeError:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except:
                print("problem writing table element",colname,values[i])

        print("writing fits table")
        self.write_fits_table()            
        if e.statmorph_flag:
            print("writing statmorph table")
            self.write_statmorph()
            self.write_fits_table()
        
        #c1 = Column(data=np.array(r30[e.objectIndex]),name='PHOTR30',unit='arcsec',description='photutils fluxfrac_radius')
        #c2 = Column(data=np.array(r50[e.objectIndex]),name='PHOTR50',unit='arcsec',description='photutils fluxfrac_radius')
        #c3 = Column(data=r90[e.objectIndex],name='PHOTR90',unit='arcsec',description='photutils fluxfrac_radius')
        #qtable.add_columns([c1,c2,c3])

        ## RF - getting an error on this part, so skipping for now...

        #qtable = e.cat[e.objectIndex].to_table(colnames)
        #print(qtable)
        #phot_table_name = self.prefix+'-CS-gr-photuil_tab.fits'
        #print("phot_table_name = ",phot_table_name)
        #qtable.write(phot_table_name,format='fits',overwrite=True)

    def write_sky(self):
        e18 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKYNOISE',description='R skynoise in 1E-17 erg/s/cm^2/arcsec^2')
        e19 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKYNOISE',description='HA skynoise in 1E-17  erg/s/cm^2/arcsec^2')
        e20 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKY',description='R sky level in ADU')
        e21 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKY',description='HA sky level in ADU')
        fields = ['R_SKYNOISE','H_SKYNOISE',\
                  'R_SKY','H_SKY']
        
        values = [self.e.im1_skynoise/1.e-17,\
                  self.e.im2_skynoise/1.e-17,\
                  self.e.sky,\
                  self.e.sky2]
        print("\nIn write_sky, skynoise = ",self.e.im1_skynoise/1.e-17,self.e.im2_skynoise/1.e-17)
        for i,colname in enumerate(fields):
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()

    def write_statmorph(self):
        #########################################################
        ## ADD STATMORPH PARAMETERS
        #########################################################

        
        # write these out to the main table
        fields = ['XCENTROID', 'YCENTROID', 'RPETRO_CIRC', 'RPETRO_ELLIP', \
                      'R20','R50','R80',\
                      'FLUX_CIRC','RHALF_CIRC', 'RMAX_CIRC', \
                      'FLUX_ELLIP', 'RHALF_ELLIP', 'RMAX_ELLIP',\
                      'GINI', 'M20','F_GM20','S_GM20','C','A','S',\
                      'SERSIC_AMP','SERSIC_RHALF','SERSIC_N',\
                      'SERSIC_XC','SERSIC_YC','SERSIC_ELLIP','SERSIC_THETA',\
                      #'SERSIC_CHISQ',\
                      #'SERSIC_FLAG',\
                      'SKY_MEAN',\
                      'SKY_MED',\
                      'SKY_STD',\
                      'SNR_PIXEL','FLAG']
                  
        values = [self.e.morph.xc_centroid,\
                  self.e.morph.yc_centroid,\
                  self.e.morph.rpetro_circ*self.pixelscale,\
                  self.e.morph.rpetro_ellip*self.pixelscale,\
                  self.e.morph.r20*self.pixelscale,\
                  self.e.morph.r50*self.pixelscale,\
                  self.e.morph.r80*self.pixelscale,\
                  self.e.morph.flux_circ,\
                  self.e.morph.rmax_circ*self.pixelscale,\
                  self.e.morph.rhalf_circ*self.pixelscale,\
                  self.e.morph.flux_ellip,\
                  self.e.morph.rmax_ellip*self.pixelscale,\
                  self.e.morph.rhalf_ellip*self.pixelscale,\
                  self.e.morph.gini,\
                  self.e.morph.m20,\
                  self.e.morph.gini_m20_bulge,\
                  self.e.morph.gini_m20_merger,\
                  self.e.morph.concentration,\
                  self.e.morph.asymmetry,\
                  self.e.morph.smoothness,\
                  self.e.morph.sersic_amplitude,\
                  self.e.morph.sersic_rhalf*self.pixelscale,\
                  self.e.morph.sersic_n,\
                  self.e.morph.sersic_xc,\
                  self.e.morph.sersic_yc,\
                  self.e.morph.sersic_ellip,\
                  self.e.morph.sersic_theta,\
                  #self.e.morph.sersic_chi2_dof,\
                  #self.e.morph.sersic_flag,\
                  self.e.morph.sky_mean,\
                  self.e.morph.sky_median,\
                  self.e.morph.sky_sigma,\
                  self.e.morph.sn_per_pixel,\
                  self.e.morph.flag]
                  
        for i,f in enumerate(fields):
            colname = 'SMORPH_'+f
            #print(colname)
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()

        ## Add Halpha values
        fields = ['XCENTROID', 'YCENTROID', 'RPETRO_CIRC', 'RPETRO_ELLIP', \
                      'R20','R50','R80',\
                      'FLUX_CIRC','RHALF_CIRC', 'RMAX_CIRC', \
                      'FLUX_ELLIP', 'RHALF_ELLIP', 'RMAX_ELLIP',\
                      'GINI', 'M20','F_GM20','S_GM20','C','A','S',\
                      'SERSIC_AMP','SERSIC_RHALF','SERSIC_N',\
                      'SERSIC_XC','SERSIC_YC','SERSIC_ELLIP','SERSIC_THETA',\
                      #'SERSIC_CHISQ',\
                      #'SERSIC_FLAG',\
                      'SKY_MEAN',\
                      'SKY_MED',\
                      'SKY_STD',\
                      'SNR_PIXEL','FLAG']
                  
        values = [self.e.morph2.xc_centroid,\
                  self.e.morph2.yc_centroid,\
                  self.e.morph2.rpetro_circ*self.pixelscale,\
                  self.e.morph2.rpetro_ellip*self.pixelscale,\
                  self.e.morph2.r20*self.pixelscale,\
                  self.e.morph2.r50*self.pixelscale,\
                  self.e.morph2.r80*self.pixelscale,\
                  self.e.morph2.flux_circ,\
                  self.e.morph2.rmax_circ*self.pixelscale,\
                  self.e.morph2.rhalf_circ*self.pixelscale,\
                  self.e.morph2.flux_ellip,\
                  self.e.morph2.rmax_ellip*self.pixelscale,\
                  self.e.morph2.rhalf_ellip*self.pixelscale,\
                  self.e.morph2.gini,\
                  self.e.morph2.m20,\
                  self.e.morph2.gini_m20_bulge,\
                  self.e.morph2.gini_m20_merger,\
                  self.e.morph2.concentration,\
                  self.e.morph2.asymmetry,\
                  self.e.morph2.smoothness,\
                  self.e.morph2.sersic_amplitude,\
                  self.e.morph2.sersic_rhalf*self.pixelscale,\
                  self.e.morph2.sersic_n,\
                  self.e.morph2.sersic_xc,\
                  self.e.morph2.sersic_yc,\
                  self.e.morph2.sersic_ellip,\
                  self.e.morph2.sersic_theta,\
                  #self.e.morph2.sersic_chi2_dof,\
                  #self.e.morph2.sersic_flag,\
                  self.e.morph2.sky_mean,\
                  self.e.morph2.sky_median,\
                  self.e.morph2.sky_sigma,\
                  self.e.morph2.sn_per_pixel,\
                  self.e.morph2.flag]
                  
        for i,f in enumerate(fields):
            colname = 'SMORPH_H'+f
            #print(colname)
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
        

    def write_rprofile_fits(self,igal,pfit,prefix=None): # MVC - model
        """ set the prefix='H' for halpha """
        self.rfit = pfit
        fields = ['R24','R25','R26','R24V','R25V',\
                  'R_F25','R_F50','R_F75',\
                  'M24','M25','M26',\
                  'F_30R24','F_R24','C30',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        d = pfit
        values = [d.iso_radii[0],d.iso_radii[1],d.iso_radii[2],d.iso_radii[3],d.iso_radii[4],\
                  d.flux_radii[0],d.flux_radii[1],d.flux_radii[2],\
                  d.iso_mag[0],d.iso_mag[1],d.iso_mag[2],\
                  d.flux_30r24,d.flux_r24,d.c30,\
                  d.petrorad,d.petroflux_erg,d.petror50,d.petror90,d.petrocon,d.petromag
                  ]
        for i,f in enumerate(fields):
            if prefix is None:
                colname = f
            else:
                colname = prefix+f
            #print(colname, values[i])
            self.table[colname][igal]=float('%.3e'%(values[i][0]))
            self.table[colname+'_ERR'][igal]=float('%.3e'%(values[i][1]))
            
    def write_hprofile_fits(self,igal,pfit,gzdist,prefix=None):
        """pass in the profile fit and array of galaxy redshift (flow corrected) """

        self.hafit = pfit
        fields = ['R16','R17','R_F25','R_F50','R_F75','M16','M17','F_30R24','F_R24','C30','R_F95R24','F_TOT',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        d = pfit
        values = [d.iso_radii[0],d.iso_radii[1],\
                  d.flux_radii[0],d.flux_radii[1],d.flux_radii[2],\
                  d.iso_mag[0],d.iso_mag[1],\
                  d.flux_30r24,d.flux_r24,d.c30,d.flux_95r24, d.total_flux,\
                  d.petrorad,d.petroflux_erg,d.petror50,d.petror90,d.petrocon,d.petromag
                  ]
        for i,f in enumerate(fields):
            if prefix is None:
                colname = 'H'+f
            else:
                colname = prefix+'H'+f

            # values still have appropriate precision in print statement
            #print(colname,values[i])
            self.table[colname][igal]=float('{:.3e}'.format(values[i][0]))
            self.table[colname+'_ERR'][igal]=float('{:.3e}'.format(values[i][1]))

        # SFR conversion from Kennicutt and Evans (2012)
        # log (dM/dt/Msun/yr) = log(Lx) - logCx
        logCx = 41.27
        #print(len(self.hafit.total_flux),len(gzdist[igal]))
        #print("hafit.total_flux = ",self.hafit.total_flux)
        #print("gzdist = ",gzdist)
        #print("luminosity distance = ",cosmo.luminosity_distance(gzdist[0]:.1f))
        L = self.hafit.total_flux*(4.*np.pi*cosmo.luminosity_distance(gzdist[0]).cgs.value**2)
        #print(L)
        detect_flag = L > 0
        self.sfr = np.zeros(len(L),'d')
        self.sfr[detect_flag] = np.log10(L[detect_flag]) - logCx
        if prefix is None:
            colname='LOG_SFR_HA'
        else:
            colname=prefix+'LOG_SFR_HA'
        #print('sfr = ',self.sfr)
        #print(self.sfr[0], self.sfr[1])
        self.table[colname][igal]=float('%.2e'%(self.sfr[0]))
        self.table[colname+'_ERR'][igal]=float('%.2e'%(self.sfr[1]))
        self.table[colname+'_FLAG'][igal]=detect_flag[0]
        # inner ssfr
        a = self.hafit.flux_30r24
        b = self.rfit.flux_30r24
        self.inner_ssfr = a[0]/b[0]
        self.inner_ssfr_err = ratio_error(a[0],b[0],a[1],b[1])
        if prefix is None:
            colname='SSFR_IN'
        else:
            colname = prefix+'SSFR_IN'
        self.table[colname][igal]=float('%.2e'%(self.inner_ssfr))
        self.table[colname+'_ERR'][igal]=float('%.2e'%(self.inner_ssfr_err))
        # outer ssfr
        c = self.hafit.flux_r24
        d = self.rfit.flux_r24
        self.outer_ssfr = (c[0] - a[0])/(d[0] - b[0])
        self.outer_ssfr_err = ratio_error(c[0] - a[0],d[0] - b[0],np.sqrt(a[1]**2 + c[1]**2),np.sqrt(b[1]**2 + d[1]**2))
        if prefix is None:
            colname='SSFR_OUT'
        else:
            colname=prefix+'SSFR_OUT'
        self.table[colname][igal]=float('%.2e'%(self.outer_ssfr))
        self.table[colname+'_ERR'][igal]=float('%.2e'%(self.outer_ssfr_err))


class galaxy():
    def __init__(self,vfid):
        self.vfid = vfid

    def get_redshift(self):
        # get redshift

        maintab = homedir+"/research/Virgo/tables-north/v2/vf_v2_main.fits"
        # read in vf_main
        mtab = Table.read(maintab)
        # get redshift
        galindex = np.arange(len(mtab))[mtab['VFID'] == self.vfid]
        
        #print(f"in galaxy class, VFID = {self.vfid},galindex = {galindex}")
        #print(f"\t vr = {mtab['vr'][galindex]}")
        self.zdist = mtab['vr'][galindex].value/3.e5

    
def fit_profiles(cutout_name_r,cutout_name_ha,prefix=None):
    if prefix is None:
        rphot_table = cutout_name_r.split('.fits')[0]+'_phot.fits'
        haphot_table = cutout_name_ha.split('.fits')[0]+'_phot.fits'
    else:
        rphot_table = cutout_name_r.split('.fits')[0]+'-'+prefix+'_phot.fits'
        haphot_table = cutout_name_ha.split('.fits')[0]+'-'+prefix+'_phot.fits'

    print("measuring r-band profiles becky")
    rfit = rprofile(cutout_name_r, rphot_table, label='R')
    rfit.becky_measurements()

    print("measuring halpha profiles becky")    
    hafit = haprofile(cutout_name_ha, haphot_table, label=r"$H\alpha$")
    
    hafit.becky_measurements()
    hafit.get_r24_stuff(rfit.iso_radii[rfit.isophotes == 24.][0][0])

    return rfit,hafit



###############################################################
###############################################################
## MAIN PROGRAM
###############################################################
###############################################################
if __name__ == '__main__':

    subdirname = sys.argv[1]
    topdir = os.getcwd()

    # move to subdirectory
    os.chdir(subdirname)


    ###################################################################
    # get galaxy properties from VFID
    ###################################################################    
    vfid = subdirname.split('-')[0]

    maintab = homedir+"/research/Virgo/tables-north/v2/vf_v2_main.fits"
    ephottab = homedir+"/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits"    
    # read in vf_main
    mtab = Table.read(maintab)
    etab = Table.read(ephottab)    
    # get redshift
    galindex = np.arange(len(mtab))[mtab['VFID'] == vfid]
        
    # need RA, DEC, radius, BA, PA, like from halphagui
    ra = mtab['RA'][galindex][0]
    dec = mtab['DEC'][galindex][0]
    print("RA,DEC = ",ra,dec)

    ###################################################################    
    # construct list of apertures from JM's phot file
    ###################################################################    
    JMapertures_arcsec = []
    for i in range(8):
        keyword = f'SMA_AP{i+1:02d}'
        #print(keyword)
        JMapertures_arcsec.append(etab[keyword][galindex][0])

    # this is only 8 apertures, so how to add additional apertures
    # could add one between each aperture
        
    ###################################################################    
    # get R and CS-gr image
    ###################################################################    
    rfile = subdirname+'-R.fits'

    hfile = subdirname+'-CS-gr.fits'    

    maskfile = subdirname+'-R-mask.fits'

    # run ephot using photwrapper.py
    # the psfs are need by statmorph - ignoring for now...
    rheader = fits.getheader(rfile)
    rwcs =  WCS(rheader)
    filter_ratio = rheader['FLTRATIO']
    

    try:
        pixelscale = np.abs(float(rheader['PIXSCAL1'])) # convert deg/pix to arcsec/pixel                        
    except KeyError:
        try:
            pixelscale = np.abs(float(rheader['CD1_1']))*3600. # convert deg/pix to arcsec/pixel
        except KeyError:
            pixelscale = np.abs(float(rheader['PC1_1']))*3600. # Siena pipeline from astronometry.net
    
    
    # setup output table
    prefix = "halpha-csgr"    
    #otab = create_output_table(prefix=prefix,virgo=False,nogui=True)
    ids = [subdirname.split('-')[0]]
    otab = output_table(1, ids=ids,prefix=prefix,pixelscale=pixelscale,rwcs=rwcs)#,virgo=False,nogui=True)    
    # add columns for photutils ephot
    otab.add_ellipse()
    # add columns for my measured radii
    otab.add_photutils()
    otab.add_statmorph()
    otab.add_pointing_params()

    # TODO - get filter for each image
    ##
    # get the halpha filter name
    ##
    if 'INT' in hfile:
        # figure out which INT filter we are using
        hheader = fits.getheader(hfile)
        hfilter = hheader['FILTER']
        if 'Ha6657' in hfilter:
            hfilter = 'intha6657'
        else:
            hfilter = 'inthalpha'
        
    else: # includes BOK, HDI, MOSAIC
        hfilter = '4'

    
    e = ellipse(rfile, image2=hfile, mask = maskfile, image_frame = None,image2_filter='4', filter_ratio=filter_ratio,psf=None,psf_ha=None,objra=ra,objdec=dec,apertures=None)
    
    print("just before run_for_gui\n")
    e.run_for_gui(runStatmorphFlag=True)
    print()
    print("done with run_for_gui")
    otab.write_ellipse_output(e)
    print('wrote output from ellipse')
    
    ###############################################
    # measure radii using custom fitting
    ###############################################
    print("running fit_profiles")
    rfit, hfit = fit_profiles(rfile,hfile)
    print("finished running fit_profiles")    
    dirname = os.path.basename(os.getcwd())
    vfid = dirname.split('-')[0]
    g = galaxy(vfid)
    g.get_redshift()
    # write output table
    igal = 0 # only doing one galaxy at a time, so set igal to zero...

    print("writing profile fits")
    otab.write_rprofile_fits(igal,rfit)

    print("writing hprofile fits")    
    #print(f"g.zdist = {g.zdist}")
    otab.write_hprofile_fits(igal,hfit,np.array(g.zdist))
    otab.write_sky()
    otab.write_fits_table()            

    print("plotting fancy_profiles")
    otab.e.plot_fancy_profiles()
    # should I run the photometry on the matched g and r images as well???
    # use our r-band as the reference, measure on legacy g

    # use our r-band as the reference, measure on legacy r
    # let's do it in a separate program
    os.chdir(topdir)


    
