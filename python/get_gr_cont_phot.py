#!/usr/bin/env python
"""
GOAL:
measure photometry of the continuum subtracted images that are created using the g-r color image


PROCEDURE:
* this will run in each cutout directory
* can be run in parallel


"""
import sys
import os
import date
from astropy.cosmology import WMAP9 as cosmo

homedir = os.getenv("HOME")
sys.path.append(homedir+"/github/halphagui/")

from halphamain import create_output_table

from photwrapper import ellipse
from fit_profile import profile, dualprofile, rprofile, haprofile, ratio_error


class output_table():
    def __init__(self,ngal,ids=None,prefix=None):
        self.ngalaxies = ngal
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
        
    def add_ellipse(self):
        #####################################################################
        # ellipse output
        # xcentroid, ycentroid, eps, theta, gini, sky_centroid, area, background_mean, source_sum, source_sum_err
        #####################################################################
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_EPS',description='axis ratio from ellipse')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_THETA', unit=u.degree,description='position angle from ellipse')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_GINI',description='gini coeff from ellipse')
        e6 = Column(np.zeros(self.ngalaxies), name='ELLIP_GINI2',description='gini coeff method 2')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_M20',description='M20 for r image')
        e6 = Column(np.zeros(self.ngalaxies), name='ELLIP_HM20',description='M20 for Halpha image ')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_AREA',description='area from ellipse')
        e8 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM', unit = u.erg/u.s/u.cm**2,description='total flux from ellipse')
        e9 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM_MAG', unit = u.mag,description='mag from ellipse')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM',description='asym from ellipse')
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM_ERR')
        e12 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM', unit=u.erg/u.s/u.cm**2,description='HA flux from ellipse')
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM_MAG', unit=u.mag,description='HA mag from ellipse')
        e14 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM',description='HA asymmetry from ellipse')
        e15 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM_ERR')
        e16 = Column(np.zeros(self.ngalaxies,'f'), name='R_SKYNOISE',description='R skynoise in erg/s/cm^2/arcsec^2')
        e17 = Column(np.zeros(self.ngalaxies,'f'), name='H_SKYNOISE',description='HA skynoise in erg/s/cm^2/arcsec^2')
        self.table.add_columns([e1,e2,e3,e4,e5,e6, e7,e8, e9, e10, e11, e12, e13,e14,e15,e16,e17])
    def add_photutils(self):
        #####################################################################
        # profile fitting using photutils geometry
        #####################################################################
        #
        # r-band parameters
        #
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
        
        
    def write_fits_table(self):
        if self.prefix is not None:
            # this is not working when running gui - need to feed in the r-band image name
            telescope,dateobs,p = get_params_from_name(self.prefix)
            for i in range(len(self.table)):
                self.table['POINTING'][i] = self.prefix
                self.table['TEL'][i] = telescope
                self.table['DATE-OBS'] = dateobs
        self.table.write(self.output_table, format='fits', overwrite=True)

    def write_ellipse_output(self,e):
        ### SAVE DATA TO TABLE
        fields = ['XCENTROID','YCENTROID','EPS','THETA','GINI','GINI2',\
                  'M20','HM20',
                  'AREA',\
                  'SUM','SUM_MAG','ASYM','ASYM_ERR',\
                  'HSUM','HSUM_MAG','HASYM','HASYM_ERR']#,'SUM_ERR']
        values = [e.xcenter, e.ycenter,e.eps, np.degrees(e.theta), \
                  e.gini,e.gini2,\
                  e.M20_1,e.M20_2,\                  
                  e.cat[e.objectIndex].area.value*self.pixelscale*self.pixelscale,\
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
        wcs = WCS(self.cutout_name_r)
        ra,dec = wcs.wcs_pix2world(e.xcenter,e.ycenter,0)
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
        r30 = e.cat.fluxfrac_radius(0.3)*self.pixelscale*u.arcsec/u.pixel
        r50 = e.cat.fluxfrac_radius(0.5)*self.pixelscale*u.arcsec/u.pixel
        r90 = e.cat.fluxfrac_radius(0.9)*self.pixelscale*u.arcsec/u.pixel

        e.cat.add_extra_property('PHOT_R30',r30)
        e.cat.add_extra_property('PHOT_R50',r50)
        e.cat.add_extra_property('PHOT_R90',r90)

        #c1 = Column(data=np.array(r30[e.objectIndex]),name='PHOTR30',unit='arcsec',description='photutils fluxfrac_radius')
        #c2 = Column(data=np.array(r50[e.objectIndex]),name='PHOTR50',unit='arcsec',description='photutils fluxfrac_radius')
        #c3 = Column(data=r90[e.objectIndex],name='PHOTR90',unit='arcsec',description='photutils fluxfrac_radius')
        #qtable.add_columns([c1,c2,c3])

        qtable = e.cat[e.objectIndex].to_table(colnames)
        
        phot_table_name = self.prefix+'-CS-gr-photuil_tab.fits')
        
        qtable.write(phot_table_name,format='fits',overwrite=True)

    def write_rprofile_fits(self,igal,pfit,prefix=None): # MVC - model
        """ set the prefix='H' for halpha """
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
            self.table[colname][igal]=float('%.2e'%(values[i][0]))
            self.table[colname+'_ERR'][igal]=float('%.2e'%(values[i][1]))
            
    def write_hprofile_fits(self,igal,pfit,gzdist):
        """pass in the profile fit and array of galaxy redshift (flow corrected) """
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

            #print(colname,values[i])
            self.table[colname][igal]=float('{:.2e}'.format(values[i][0]))
            self.table[colname+'_ERR'][igal]=float('{:.2e}'.format(values[i][1]))

        # SFR conversion from Kennicutt and Evans (2012)
        # log (dM/dt/Msun/yr) = log(Lx) - logCx
        logCx = 41.27
        print(len(self.hafit.total_flux),len(gzdist[igal]))
        L = self.hafit.total_flux*(4.*np.pi*cosmo.luminosity_distance(gzdist[igal]).cgs.value**2)
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
        maintab = home+"/research/Virgo/tables-north/v2/vf_v2_main.fits"
        # read in vf_main
        mtab = Table.read(maintab)
        # get redshift
        galindex = mtab['VFID'] == self.vfid
        self.gzdist = mtab['vr'][galindex]/3.e5

    
def fit_profiles(cutout_name_r,cutout_name_ha,prefix=None):
    if prefix is None:
        rphot_table = cutout_name_r.split('.fits')[0]+'-phot.fits'
        haphot_table = cutout_name_ha.split('.fits')[0]+'-phot.fits'
    else:
        rphot_table = cutout_name_r.split('.fits')[0]+'-'+prefix+'-phot.fits'
        haphot_table = cutout_name_ha.split('.fits')[0]+'-'+prefix+'-phot.fits'

    rfit = rprofile(cutout_name_r, rphot_table, label='R')
    rfit.becky_measurements()
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

    # get R and CS-gr image
    rfile = subdirname+'-R.fits'

    hfile = subdirname+'-CS-gr.fits'    

    maskfile = subdirname+'-R-mask.fits'

    
    # setup output table
    prefix = "halpha-csgr-"    
    otab = create_output_table(prefix=prefix,virgo=False,nogui=True)
    # add columns for photutils ephot
    otab.add_ellipse()
    # add columns for my measured radii
    otab.add_photutils()
    
    # run ephot using photwrapper.py
    # the psfs are need by statmorph - ignoring for now...
    e = ellipse(rfile, image2=hfile, mask = maskfile, image_frame = None,image2_filter='ha4', filter_ratio=self.filter_ratio,psf=None,psf_ha=None)
    e.run_for_gui()

    otab.write_ellipse_output(e)

    
    # measure radii
    rfit, hfit = fit_profiles(rfile,hfile)

    vfid = dirname.split('-')[0]
    g = galaxy(vfid)
    g.get_redshift()
    # write output table
    igal = 0 # only doing one galaxy at a time, so set igal to zero...
    otab.write_rprofile_fits(igal,rfit)

    otab.write_hprofile_fits(igal,hfit,np.array(g.zdist))

    otab.write_fits_table()            


    # should I run the photometry on the matched g and r images as well???
    # use our r-band as the reference, measure on legacy g

    # use our r-band as the reference, measure on legacy r
    # let's do it in a separate program
    os.chdir(topdir)


    
    
