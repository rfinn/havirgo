#!/usr/bin/env python


"""
This is my subtract_continuum.py program with Gautam's code to optimize the scale factor used for the r-band image.

It creates an output image -CS-gr-auto.fits

GOAL:
* subtract the continuum using color-correction from legacy images

USEAGE:

python subtract_continuum_auto.py dirname

this assumes the file structure and naming convention used in the Virgo Filament Survey, with 

dirname-R.fits
dirname-Ha.fits
 
mask image:
dirname-R-mask.fits

subdirectory called
legacy/

that has g and r images that have been reprojected onto the halpha image pixel scale


PROCEDURE:
* create a g-r image
   2.5 log10(flux_r/flux_g)

* calculate rms in r-band outside the masked regions

* apply color correction to continuum-subtracted images
  
* determine best scale factor in the range from 0.8 - 1.3

FILTER TRANSFORMATIONS FROM MATTEO

---------------------------------------------------------------------
Number of stars meeting the color and brightness cuts: 54594
---------------------------------------------------------------------
Best fit linear    Ha4 - KPHr = -0.0465 * (PS1_g-PS1_r) + 0.0012
Best fit quadratic Ha4 - KPHr = -0.0156*(PS1_g-PS1_r)^2 + -0.0262*(PS1_g-PS1_r) + -0.0048
---------------------------------------------------------------------
Initial bias/scatter:             -0.0265/0.0158
Corr linear fit              :    0.0000/0.0130
Corr quadratic fit           :    0.0000/0.0130
---------------------------------------------------------------------

---------------------------------------------------------------------
Number of stars meeting the color and brightness cuts: 54692
---------------------------------------------------------------------
Best fit linear    Ha4 - KPSr = -0.1804 * (PS1_g-PS1_r) + 0.0158
Best fit quadratic Ha4 - KPSr = -0.0066*(PS1_g-PS1_r)^2 + -0.1718*(PS1_g-PS1_r) + 0.0133
---------------------------------------------------------------------
Initial bias/scatter:             -0.0916/0.0375
Corr linear fit              :    0.0000/0.0138
Corr quadratic fit           :    0.0000/0.0138
---------------------------------------------------------------------

---------------------------------------------------------------------
Number of stars meeting the color and brightness cuts: 54712
---------------------------------------------------------------------
Best fit linear    Intha - INTSr = -0.2334 * (PS1_g-PS1_r) + 0.0711
Best fit quadratic Intha - INTSr = 0.0160*(PS1_g-PS1_r)^2 + -0.2541*(PS1_g-PS1_r) + 0.0772
---------------------------------------------------------------------
Initial bias/scatter:             -0.0682/0.0491
Corr linear fit              :    0.0000/0.0193
Corr quadratic fit           :    0.0000/0.0192
——————————————————————————————————
"""
import sys
import os
import numpy as np
from astropy.io import fits
from astropy import stats, convolution
from astropy import wcs
from astropy.table import Table

import glob
from reproject import reproject_interp

import warnings
warnings.filterwarnings('ignore')


######################################################################
###  FILTER DEFINITIONS
######################################################################
# Halpha filter width in angstrom
filter_width_AA = {'BOK':80.48,'HDI':80.48,'INT':95,'MOS':80.48,'INT6657':80}

# central wavelength in angstroms
filter_lambda_c_AA = {'BOK':6620.52,'HDI':6620.52,'INT':6568,'MOS':6620.52,'INT6657':6657}


# integral of filter transmission
# calculated in Halpha-paper1.ipynb
filter_Rlambda = {"KPNO_Ha+4nm": 78.58, "WFC_Ha": 84.21, "WFC_Ha6657": 71.96,\
                  "KPNO_R" : 1341.54, "KPNO_r" : 1283.47, "BASS_r": 1042.18, "WFC_r": 1097.07}

# from oversubtraction of continuum, a la Gavazzi+2006
# take telescope as the key
halpha_continuum_oversubtraction = {'BOK':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["BASS_r"]),\
                            'HDI':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["KPNO_r"]),\
                            'INT':(1 +filter_Rlambda["WFC_Ha"]/filter_Rlambda["WFC_r"]),\
                            'MOS':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["KPNO_R"]),\
                            'INT6657':(1 +filter_Rlambda["WFC_Ha6657"]/filter_Rlambda["WFC_r"])}

def getEllipseFocii(xcent, ycent, a, ba, pa):
    pa_rad = pa * np.pi/180.
    b = a*ba
    c = np.sqrt(a**2-b**2)
    x1, y1 = c*np.cos(np.pi/2+pa_rad) + xcent, c*np.sin(np.pi/2+pa_rad) + ycent
    x2, y2 = c*np.cos(3*np.pi/2+pa_rad) + xcent, c*np.sin(3*np.pi/2+pa_rad) + ycent
    return x1, y1, x2, y2

def getEllipseCriterion(x, y, x1, y1, x2, y2, a):
    tot_dist = np.sqrt((x-x1)**2+(y-y1)**2) + np.sqrt((x-x2)**2+(y-y2)**2)
    cond = tot_dist <= 2*a
    return cond

def getEllipseAll(x, y, xcent, ycent, a, ba, pa):
    x1, y1, x2, y2 = getEllipseFocii(xcent, ycent, a, ba, pa)
    cond = getEllipseCriterion(x, y, x1, y1, x2, y2, a)
    return cond

def getPixLength(fits_image):
    hdr = fits_image[0].header
    try:
        cd1_1, cd1_2, cd2_1, cd2_2 = hdr['CD1_1'], hdr['CD1_2'], hdr['CD2_1'], hdr['CD2_2']
        pixlength_x = 3600.0 * np.sqrt(cd1_1**2 + cd2_1**2)
        pixlength_y = 3600.0 * np.sqrt(cd1_2**2 + cd2_2**2)
        return pixlength_x, pixlength_y
    except:
        return 3600.0 * abs(hdr['CDELT1']), 3600.0 * abs(hdr['CDELT2'])


def filter_transformation(telescope,rfilter, gr_col):
    
    """
    use Matteo's linear fits to transform r to Halpha 

    I need to get more info on what these are - are they in mag or flux?

    QFM: is it ok that I am using these transformations on legacy g-r
    when they were derived for panstarrs g-r
    """
    if (telescope == 'BOK') :
        # need to get updated transformation from matteo that is using the
        # BASS r filter, which seems to havesignificantly lower transmission
        # than other r-band filters
        #Ha4_KPSr = -0.1804 * (gr_col) + 0.0158

        # FROM MATTEO 4/15/2024
        #Number of stars meeting the color and brightness cuts: 54842
        #---------------------------------------------------------------------
        #Best fit linear    Ha4 - BASSr = -0.1274 * (PS1_g-PS1_r) + 0.0151
        #Best fit quadratic Ha4 - BASSr = -0.0230*(PS1_g-PS1_r)^2 + -0.0976*(PS1_g-PS1_r) + 0.0063
        #---------------------------------------------------------------------
        #Initial bias/scatter:             -0.0608/0.0276
        #Corr linear fit              :    -0.0000/0.0125
        #Corr quadratic fit           :    0.0000/0.0125
        #---------------------------------------------------------------------
        #ha_r = -0.1804*gr_col + 0.0158
        ha_r = -0.1274*gr_col + 0.0151 # Ha4 vs BASS r
    
    elif ((telescope == 'HDI') and (rfilter == 'r')):
        
        #Ha4_KPSr = -0.1804 * (gr_col) + 0.0158
        ha_r = -0.1804*gr_col + 0.0158
    elif telescope == 'INT': # should be another case for the redder halpha, right?
        #Intha_INTSr = -0.2334 * (gr_col) + 0.0711
        ha_r = -0.2334 * (gr_col) + 0.0711
    elif (telescope == 'MOS') | ((telescope == 'HDI') and (rfilter == 'R')):
        #Ha4_KPHr = -0.0465 * (gr_col) + 0.0012
        ha_r = -0.0465 * (gr_col) + 0.0012
    else:
        print(f"HEY - DID NOT FIND A MATCH FOR TELESOPE {telescope} AND FILTER {rfilter}")

    # replace the nans with zeros
    ha_r_no_nan = np.nan_to_num(ha_r)
    return ha_r_no_nan

def getCorrelation(Halpha, cont):
    return filter2D(Halpha, ddepth=-1, kernel=cont)
     
            
def get_gr(gfile,rfile,mask=None):
    
    """ take g and r filenames, return g-r data and save g-r color image """
    g = fits.open(gfile)
    r = fits.open(rfile)
    data_g = g[0].data
    data_r = r[0].data
    g.close()
    r.close()
    
    ###
    # the following is from Matteo Fossati
    ###
    
    # get noise in the image    
    #stat is a tuple of mean, median, sigma
    #print('\nIn get_gr \nComputing median values for g and r images')
    stat_r = stats.sigma_clipped_stats(data_r,mask=mask)
    print('Subtracting {0:3.2e} from r-band image'.format(stat_r[1]))

    # this is not going to mask out the galaxy, so the sky values will likely be skewed
    # I am going to assume that the legacy images don't need another round of sky subtraction???
    #data_r -= stat_r[1]
    stat_g = stats.sigma_clipped_stats(data_g,mask=mask)
    print('Subtracting {0:3.2e} from g-band image'.format(stat_g[1]))
    #data_g -= stat_g[1]

    # create a mask, where SNR > 10    
    usemask = (data_g>3*stat_g[2])

    # calculate the g-r color 
    gr_col = -2.5*np.log10(data_g/data_r)

    # TODO - should add masking here - we don't want stars to be in our g-r image, right?
    gr_col[mask] = np.nan
    print('Smoothing images for color calculation')
    # changing convolution size from 20 to 10 b/c I'm wondering if it's blurring the color
    # gradients too much - specific example is

    # Testing to remove smoothing on M109
    gr_col = convolution.convolve_fft(gr_col, convolution.Box2DKernel(10), allow_huge=True, nan_treatment='interpolate')

    # set the pixel with SNR < 10 to nan - don't use these for color correction
    gr_col[np.logical_not(usemask)] = np.nan
    
    # save gr color image
    hdu = fits.PrimaryHDU(gr_col, header=r[0].header)
    outimage = rfile.replace('r-ha.fits','gr-ha-smooth.fits')
    #print(f"name for g-r image is {outimage}")
    print(f"writing g-r color image to {outimage}")
    hdu.writeto(outimage, overwrite=True)
    #hdu.close()
    return gr_col

def plot_image(data):
    ###########################################################
    # show the continuum subtracted image
    ###########################################################
    from matplotlib import pyplot as plt
    #from scipy.stats import scoreatpercentile
    from astropy.visualization import simple_norm
    
    plt.figure()
    norm = simple_norm(data, stretch='asinh',max_percent=99,min_percent=.5)
    plt.imshow(data, norm=norm,origin='lower',interpolation='nearest')#,vmin=v1,vmax=v2)
    #plt.show()


# def subtract_continuum(Rfile, Hfile, gfile, rfile, mask=None,overwrite=False,testing=False):
#     """
#     reproject infile to reffile image

#     PARAMS:
#     Rfile : r-band image taken with halpha, to be used for continuum
#     Hfile : halpha image filename
#     gfile : g-band filename to be used for calculating g-r color (legacy image)
#     rfile : r-band filename to be used for calculating g-r color (legacy image)
#     outname: output name 


#     RETURN:
#     nothing, but save the CS subtracted image that uses the g-r color in the current directory
    
#     """
#     outname = Hfile.replace('Ha.fits','CS-gr.fits')
#     if os.path.exists(outname) & (not overwrite):
#         print("continuum-subtracted image exists - not redoing it")
#         return

#     outimage = rfile.replace('r-ha.fits','gr-ha-smooth.fits')
#     print(f"g-r image = {outimage}")
#     if os.path.exists(outimage):
#         print("found g-r image.  not remaking this")
#         hdu = fits.open(outimage)
#         gr_col = hdu[0].data
#         hdu.close()
#     else:
#         gr_col = get_gr(gfile,rfile,mask=mask)

#     # usemask should be all the values in the color image than are
#     # not equal to np.nan
#     usemask = ~np.isnan(gr_col) # these are the good values in the g-r color

    
#     # this should be the text describing the galaxy
#     # like : VFID0569-NGC5989-INT-20190530-p002 
#     fileroot = Rfile.replace('-R.fits','')

#     # read in *our* r-band and halpha images
#     hhdu = fits.open(Hfile)
#     rhdu = fits.open(Rfile)

#     # get photometric ZP for each image
#     rZP = rhdu[0].header['PHOTZP']
#     hZP = hhdu[0].header['PHOTZP']
    
#     # get filter names
#     rfilter = rhdu[0].header['FILTER']
#     hfilter = hhdu[0].header['FILTER']    

#     # TODONE - get the pixel scale in the halpha image
#     # use the WCS function
#     wcs_NB = wcs.WCS(Hfile)
#     pscale_NB = wcs.utils.proj_plane_pixel_scales(wcs_NB)*3600.
    
#     # get the telescope name from the directory/filename
#     tels = ['BOK','INT','HDI','MOS']
#     for t in tels:
#         if t in dirname:
#             telescope = t
#             print('telescope = ',t,dirname)
#             break
#     ##
#     # The following is from Matteo Fossati
#     ##

#     print('Generate NET image')

#     # TODO - revisit this and examine the masking.
#     # the mask we are currently using does not mask the central galaxy

    
#     # TODONE - subtract sky from r-band image
#     print('Computing median values for r and halpha images')
#     stat_r = stats.sigma_clipped_stats(rhdu[0].data,mask=mask)
#     #print('Subtracting {0:3.2f} from r-band image'.format(stat_r[1]))
#     # do I save the r-band image with new sky subtraction???
#     data_r = rhdu[0].data #- stat_r[1]

#     # sky subtracted r-band image
#     skysub_r_name = Rfile.replace('-R.fits','-R-sky.fits')
#     hdu = fits.PrimaryHDU(data_r, header=rhdu[0].header)
#     hdu.writeto(skysub_r_name, overwrite=True) #sky-subtracted r-band image - use this for photometry

    
#     # TODONE - subtract sky from Halpha image
#     # not subtracting sky for now - can try after we get CS to work
#     stat_h = stats.sigma_clipped_stats(hhdu[0].data,mask=mask)
#     #print('Subtracting {0:3.2f} from halpha image'.format(stat_h[1]))
#     data_NB = hhdu[0].data #- stat_h[1]


#     ##
#     # this comments are from matteo's program
#     ##
#     # Generate the r band mag image and the r band calibrated to Halpha wave
#     # This works only for positive flux pixels. Take this into account

#     # TODONE - change ZP - get this from image header
#     # but this will only work for positive flux values - should be ok b/c it's before cont sub

#     # this creates nans where counts are negative
#     mag_r = -2.5*np.log10(data_r) + rZP
#     mag_NB = -2.5*np.log10(data_r) + hZP    

#     # I think we should instead just stay in counts, and scale the counts in both to ZP of 30
#     # now calc fluxes using the same ZP
#     data_r_ZP30 = 10.**(-0.4*(mag_r-30))
#     data_NB_ZP30 = 10.**(-0.4*(mag_NB-30))    

#     # Transform the mag_r image to the observed Halpha filter
#     #
#     # TODONE - change the color conversion -
#     # need to use different conversion for each Halpha/r combo
#     # 
#     # BUT: for pixels which have a nan in the g-r image
#     # (r-band SNR in legacy image < 10)
#     # then mag r to Ha is just mag r
#     #
#     # this doesn't seem right.
#     # seems like it should be the mag r scaled but some default value

    
#     # QFM: these filter transformations are derived for panstarrs colors
#     # is this a problem that we are using the legacy g-r color instead?


#     # The r-band mag has to be scaled to the same ZP as the halpha image


#     # mag_r in AB mags
#     # g-r color is in AB mags
#     # mag_r_to_Ha and mag_r should be very similar - check this
#     mag_r_to_Ha = mag_r + filter_transformation(telescope,rfilter, gr_col)

#     # QFM: what happens to mag_r_to_Ha -
#     # this is the array that has the correct filter transformation, I think

#     ##
#     # the following is the orginal comment from matteo's code
#     ##
#     #Back to calibrated flux units

#     # RF - changing this - will see if I'm doing this correctly
#     #data_r_to_Ha = np.copy(data_r) # this is the sky-subtracted r-band data
#     data_r_to_Ha = data_r_ZP30

#     # smooth the r-band image before subtracting from halpha
#     # QFM : why are we doing this?  I usually do a straight image subtraction
    
#     # QFM: in your code, the following line feeds in data_r rather than data_r_to_Ha
#     # I changed it to data_r_to_Ha.  Is this wrong?
#     # can change smoothing to change to 1-2 psf size
#     data_r_to_Ha = convolution.convolve_fft(data_r_to_Ha, convolution.Box2DKernel(5), allow_huge=True, nan_treatment='interpolate')
    
    
#     # so I don't understand what this line is doing - converting to flux?
#     # QFM: which photometric ZP should I use here?
#     # still the rZP I think - b/c it did a color transformation
#     # but did not scale it, right?

#     # QFM: I think that in the following line, I think the
#     # default values should be adjusted by this amount


#     # this is from matteo's code
#     # data_r_to_Ha[usemask] = 10**(-0.4*(mag_r_to_Ha[usemask]-30))

#     data_r_to_Ha[usemask] = 10**(-0.4*(mag_r_to_Ha[usemask]-30))    

#     # save output to check that color adjustment is working ok
#     hdu = fits.PrimaryHDU(10**(-0.4*(mag_r_to_Ha-rZP)), header=rhdu[0].header)
#     gr_r_name = Rfile.replace('-R.fits','-R-gr-trans.fits')    
#     hdu.writeto(gr_r_name, overwrite=True) #sky-subtracted r-band image - use this for photomet
    

#     # QFM (question for Matteo)
#     # in the above eqn, why are we mixing smoothed and unsmoothed images?
#     # or are we only using the smoothed values for where the pixels are masked?
#     # Rose's answer: I had previously thought that usemask was bad values
#     # but it's actually good.  we are using a color transformation that is smoothed
    
#     # Matteo Comment: Go to cgs units
#     # TODO: check this conversion - need to adjust for my ZP, these are calculated for ZP=30
#     # can choose a ZP to convert to cgs to convert from mag/pix to erg/cm^2/s/Hz
#     # can use ZP=48.6 to be consistent with

#     # apply ZP to NB to get to mag_NB
#     # then convert to data_NB = 10**(-0.4*(mag_r_to_Ha[usemask]-this should be ZP that you choose
#     # - leave it at 30 to be consistent with matteo's numbers below
#     # 10^(-0.4*48.6) to get factor that is now 1E-12
    
#     fnu_NB  = 3.631E3*data_NB_ZP30*1E-12 # now in units of fnu
#     # TODONE - change the filter EWs - need a dictionary for each Halpha filter
#     # factor of 1E18 to avoid large # in image - image is in units of 1E18
#     flam_NB = 2.99792458E-5*fnu_NB/(filter_lambda_c_AA[telescope]**2) *1E18

#     # continuum image - but why are we using the smoothed image?
#     # Rose's answer: same as above.  we smoothed the g-r image so that we apply
#     # a smoother color correction
#     cnu_NB  = 3.631E3*data_r_to_Ha*1E-12
#     clam_NB = 2.99792458E-5*cnu_NB/(filter_lambda_c_AA[telescope]**2) *1E18 # change central wavelength

#     # TODONE - change width of the filter
#     # QFM : in his code, there is a multiplicative factor of 1.03 on clam_NB
#     # need to multiply by width of filter to convert from flux/A to flux
#     # image is  the average flux within the filter
#     # can adjust factor of 1.03 to scale the continuum
#     # in vestige, they check the star subtraction and then adjust the factor to make the stars go away
#     # even with same telescope/filter, this factor can vary
#     flam_net = filter_width_AA[telescope]*(flam_NB-clam_NB) # matteo comment: 106 is the width of the filter


#     # TODO - I would like to save a version in AB mag for compatibility with my photometry programs
#     # QFM - is this just (data_NB - data_r_to_Ha)?
#     NB_ABmag = (data_NB_ZP30 - data_r_to_Ha)

#     # add continuum scale factor to the NB image header

#     hdu = fits.PrimaryHDU(NB_ABmag, header=hhdu[0].header)

#     # outname is *CS-gr.fits
#     hdu.writeto(outname, overwrite=True) #NB image in F_lambda units, before


#     # The rest are different version of the CS image that matteo saves
#     # don't know if I need all of this...
    
#     # DONE: TODO - change output image name
#     # this is the new net NB image
#     hdu = fits.PrimaryHDU(flam_NB, header=hhdu[0].header)
#     hdu.writeto(fileroot+'-net-new.fits', overwrite=True) #NB image in F_lambda units, before
#     # DONE: TODO - change output image name
#     hdu = fits.PrimaryHDU(clam_NB, header=hhdu[0].header)
#     hdu.writeto(fileroot+'-cont-new.fits', overwrite=True)
#     #hdu.close()
    
#     #Calculate clipped statistic
#     #stat is a tuple of mean, median, sigma
#     stat = stats.sigma_clipped_stats(flam_net,mask=mask)

#     # RF - I implemented the sky subtraction of each cutout image in halphagui
#     # would rather keep it there b/c it masks out the central galaxy
#     #flam_net -= stat[1]

#     print('Unbinned SB limit 1sigma {0:3.1f} e-18'.format(stat[2]/(pscale_NB[0]**2)))

#     # DONE: TODO - change output image name
#     # this is the continuum-subtracted image
#     hdu = fits.PrimaryHDU(flam_net, header=hhdu[0].header)
#     hdu.writeto(fileroot+'-net-flux.fits', overwrite=True)
#     #hdu.writeto(outname, overwrite=True)    
#     #hdu.close()

#     # convert image to surface brightness units
#     # DONE: TODO - change pixel scale
#     sblam_net = flam_net/(pscale_NB[0]**2)
    
#     # DONE: TODO - change output image name
#     hdu = fits.PrimaryHDU(sblam_net, header=hhdu[0].header)
#     hdu.writeto(fileroot+'-net-sb.fits', overwrite=True)
#     #hdu.close()

#     print('Smoothing net image')

#     flam_net_smooth = convolution.convolve_fft(flam_net, convolution.Box2DKernel(15), allow_huge=True, nan_treatment='interpolate')

#     hdu = fits.PrimaryHDU(flam_net_smooth, header=hhdu[0].header)
#     # DONE: TODO - change output image name
#     hdu.writeto(fileroot+'-net-smooth.fits', overwrite=True)
#     #hdu.close()
    
#     stat_sm = stats.sigma_clipped_stats(flam_net_smooth,mask=mask)
#     # TODO - add this to image header

#     print('Smoothed {1}x{1} SB limit 1sigma {0:3.1f} e-18'.format(stat_sm[2]/(pscale_NB[0]**2), 15))

#     # close hdu files
#     hhdu.close()
#     rhdu.close()

#     imlist = [mag_r_to_Ha]
#     return imlist
 

if __name__ == '__main__':

    # directory to analyze is specified on the command line
    # this makes the program easy to run with gnu parallel
    dirname = sys.argv[1]

    if len(sys.argv) > 2:
        coloruse = int(sys.argv[2])
    else:
        coloruse = 0

    #print("\nIn subtract continuum, continuum scale factor = ",contscale,"\n")
    # get current directory
    topdir = os.getcwd()

    # move to subdirectory specified in the command line
    os.chdir(dirname)

    # get prefix from the directory name
    t = dirname.split('-')
    prefix = t[0]+'-'+t[1]
    vfid = t[0]

    # get the halpha filter correction from the vf_v2_halpha.fits table
    # need to do this until I update halphamain to add it to the image header
    tabledir = os.getenv("HOME")+'/research/Virgo/tables-north/v2/'
    vhalpha = Table.read(tabledir+'vf_v2_halpha.fits')

    # get ellipse data to use when optimizing continuum subtraction
    ellip_dat = fits.getdata(tabledir+'vf_v2_legacy_ephot.fits')
    sma, pa, ba = ellip_dat['sma_moment'][gindex], ellip_dat['pa_moment'][gindex], ellip_dat['ba_moment'][gindex]
    
    gindex = np.arange(len(vhalpha))[vhalpha['VFID'] == vfid][0]

    # correction for the variation in the halpha filter transmission
    halpha_filter_cor = vhalpha['FILT_COR'][gindex]
    print(f"{vfid}: halpha filter correction = {halpha_filter_cor:.3f}")
    if halpha_filter_cor == 0:
        print("resetting filter correction to 1")
        halpha_filter_cor = 1

    # table with extinction values
    vext = Table.read(tabledir+'vf_v2_extinction.fits')
    # I was multiplying this without first converting from magnitude to flux ratio    
    halpha_extinction_correction = 10.**(vext['A(R)_SandF'][gindex]/2.5)
    
    
    #print("\nhalpha extinction correction = ",halpha_extinction_correction)
    # define the file names
    Rfile = dirname+'-R.fits' # r-band image taken with same telescope as halpha
    Hfile = dirname+'-Ha.fits'  # halpha image

    # get legacy images that are reprojected to the halpha image
    # these are in the legacy subdirectory
    legacy_path = os.path.join('legacy',vfid+'*r-ha.fits')
    rfiles = glob.glob(legacy_path)
    #print(rfiles)
    if len(rfiles) < 1:
        print("problem getting r-ha.fits legacy image",len(rfiles))
        sys.exit()
    else:
        rfile = rfiles[0] # legacy r-band image
        
    # legacy g-band image, shifted to match halpha footprint and pixel scale
    gfiles = glob.glob(os.path.join('legacy',vfid+'*g-ha.fits'))
    if len(gfiles) < 1:
        print("problem getting g-ha.fits legacy image")
        sys.exit()
    else:
        gfile = gfiles[0] # legacy r-band image

    # define the mask file
    maskfile = Rfile.replace('-R.fits','-R-mask.fits')
    if not os.path.exists(maskfile):
        print("WARNING: no mask found")
        mask = None
    else:
        mask = fits.getdata(maskfile)
        mask = mask > 0

    # call the main function to subtract the continuum
    #imlist = subtract_continuum(Rfile, Hfile, gfile, rfile,mask=mask,overwrite=False)


    ## MOVING TO MAIN PROGRAM FOR DEBUGGING
    overwrite=True
    """
    reproject infile to reffile image

    PARAMS:
    Rfile : r-band image taken with halpha, to be used for continuum
    Hfile : halpha image filename
    gfile : g-band filename to be used for calculating g-r color (legacy image)
    rfile : r-band filename to be used for calculating g-r color (legacy image)
    outname: output name 


    RETURN:
    nothing, but save the CS subtracted image that uses the g-r color in the current directory
    
    """
    outname = Hfile.replace('Ha.fits','CS-gr-auto.fits')
    if os.path.exists(outname) & (not overwrite):
        print("continuum-subtracted image exists - not redoing it")
        #return

    outimage = rfile.replace('r-ha.fits','gr-ha-smooth.fits')
    #print(f"g-r image = {outimage}")
    if os.path.exists(outimage) & (not overwrite):
        print("found g-r image.  not remaking this")
        hdu = fits.open(outimage)
        gr_col = hdu[0].data
        hdu.close()
    else:
        gr_col = get_gr(gfile,rfile,mask=mask)

    # usemask should be all the values in the color image than are
    # not equal to np.nan
    usemask = ~np.isnan(gr_col) # these are the good values in the g-r color

    
    # this should be the text describing the galaxy
    # like : VFID0569-NGC5989-INT-20190530-p002 
    fileroot = Rfile.replace('-R.fits','')

    # read in *our* r-band and halpha images
    hhdu = fits.open(Hfile)
    rhdu = fits.open(Rfile)

    # get photometric ZP for each image
    rZP = rhdu[0].header['PHOTZP']
    hZP = hhdu[0].header['PHOTZP']
    
    # get filter names
    rfilter = rhdu[0].header['FILTER']
    hfilter = hhdu[0].header['FILTER']

    # get the r-band scale factor from the r-band header
    rscale = rhdu[0].header['FLTRATIO']

    # TODONE - get the pixel scale in the halpha image
    # use the WCS function
    wcs_NB = wcs.WCS(Hfile)
    pscale_NB = wcs.utils.proj_plane_pixel_scales(wcs_NB)*3600.
    
    # get the telescope name from the directory/filename
    tels = ['BOK','INT','HDI','MOS']
    for t in tels:
        if t in dirname:
            telescope = t
            print('telescope = ',t,dirname)
            break
    ##
    # The following is from Matteo Fossati
    ##

    print('\nGenerate NET image\n')

    # TODO - revisit this and examine the masking.
    # the mask we are currently using does not mask the central galaxy

    
    # TODONE - subtract sky from r-band image
    print('Computing median values for r and halpha images')
    print("subtracting these values from the image...")    
    #print("currently, I am not subtracting these, so check values...")
    stat_r = stats.sigma_clipped_stats(rhdu[0].data,mask=mask)
    print('Subtracting {0:3.2e} from r-band image'.format(stat_r[1]))
    # do I save the r-band image with new sky subtraction???
    data_r = rhdu[0].data - stat_r[1]
    data_r_to_Ha = data_r * rscale

    # sky subtracted r-band image
    skysub_r_name = Rfile.replace('-R.fits','-R-sky.fits')
    hdu = fits.PrimaryHDU(data_r, header=rhdu[0].header)
    hdu.writeto(skysub_r_name, overwrite=True) #sky-subtracted r-band image - use this for photometry

    
    # TODONE - subtract sky from Halpha image
    # not subtracting sky for now - can try after we get CS to work
    stat_h = stats.sigma_clipped_stats(hhdu[0].data,mask=mask)
    print('Subtracting {0:3.2e} from halpha image'.format(stat_h[1]))
    data_NB = hhdu[0].data - stat_h[1]


    ##
    # this comments are from matteo's program
    ##
    # Generate the r band mag image and the r band calibrated to Halpha wave
    # This works only for positive flux pixels. Take this into account

    # TODONE - change ZP - get this from image header
    # but this will only work for positive flux values - should be ok b/c it's before cont sub
    mag_r = -2.5*np.log10(data_r) + rZP
    mag_NB = -2.5*np.log10(data_NB) + hZP    

    # now calc fluxes using the same ZP
    data_r_ZP30 = 10.**(-0.4*(mag_r-30))
    data_NB_ZP30 = 10.**(-0.4*(mag_NB-30))    

    # Transform the mag_r image to the observed Halpha filter
    #
    # TODONE - change the color conversion -
    # need to use different conversion for each Halpha/r combo
    # 
    # BUT: for pixels which have a nan in the g-r image
    # (r-band SNR in legacy image < 10)
    # then mag r to Ha is just mag r
    #
    # this doesn't seem right.
    # seems like it should be the mag r scaled but some default value

    
    # QFM: these filter transformations are derived for panstarrs colors
    # is this a problem that we are using the legacy g-r color instead?


    # The r-band mag has to be scaled to the same ZP as the halpha image

    if coloruse:
        # mag_r in AB mags
        # g-r color is in AB mags
        # mag_r_to_Ha and mag_r should be very similar - check this
        # mag_r_to_Ha = mag_r + filter_transformation(telescope,rfilter, gr_col)

        # going to stay in counts to avoid nans
        # create an image with delta needed to correct for color term
        delta_mag = filter_transformation(telescope,rfilter, gr_col)

        # convert to flux units
        delta_flux = 10.**(-0.4*delta_mag)

        # use the color correction for pixels with sufficient SNR
        data_r_to_Ha[usemask] = data_r_to_Ha[usemask] * delta_flux[usemask]
    
    
    # QFM: what happens to mag_r_to_Ha -
    # this is the array that has the correct filter transformation, I think

    ##
    # the following is the orginal comment from matteo's code
    ##
    #Back to calibrated flux units

    # RF - changing this - will see if I'm doing this correctly
    #data_r_to_Ha = np.copy(data_r) # this is the sky-subtracted r-band data
    #data_r_to_Ha = data_r_ZP30

    # smooth the r-band image before subtracting from halpha
    # QFM : why are we doing this?  I usually do a straight image subtraction
    
    # QFM: in your code, the following line feeds in data_r rather than data_r_to_Ha
    # I changed it to data_r_to_Ha.  Is this wrong?
    # can change smoothing to change to 1-2 psf size

    # skipping this convolution for now
    #data_r_to_Ha = convolution.convolve_fft(data_r_to_Ha, convolution.Box2DKernel(5), allow_huge=True, nan_treatment='interpolate')
    
    
    # so I don't understand what this line is doing - converting to flux?
    # QFM: which photometric ZP should I use here?
    # still the rZP I think - b/c it did a color transformation
    # but did not scale it, right?

    # QFM: I think that in the following line, I think the
    # default values should be adjusted by this amount


    # this is from matteo's code
    # data_r_to_Ha[usemask] = 10**(-0.4*(mag_r_to_Ha[usemask]-30))

    # did this already in ADU/flux units
    #data_r_to_Ha[usemask] = 10**(-0.4*(mag_r_to_Ha[usemask]-30))    

    # save output to check that color adjustment is working ok
    #hdu = fits.PrimaryHDU(data_r_to_Ha/rscale, header=rhdu[0].header)
    #gr_r_name = Rfile.replace('-R.fits','-R-gr-trans.fits')    
    #hdu.writeto(gr_r_name, overwrite=True) #sky-subtracted r-band image - use this for photomet


    # scale both data arrays to ZP = 30?
    # revisit after 

    # QFM (question for Matteo)
    # in the above eqn, why are we mixing smoothed and unsmoothed images?
    # or are we only using the smoothed values for where the pixels are masked?
    # Rose's answer: I had previously thought that usemask was bad values
    # but it's actually good.  we are using a color transformation that is smoothed
    
    # Matteo Comment: Go to cgs units
    # TODO: check this conversion - need to adjust for my ZP, these are calculated for ZP=30
    # can choose a ZP to convert to cgs to convert from mag/pix to erg/cm^2/s/Hz
    # can use ZP=48.6 to be consistent with

    # apply ZP to NB to get to mag_NB
    # then convert to data_NB = 10**(-0.4*(mag_r_to_Ha[usemask]-this should be ZP that you choose
    # - leave it at 30 to be consistent with matteo's numbers below
    # 10^(-0.4*48.6) to get factor that is now 1E-12
    
    fnu_NB  = 3.631E3*data_NB*1E-12 # now in units of fnu
    # TODONE - change the filter EWs - need a dictionary for each Halpha filter
    # factor of 1E18 to avoid large # in image - image is in units of 1E18
    flam_NB = 2.99792458E-5*fnu_NB/(filter_lambda_c_AA[telescope]**2) *1E18

    # continuum image - but why are we using the smoothed image?
    # Rose's answer: same as above.  we smoothed the g-r image so that we apply
    # a smoother color correction
    cnu_NB  = 3.631E3*data_r_to_Ha*1E-12
    clam_NB = 2.99792458E-5*cnu_NB/(filter_lambda_c_AA[telescope]**2) *1E18 # change central wavelength

    #############################################################
    ## GAUTAM'S CODE FOR OPTIMIZING SCALE FACTOR
    #############################################################    

    rnu = cnu_NB / rscale
    rflux = 2.99792458E-5*rnu/(filter_lambda_c_AA_R[telescope]**2) * 1E18 * filter_Rlambda_weff[telescope]

    xpix, ypix = np.arange(wcs_NB._naxis[0]), np.arange(wcs_NB._naxis[1])
    xx, yy = np.meshgrid(xpix, ypix, indexing='xy')
    cond_ellipse = getEllipseAll(xx, yy, xpix[len(xpix)//2], ypix[len(ypix)//2], int(sma/np.average(pscale_NB)), ba, pa)

    contscale_list = np.linspace(0.8, 1.3, 51)
    corvals = np.zeros_like(contscale_list)
    flamenorm = np.zeros_like(contscale_list)
    flamenegs = np.zeros_like(contscale_list)
    overall = np.zeros_like(contscale_list)
    for i, cval in enumerate(contscale_list):
        flam_test = filter_width_AA[telescope]*(halpha_continuum_oversubtraction[telescope]*halpha_filter_cor*flam_NB-cval*clam_NB)
        flame = flam_test[cond_ellipse]
        corvals[i] = np.linalg.norm(getCorrelation(flame, rflux[cond_ellipse]))
        flamenorm[i] = np.linalg.norm(flame, 1)
        flamenegs[i] = np.linalg.norm(flame[flame<0], 1)
        # overall[i] = corvals[i] * flameneg
        
        if i%10==0: 
            # clipped = stats.sigma_clip(flam_test, sigma=50)
            # out = clipped.data[clipped.mask]
            print(f"For contscale={cval:0.2f}, correlation value = {corvals[i]:0.3f}")
            # print(f"For contscale={cval:0.2f}, fraction of negative H-alpha values is {flame[flame<0].size/flame.size}")
            print(f"For contscale={cval:0.2f}, total negative norm over total norm of H-alpha is {flamenegs[i]/flamenorm[i]}")
            # print(f"For contscale={cval:0.2f}, total H-alpha norm over total r-band norm is {np.linalg.norm(flame)/np.linalg.norm(rflux[cond_ellipse])}")
            print(f"For contscale={cval:0.2f}, total norm of H-alpha is {flamenorm[i]} \n")
            # flam_plot = flam_test * 1
            # flam_plot[~cond_ellipse] = 0
            # norm = simple_norm(flam_plot, stretch='asinh',max_percent=99.0,min_percent=1.0)
            # plt.imshow(flam_plot, cmap='viridis', norm=norm)
            # plt.show()
            # breakpoint()
    for i, fl in enumerate(flamenorm):
        overall[i] = 0.5*corvals[i]/corvals.max() + 0.5*flamenorm[i]/flamenorm.max() + flamenegs[i]/flamenegs.max()
    indcor1, indcor2, indcor = np.argmin(corvals), np.argmin(flamenorm), np.argmin(overall)
    cs_use1, cs_use2 = contscale_list[indcor1], contscale_list[indcor2]
    cs_use = contscale_list[indcor]
    print(f"Minimum correlation / total norm / overall for contscale = {cs_use1:0.2f} / {cs_use2:0.2f} / {cs_use:0.2f}")
    # cs_use = (cs_use1 + cs_use2) / 2.0
    
    # TODONE - change width of the filter
    # QFM : in his code, there is a multiplicative factor of 1.03 on clam_NB
    # need to multiply by width of filter to convert from flux/A to flux
    # image is  the average flux within the filter
    # can adjust factor of 1.03 to scale the continuum
    # in vestige, they check the star subtraction and then adjust the factor to make the stars go away
    # even with same telescope/filter, this factor can vary
    flam_net = filter_width_AA[telescope]*(halpha_continuum_oversubtraction[telescope]*halpha_filter_cor*flam_NB-contscale*clam_NB) # matteo comment: 106 is the width of the filter

    # correct for extinction
    flam_net = flam_net * halpha_extinction_correction

    # TODONE - I would like to save a version in AB mag for compatibility with my photometry programs
    NB_ABmag = (halpha_continuum_oversubtraction[telescope]*halpha_filter_cor*data_NB - contscale*data_r_to_Ha)

    # correct for extinction
    NB_ABmag = NB_ABmag * halpha_extinction_correction    


    # this should still be good to use the Halpha ZP
    hhdu[0].header['CONSCALE']=(float(f'{contscale:.3f}'),'Continuum scale factor')    
    hdu = fits.PrimaryHDU(NB_ABmag, header=hhdu[0].header)

    # outname is *CS-gr.fits
    hdu.writeto(outname, overwrite=True) #NB image in F_lambda units, before
    #plot_image(hdu.data)

    ############################################################
    ## COMMENTING OUT WRITING UNTIL I UNDERSTAND WHAT THESE ARE
    ############################################################    
    

    # The rest are different version of the CS image that matteo saves
    # don't know if I need all of this...
    
    # DONE: TODO - change output image name
    # this is the new net NB image

    # add new field to header to track continuum scale factor

    #hdu = fits.PrimaryHDU(flam_NB, header=hhdu[0].header)
    #hdu.writeto(fileroot+'-net-new.fits', overwrite=True) #NB image in F_lambda units, before
    # DONE: TODO - change output image name
    #hdu = fits.PrimaryHDU(clam_NB, header=hhdu[0].header)
    #hdu.writeto(fileroot+'-cont-new.fits', overwrite=True)
    #hdu.close()
    
    #Calculate clipped statistic
    #stat is a tuple of mean, median, sigma
    #stat = stats.sigma_clipped_stats(flam_net,mask=mask)

    # RF - I implemented the sky subtraction of each cutout image in halphagui
    # would rather keep it there b/c it masks out the central galaxy
    #flam_net -= stat[1]

    #print('Unbinned SB limit 1sigma {0:3.2e} e-18'.format(stat[2]/(pscale_NB[0]**2)))

    # DONE: TODO - change output image name
    # this is the continuum-subtracted image
    #hdu = fits.PrimaryHDU(flam_net, header=hhdu[0].header)
    #hdu.writeto(fileroot+'-net-flux.fits', overwrite=True)
    #hdu.writeto(outname, overwrite=True)    
    #hdu.close()

    # convert image to surface brightness units
    # DONE: TODO - change pixel scale
    #sblam_net = flam_net/(pscale_NB[0]**2)
    
    # DONE: TODO - change output image name
    #hdu = fits.PrimaryHDU(sblam_net, header=hhdu[0].header)
    #hdu.writeto(fileroot+'-net-sb.fits', overwrite=True)
    #hdu.close()

    #print('Smoothing net image')
    # QFM: why are we doing this?
    #flam_net_smooth = convolution.convolve_fft(flam_net, convolution.Box2DKernel(15), allow_huge=True, nan_treatment='interpolate')

    #hdu = fits.PrimaryHDU(flam_net_smooth, header=hhdu[0].header)
    # TODONE - change output image name
    #hdu.writeto(fileroot+'-net-smooth.fits', overwrite=True)
    #hdu.close()
    
    #stat_sm = stats.sigma_clipped_stats(flam_net_smooth,mask=mask)
    # TODO - add this to image header

    #print('Smoothed {1}x{1} SB limit 1sigma {0:3.2e} e-18'.format(stat_sm[2]/(pscale_NB[0]**2), 15))

    # close hdu files
    hhdu.close()
    rhdu.close()
    os.chdir(topdir)



    # move back to the top directory
    os.chdir(topdir)

    
