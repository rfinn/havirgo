#!/usr/bin/env python

import astropy.units as u
from radio_beam import Beam
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
import numpy as np


def get_convolution_kernel(hafwhm, beamsize, pixelscale):
    """
    Create and return a convolution kernel that aligns with radio beam but has size needed to transform the halpha image to radio

    INPUT
    - hafwhm = halpha image FWHM, in units of arcsec 
    - beamsize = list containing beam size parameters as [a_arcsec, b_arcsec, PA_deg]
    - pixelscale = pixel scale of halpha image in units of arcsec
    RETURN:
    - ellip_kernel
    """

    ka_arcsec = np.sqrt(beamsize[0]**2 - hafwhm**2)
    kb_arcsec = np.sqrt(beamsize[1]**2 - hafwhm**2)

    my_beam = Beam(ka_arcsec*u.arcsec, kb_arcsec*u.arcsec, beamsize[2]*u.deg)

    pix_scale = pixelscale * u.arcsec
    ellip_kernel = my_beam.as_kernel(pix_scale)
    return ellip_kernel


def convolve_image(image, imagefwhm, radiobeam_a, radiobeam_b, radiobeam_PA):
    """  
    INPUT:
    image
    imagefwhm
    radiobeam_a
    radiobeam_b
    radiobeam_PA

    RETURN
    * output image name
    * convolved image data
    """
    # read in input image
    imhdu = fits.open(image)
                       
    # get pixel scale of input image
    imwcs = WCS(imhdu[0].header)
    pscale = wcs.utils.proj_plane_pixel_scales(imwcs) # in deg per pixel
    pixel_scale = pscale[0]*3600
    
    # get FWHM of input image
    im_fwhm = float(imagefwhm)
    
    # get beam
    radio_beam = [float(radiobeam_a),float(radiobeam_b),float(radiobeam_PA)]

    # get convolution kernel
    mykernel = get_convolution_kernel(im_fwhm, radio_beam, pixel_scale)
    
    # convolve halpha with elliptical convolution kernel
    convolved_data = convolve(imhdu[0].data, mykernel)

    # write out convolved image    
    chdu = fits.PrimaryHDU(convolved_data, header=imhdu[0].header)
    chdu.writeto('c'+image, overwrite=True)

    return 'c'+image, convolved_data 

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ="Convolve input optical image to match the elliptical psf of a radio image.")
    parser.add_argument('--image', dest = 'image', default =None, help = 'image filename to convolve.')
    parser.add_argument('--imagefwhm', dest = 'imagefwhm', default = None, help = 'image fwhm in arcsec')
    
    parser.add_argument('--radiobeam_a', dest = 'radiobeam_a', default = None,  help = 'radio beam semi-major axis a_arcsec')
    parser.add_argument('--radiobeam_b', dest = 'radiobeam_b', default = None,  help = 'radio beam semi-minor axis b_arcsec')
    parser.add_argument('--radiobeam_PA', dest = 'radiobeam_PA', default = None,  help = 'radio beam PA in deg')    
    
    
    args = parser.parse_args()    

    radiobeam_a, radiobeam_b, radiobeam_PA = float(args.radiobeam_a),float(args.radiobeam_b),float(args.radiobeam_PA)

    
    t = convolve_image(args.image, float(args.imagefwhm), radiobeam_a, radiobeam_b, radiobeam_PA)

    
