#!/usr/bin/env

"""
One cutout from INT of VFID5855 has zeros where there is not coverage, as opposed to nans.

Could go back to square one to fix, but I would rather try to just mask these out.


"""

import sys

from astropy.io import fits

imfile = sys.argv[1]
hdu = fits.open(imfile)

# edges have values of zero
mask = hdu[0].data == 0

hdu.close()

# save as a mask
outimage = imfile.replace('.fits','-edgemask.fits')
hdu = fits.PrimaryHDU(1*mask, header=hdu[0].header)
hdu.writeto(outimage, overwrite=True) #sky-subtracted r-band image - use this for photometry
