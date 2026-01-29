#!/usr/bin/env python

"""
GOAL:
update headers of Becky's images

add:
* RA
* DEC
* rough WCS

Tables from papers are in ~/research/Virgo/koopman-images/paper-tables/

Cluster data:
KKY01 table 1 - RA and DEC (1950)
KKY01 table 3 - gives detector
KKY01 table 4 - gives pixel scale for each detector


isolated data:
KK06 table 2 - name, RA and DEC (1950)
KK06 table 3 - name, gives detector
KK06 table 4 - gives pixel scale for each detector


method:
* get RA and DEC
* get pixelscale
* get filter
* update header

"""

import argparse
import os
from astropy.io import fits
import sys

homedir = os.getenv("HOME")
tabledir = os.path.join(homedir,'research/Virgo/koopmann-images/paper-tables/')

def get_coords(galname):
    """
    INPUT:
    * galname: e.g. NGC4178, IC3392

    PROCEDURE:
    read in KKY01 table 1
    find line starting with gal name
    get coord strings
    convert to (RA,DEC) in J2000 deg

    RETURN:
    ra: in deg (J2000)
    dec: in deg (J2000)
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.coordinates import FK5
    from astropy.time import Time

    foundMatch = False
    
    # open file
    input = open(os.path.join(tabledir,'KKY01-table1.txt'),'r')
    # search for line starting with gal name
    # when gal is found, parse strings that hold ra and dec
    for line in input:
        if line.startswith('NGC') | line.startswith('IC') | line.startswith('UGC'):
            t = line.split()
            testname = t[0]+t[1]
            if testname == galname: # t[1] is the NGC number
                ra = ':'.join(t[2:5])
                dec = ':'.join(t[5:8])
                c = SkyCoord(ra,dec,unit=(u.hourangle,u.deg),frame=FK5(equinox=Time('J1950')))
                foundMatch = True
                break
    input.close()

    # if no match found in cluster table, look in isolated galaxies table
    if not foundMatch: 
        # open file
        input = open(os.path.join(tabledir,'KK06-table2.txt'),'r')
        # search for line starting with gal name
        # when gal is found, parse strings that hold ra and dec
        for line in input:
            if line.startswith('NGC') | line.startswith('IC') | line.startswith('UGC'):
                t = line.split()
                
                testname = t[0]+t[1]
                if testname == galname: # t[1] is the NGC number
                    ra = ':'.join(t[2:5])
                    dec = ':'.join(t[5:8])
                    print(f"ra = {ra}, dec = {dec}")
                    c = SkyCoord(ra,dec,unit=(u.hourangle,u.deg),frame=FK5(equinox=Time('J1950')))
                    foundMatch = True
                    break
        input.close()

    if foundMatch:
        # transform to J2000
        c2000 = c.transform_to(FK5(equinox='J2000'))
        ra, dec = c2000.ra.value, c2000.dec.value
    else:
        ra, dec = None, None

    return ra, dec


def get_instrument(galname):
    """
    INPUT:
    * galname: e.g. NGC4178, IC3392

    PROCEDURE:

    RETURN:
    * instrument
    """
    foundMatch = False

    if galname in ['NGC4298', 'NGC4302', 'NGC4567', 'NGC4568', 'NGC4647', 'NGC4649']:
        return 'TEK1'
    elif galname in ['NGC4383', 'UGC7504', 'NGC4606', 'NGC4607', 'NGC4694']:
        return 't2ka'
    # open file
    input = open(os.path.join(tabledir,'KKY01-table3.txt'),'r')
    # search for line starting with gal name
    # when gal is found, parse strings that hold ra and dec
    for line in input:
        if line.startswith('NGC') | line.startswith('IC') | line.startswith('UGC'):
            t = line.split()
            testname = t[0]+t[1]
            if '/' in testname:
                testname = t[0]+t[1].split('/')[0]
            if testname == galname: # t[1] is the NGC number
                instrument = t[5].split('/')[1]
                foundMatch = True
                break
    input.close()

    # if no match found in cluster table, look in isolated galaxies table
    if not foundMatch: 
        # open file
        input = open(os.path.join(tabledir,'KK06-table3.txt'),'r')
        # search for line starting with gal name
        # when gal is found, parse strings that hold ra and dec
        for line in input:
            if line.startswith('NGC') | line.startswith('IC') | line.startswith('UGC'):
                t = line.split()
                testname = t[0]+t[1]
                if testname == galname: # t[1] is the NGC number
                    instrument = t[5].split('/')[1]
                    foundMatch = True
                    break
        input.close()

    if not foundMatch:
        instrument = None
    if instrument == 'S2KB':
        instrument = 's2kb'
    elif instrument == 't2kA':
        instrument = 't2ka'

    return instrument

def get_pixel_scale(instrument):
    """
    get pixel scale from file

    KKY01 table 4 - gives pixel scale for each detector
    KK06 table 4 - gives pixel scale for each detector

    RETURN
    pixelScale : pixel scale in arcsec per pixel

    """

    # making a dictionary in case that's easier
    pixel_dictionary = {
        'TI2': 0.86,
        'TEK1': 0.77,
        't2ka': 0.68,
        'TEK2K': 0.40,
        'TEK1K': 0.40,
        'TEK1K-1': 0.40,
        's2kb': 0.20}
    
    foundMatch = False
    try:
        pixelScale = pixel_dictionary[instrument]
        foundMatch = True
    except KeyError:
        print(f"WARNING: did not find pixel scale for instrument:{instrument}")
        sys.exit()
    
    return pixelScale

def get_filters(galname):
    """ match galaxy to its corresponding r and halpha filter """
    # TODO: read in tables from KKY01 and KK06 and match halpha and r filter, or create a dictionary
    
    pass



def get_filter_props(filter):
    """ get filter center and width  """
    filter_wave_dwave = {'Halpha1':[6563,80],\
                             'Halpha2':[6608,76],\
                             'Halpha3':[6573,68],\
                             'Halpha4':[6618,74],\
                             'Halpha5':[6563,78],\
                             'Halpha6':[6606,75],\
                             'R':[6425,1540],\
                             'nmR':[6470,1110],\
                             'sR':[7024,380]
                             }
    return filter_wave_dwave[filter]

def get_zp(galname):
    """  get PHOTZP based on flux ZP and filter width """
    # CONSTANTS
    c = 3e10 # speed of light in cm/s
    f0 = 1e-18 # flux zp in erg/s/cm^2

    ##################################################
    # get filters
    ##################################################    
    rfilter, hafilter = get_filters(galname)

    ##################################################
    # get halpha ZP
    ##################################################    
    # get filter width
    hcenter_A, hwidth_A = get_filter_props(hafilter)
    hcenter_cm = hcenter_A * 1e-8
    hwidth_cm = hwidth_A * 1e-8    
    # calc ZP
    HZP = f0 * c * hwidth_cm/hcenter_cm**2

    ##################################################
    # get r-band ZP
    ##################################################    
    # get filter width
    rcenter_A, rwidth_A = get_filter_props(rfilter)
    rcenter_cm = rcenter_A * 1e-8
    rwidth_cm = rwidth_A * 1e-8    
    # calc ZP
    RZP = f0 * c * rwidth_cm/rcenter_cm**2

    return RZP, HZP

    
if __name__ == '__main__':

    topdir = os.getcwd()
    dirlist = os.listdir()
    dirlist.sort()
    for d in dirlist: # loop through directories
        if d.startswith('sn'): # not sure what these are - might be SN observations?
            continue
        if os.path.isdir(d):

            # get list of files
            # construct galname
            # n4178 -> NGC4178
            # ic3392 -> IC3398
            if d.startswith('n'):
                galname = 'NGC'+d[1:]
            elif d.startswith('ic'):
                galname = 'IC'+d[2:]
            elif d.startswith('u'):
                galname = 'UGC'+d[3:]
            elif d.startswith('i'): # one duplicate observation in isolated sample has "i" instead of "ic"
                galname = 'IC'+d[1:]

            # weird inconsistency
            if galname == 'NGC4411':
                galname = 'NGC4411B'
            print()            
            print("#################################")
            print(f"###  {d}-{galname}   ###########")
            print("#################################")        

            # get info
            ra, dec = get_coords(galname)
            instrument = get_instrument(galname)
            if instrument is None:
                print(f"WARNING: could not find instrument for {galname}, {d}")
                sys.exit()
            pixelScale = get_pixel_scale(instrument)
            pixelScaleDeg = float(pixelScale)/3600 # convert from arcsec/pix to deg/pix

            os.chdir(d) # move to directory            
            filelist = os.listdir() # get list of fits images
            for f in filelist:
                if f.startswith('h'):
                    continue
                if os.path.isfile(f) & ('.fits' in f):
                    hdu = fits.open(f)

                    # get size of image (naxis1, naxis2)                
                    naxis1 = hdu[0].header['NAXIS1']
                    naxis2 = hdu[0].header['NAXIS2']                

                    # build WCS


                    # add wcs info to header
                    hdu[0].header.set('CRVAL1', ra)
                    hdu[0].header.set('CRVAL2', dec)
                    hdu[0].header.set('CTYPE1', 'RA---TAN')
                    hdu[0].header.set('CTYPE2', 'DEC--TAN')


                    hdu[0].header.set('CRPIX1', naxis1//2)
                    hdu[0].header.set('CRPIX2', naxis2//2)


                    hdu[0].header.set('CD1_1', -1*pixelScaleDeg)
                    hdu[0].header.set('CD2_2', pixelScaleDeg)
                    hdu[0].header.set('CD1_2', 0)
                    hdu[0].header.set('CD2_1', 0)                


                    # add instrument
                    hdu[0].header.set('INSTRUME', instrument)

                    # add filter
                    if 'ha' in f:
                        hdu[0].header.set('FILTER', 'Ha4')
                    else:
                        hdu[0].header.set('FILTER', 'R')


                    # prepend an 'h' to image name and save
                    outfile = 'h'+f
                    print(f"\tWriting {f} -> {outfile}")
                    fits.writeto(outfile, hdu[0].data, hdu[0].header, overwrite=True)



            os.chdir(topdir)
            #break
