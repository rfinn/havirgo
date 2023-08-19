#!/usr/bin/env python

'''
GOAL:
* create index web page that points to html pages for all cutouts 

USAGE:
* run from html-dev/cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)

'''

import os
import numpy as np
import glob

from astropy.io import fits
homedir = os.getenv("HOME")

###########################################################
####  FUNCTIONS
###########################################################

    
def write_coadd_prop_table(html,filter,zp,fwhm_arcsec):
    html.write('<h3>Image Characteristics</h3>\n')
    html.write('<table>\n')
    html.write('<tr>')
    html.write('<th ">Filter</th>\n')
    html.write('<th ">ZP<br />mag</th>\n')
    html.write('<th ">PSF FWHM <br />arcsec</th>\n')
    html.write('</tr>\n')
    html.write('<tr><td>{}</td><td>{:.2f}</td><td>{:.2f}</td>\n'.format(filter,zp,fwhm_arcsec))
    html.write('</tr>\n')
    html.write('</table>\n')

def write_table(html,images,labels):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for i in images:
        html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
    html.write('</tr>\n')            
    html.write('</table>\n')

def write_text_table(html,labels,data):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for d in data:
        html.write('<td>{}</td>'.format(d))
    html.write('</tr>\n')            
    html.write('</table>\n')


def legacy_link(ra,dec):
    return "https://www.legacysurvey.org/viewer?ra={:.4f}&dec={:.4f}&layer=ls-dr9&zoom=13".format(ra,dec)    
###########################################################
####  CLASSES
###########################################################
    
class check_cutouts():

    def __init__(self,gallist,outdir,co=False):
        ''' pass in instance of cutout_dir class and output directory '''

        self.coflag = co
        if self.coflag:
            outfile = os.path.join(outdir,'indexco.html')
        else:
            outfile = os.path.join(outdir,'index.html')
        self.outdir = outdir

        #print('inside build html')
        #print('coutdir = ',coutdir)
        #print('outfile = ',outfile)        
        self.galnames = gallist
        self.check_gal_list()
    def check_gal_list(self):
        for i,g in enumerate(self.galnames):
            #print(g)
            vfid = g.split('-')[0]
            vfindex = vfindices[vfmain['VFID'] == vfid][0]
            ra = vfmain['RA'][vfindex]
            dec = vfmain['DEC'][vfindex]
            if self.coflag & ~vfmain['COflag'][vfindex]:
                continue
            ids.append(vfmain['VFID'][vfindex])
            #print(vfindex)
            # get legacy jpg name

            jpg_path = os.path.join(self.outdir,g)
            search_path = os.path.join(jpg_path,'*legacy*.jpg')
            #print(search_path)
            #legacy_jpg = glob.glob(search_path)[0]            
            try:
                #print()
                #print("looking for legacy image")
                #print(glob.glob(search_path))
                legacy_jpg = glob.glob(search_path)[0]
                legacy_flag = True                
            except:
                legacy_flag = False
                legacy_jpg = None
                print('WARNING: no legacy image for ',g)
                print("trying to download")
                #print('\t Skipping galaxy for now')
                #continue

                # try to rebuild website
                # download cutouts
                cutout_dir = g
                current_dir = os.getcwd()

                # move to cutouts data dir
                os.chdir("/data-pool/Halpha/halphagui-output-20230818/cutouts/")
                s = f"python ~/github/HalphaImaging/python3/generate_all_cutout_plots.py --onegal {cutout_dir}"
                os.system(s)

                # build webpage
                s = f"python ~/github/havirgo/python/build_web_cutouts2.py --oneimage {cutout_dir}"
                os.system(s)

                os.chdir(current_dir)
                try:
                    legacy_jpg = glob.glob(search_path)[0]
                    legacy_flag = True                
                except:
                    legacy_flag = False
                    legacy_jpg = None
                    print('ERROR: no legacy image for ',g)
                    

if __name__ == '__main__':
    # work through coadd directory
    #global vmain

    VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
    VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'    
    vfmain = fits.getdata(VFMAIN_PATH)

    # updating for v2 catalogs
    VFFIL_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main_filament_membership_allgalaxies.fits'
    VFFIL_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_environment.fits'    
    vffil = fits.getdata(VFFIL_PATH)

    VFMAGPHYS_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_magphys_10-Jul-2023.fits'    
    vfmagphys = fits.getdata(VFMAGPHYS_PATH)
    
    outdir = homedir+'/research/Virgo/html-dev/cutouts/'
    outdir = '/data-pool/Halpha/html_dev/cutouts/'    

    # this should contain a list of all the galaxy folders
    flist1 = os.listdir(outdir)
    flist1.sort()
    galnames=[]
    for i,subdir in enumerate(flist1): # loop through list
        

        #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
        if (os.path.isdir(subdir)) & (subdir.startswith('VF')):
            #print('adding ',subdir)
            galnames.append(subdir)
    print('number of subdirectories = ',len(galnames))
    h = check_cutouts(galnames,outdir)
