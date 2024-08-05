#!/usr/bin/env python

"""
GOAL:
* to help speed up review of continuum-subtracted image

PROCEDURE:
* run subtract_continuum.py multiple times until the user is happy with scale factor
* display r, halpha, and CS in ds9 for each time

"""

import os
import sys
import numpy
import glob
from matplotlib import pyplot as plt


def display_results(d):
    rimage = os.path.join(d,d+'-R.fits')
    haimage = os.path.join(d,d+'-CS-gr.fits')
    images = [rimage,haimage]
    for i,im in enumerate(images):
        s = f"xpaset -p ds9 frame {i+1}"
        #print(s)
        os.system(s)
        
        s = f"xpaset -p ds9 fits {im}"
        os.system(s)
        os.system("xpaset -p ds9 scale zscale")
        #os.system("xpaset -p ds9 zoom to fit")        
        

if __name__ == '__main__':
    os.chdir('/media/rfinn/hdata2/halphagui-output-20240522/cutouts/')
    dirlist = os.listdir()

    
    ngal = 0
    for d in dirlist:
        if d.startswith('VFID') and os.path.isdir(d):
            ngal += 1
            
    print(f"only {ngal} galaxies to go!")
    for j,d in enumerate(dirlist):
        if d.startswith('VFID') and os.path.isdir(d):
            print()
            print('############################')
            print(f"{d} ({j}/{ngal})")
            print('############################')
            display_results(d)
            scale = 1.1
            adjust_scale = True
            finished = False
            while adjust_scale:
                continue_query = True
                while continue_query:
                    t = input("enter a scale factor or q to quit this galaxy, and x to exit\n")
                    print(f"you entered {t}")
                    if t.find('q') > -1:
                        continue_query = False
                        adjust_scale = False
                        print()
                        print(f"{d}: scale = {scale:.2f}")
                        print()
                        finished = True
                    elif t.find('x') > -1:
                        print("exiting program")
                        sys.exit()
                    else:
                        try:
                            scale = float(t)
                            continue_query = False
                        except ValueError:
                            print("did not understand input, please try again")
                            print()
                if finished:
                    # TODO : move directory to FINISHED_CS
                    # print out current scale values
                    os.system(f"gethead CONSCALE {d}/{d}-CS-gr*.fits")
                    
                    print(f"moving {d} to FINISHED_CS")
                    os.system(f"mv {d} FINISHED_CS/.")
                    c = input("press any key to continue\n")
                    pass
                else:
                    command_string = f"python ~/github/havirgo/python/subtract_continuum.py {d} {scale}"
                    #print()
                    #print(command_string)
                    #print()
                    os.system(command_string)
                    plt.show()
                    display_results(d)

                    # print out current scale values
                    os.system(f"gethead CONSCALE {d}/{d}-CS-gr*.fits")

            
            

                
