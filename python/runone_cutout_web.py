#!/usr/bin/env python

"""
GOAL:
This program will:
* download the cutout images
* build a cutout webpage
* remake the cutout index page


RATIONALE:
for whatever reason, getting the cutout images fails in many cases.
Therefore, when I make a new version of the cutouts/data, the most time-consuming
part is having to go back and get images for the cutouts which failed.

This script will streamline the process in one program.

NOTE: 

I thought I already had a program to do this, but I can't find it.

USAGE:
Run from the cutouts directory, like: /data-pool/Halpha/halphagui-output-20230814/cutouts/

To run, type:

python ~/github/havirgo/runone_cutout_web.py cutout_dirname

Using a specific directory name:

python ~/github/havirgo/runone_cutout_web.py cutout_dirname VFID0483-NGC6306-BOK-20220424-VFID0483 

"""

import sys
import os

startdir = os.getcwd()
cutout_dir = sys.argv[1]

if len(sys.argv) > 2:
    syncfiles=False
else:
    syncfiles=True

# download cutouts
s = f"python ~/github/HalphaImaging/python3/generate_all_cutout_plots.py --onegal {cutout_dir}"
os.system(s)

# build webpage
s = f"python ~/github/havirgo/python/build_web_cutouts2.py --oneimage {cutout_dir}"
os.system(s)


# build cutout index - I am assuming this is running on draco
os.chdir('/data-pool/Halpha/html_dev/cutouts/')
os.system("python ~/github/havirgo/python/build_cutout_index.py")

# sync to web
if syncfiles:
    print()
    print("syncing to facultyweb - enter Siena password when prompted...")
    print()
    os.chdir("/data-pool/Halpha/html_dev/")
    os.system("rsync -avz cutouts facultyweb.siena.edu:public_html/virgo/. ")

# return to original directory
os.chdir(startdir)
