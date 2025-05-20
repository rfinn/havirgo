#!/usr/bin/env python

"""
conscale_factors.txt has 867 lines, whereas cutout directory has 868 subdirectories.

want to find out which one is missing and why

run this from cutouts directory

/data-pool/Halpha/halphagui-output-20240522/cutouts
"""
import os
infile = open('conscale_factors.txt','r')
cdirnames = []
for line in infile:
    t = line.split()
    # first column is the dirname
    cdirnames.append(t[0])
infile.close()

# now scan through directories

dirlist = os.listdir()

for d in dirlist:
    if os.path.isdir(d):
        if d in cdirnames:
            continue
        else:
            print(f"WARNING:{d} is not in conscale_factors.txt")
