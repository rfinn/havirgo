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
* query NED using NGC/obj name
* get RA and DEC
* get pixelscale

"""

import argparse


# give user option to input pixelscale


# get OBJECT from image header ?


# convert to NED friendly name

# get coords from NED

# build WCS  
