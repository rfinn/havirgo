# sort files into separate subdirectories

import glob
import os

# get list of files

filelist = os.listdir()
# if name starts with n, take first 5 characters
# if name starts with ic, take first 6 characters
# if name starts with sn, take first 6 characters

allprefix = []
for fname in filelist:
    #print(fname)
    if fname.find('.fits') > -1:
        if fname.startswith('n'):
            allprefix.append(fname[:5])
        elif fname.startswith('ic'):
            allprefix.append(fname[:6])
        elif fname.startswith('sn'):
            allprefix.append(fname[:6])
        elif fname.startswith('u'):
            allprefix.append(fname[:5])
        elif fname.startswith('i'):
            allprefix.append(fname[:5])

#print(allprefix)
objlist = set(allprefix)
print(f"number of unique objects = {len(objlist)}")

for gal in objlist:
    # create subdirectories
    if not os.path.exists(gal):
        os.mkdir(gal)
    for fname in filelist:
        if os.path.isdir(fname):
            continue
        if fname.startswith(gal):
            # move images
            print(os.path.join(gal,fname))
            os.rename(fname,os.path.join(gal,fname))
        
