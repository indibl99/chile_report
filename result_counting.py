#!/usr/bin/python 

import os
from os import path
import sys
#print '\n'.join(sys.path)

import numpy as np
import scipy.optimize as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys
#from astroML.datasets.tools.sdss_fits import log_OIII_Hb_NII
#from astroML.datasets.tools.sdss_fits import log_OIII_Hb_SII
#from astroML.datasets.tools.sdss_fits import log_OIII_Hb_OI
import shutil
from shutil import copyfile
from collections import namedtuple
from peakutils import indexes
from subprocess import call
import fileinput
import warnings
import pickle as pl
from distutils.dir_util import copy_tree


rootdir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/NII/AGN'
finaldir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/NII/AGN_final'
noisedir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/tooNoisy/'
#outdir = '/Users/plira/india_temp/variable_spectra/auto/output_spectra/'
outdir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/updated_spectra/'
indir = '/Users/plira/india_temp/variable_spectra/auto/input_spectra/'


#this works! DO NOT DELETE
"""for subdir in os.listdir(rootdir):
    if subdir=='.DS_Store':
        continue
    elif subdir[-6:]=='update':
        print(subdir)
        copy_tree(rootdir + '/' + subdir, finaldir + '/' + subdir[:-7] + '/')
    elif subdir[-8:]=='redacted':
        print(subdir)
    else:
        copy_tree(rootdir + '/' + subdir, finaldir + '/' + subdir + '/')"""

#totalOut = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_tooNoisy.txt', 'w')
#out = '/Users/plira/india_temp/variable_spectra/auto/cumulative_tooNoisy.txt'
updateDir = []
numDirs = 0
inDir = []
inDirs = 0

for subdir in os.listdir(indir):
    #print(subdir)
    if subdir=='.DS_Store':
        continue
    else:
        if subdir[-8:]=='redacted':
            print("redacted: " + subdir)
            continue
        else:
            inDir.append(subdir)

        inDirs = inDirs + 1
        #print(numDirs)
        
for subdir in os.listdir(outdir):
    #print(subdir)
    if subdir=='.DS_Store':
        continue
    else:
        if subdir[-8:]=='redacted':
            print("redacted: " + subdir)
            continue
        if subdir[-6:]=='update':
            updateDir.append(subdir[:-7])
        else:
            updateDir.append(subdir)

noisySpec = []        
for subdir in os.listdir(noisedir):
    #print(subdir)
    if subdir=='.DS_Store':
        continue
    else:
        if subdir[-8:]=='redacted':
            print("redacted: " + subdir)
            continue
        if subdir[-6:]=='update':
            noisySpec.append(subdir[:-7])
        else:
            #updateDir.append(subdir)
            noisySpec.append(subdir)

print("total: " + str(len(inDir)))
print("total: " + str(len(updateDir)))
print("total: " + str(len(noisySpec)))


updateDir.extend(noisySpec)
uniqueDir = updateDir

for item in uniqueDir:
    print(item)


"""
#print(str(np.shape(updateDir)))      
for elem in updateDir:
    #if elem[-6:]=='update':
    if 'update' in elem:
        #print(elem)
        updateDir.remove(elem)
        #elem = elem[:-7]
        updateDir.append(elem[:-7])
        #print(elem)
        #print(elem)
#print("updateDir: " + str(len(updateDir)))
#print(updateDir)
#uniqueDir = updateDir
#uniqueDir = list(set(uniqueDir))
seen = set()
uniqueDir = []

for item in updateDir:
    #print(item + '\n')
    if item not in seen:
        seen.add(item)
        uniqueDir.append(item)
"""
print("total: " + str(len(uniqueDir)))

    

numDup = 0

"""Broadline AGNs
spec-0270-51909-0050 (check!)

uniqueDir = updateDir
i = 1

for item in uniqueDir:
    #print(item)
    #if uniqueDir[i]==uniqueDir[i-1]:
    seen = 0
    for checkAgain in uniqueDir:
        if checkAgain==item:
            seen = seen + 1
    if seen >=2:
        print("found duplicate between noisy and update dir!")
        print(item)
    i = i + 1
"""    
