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

id_key = '/Users/plira/india_temp/variable_spectra/auto/stacking/luminosities.txt'
luminosityDir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/spectral_classifications_final_renamed_internal_ID/updated_spectra'

# Creates a list containing 2 lists, each of 1 items, all set to 0
w, h = 0, 2;
key = [[0 for x in range(w)] for y in range(h)]


graph_data = open(id_key, 'r').read()
lines = graph_data.split('\n')
for line in lines:
    temp = line.split()
    #print(temp)
    #print('\n')
    print("number array elements " + str(len(temp)))
    upperLim = len(temp) / 4
    print('upperLim: ' + str(upperLim))
    numSpecs = 1
    for i in xrange(0,len(temp),2):
        key[0].append(temp[i])
        key[1].append(temp[i+1])
        numSpecs = numSpecs + 1
#print('numSpecs: ' + str(numSpecs))

specDir = []
mags = []
#bin1 = []
#bin2 = []
#bin3 = []
bin1 = [[0 for x in range(0)] for y in range(2)]
bin2 = [[0 for x in range(0)] for y in range(2)]
bin3 = [[0 for x in range(0)] for y in range(2)]

for gmag in key[1]:
    mags.append(float(gmag))

numSpecs = 0
for filename in os.listdir(luminosityDir):
    if filename=='.DS_Store':
        continue
    else:
        #print(filename[:-21])
        specDir.append(filename)
        numSpecs = numSpecs + 1
print('numSpecs: ' + str(numSpecs))

spec2 = []
mags2 = []
g1 = 0
g2 = 0
g3 = 0
for spec in specDir:
    if spec[:-21] in key[0]:
        ind = key[0].index(spec[:-21])
        #print('luminosity of ' + str(key[0][ind]) + ' is ' + str(key[1][ind]))
        spec2.append(key[0][ind])
        mags2.append(float(key[1][ind]))
        
        if np.abs(float(key[1][ind])) > 20.0:
            g3 = g3 + 1
            bin3[1].append(float(key[1][ind]))
            bin3[0].append(spec)
        elif np.abs(float(key[1][ind])) <= 20.0 and np.abs(float(key[1][ind])) > 19:
            g2 = g2 + 1
            bin2[1].append(float(key[1][ind]))
            bin2[0].append(spec)
        else:
            g1 = g1 + 1
            bin1[1].append(float(key[1][ind]))
            bin1[0].append(spec)
            
#print('max magnitude: ' + str(np.amax(mags2)))
#print('min magnitude: ' + str(np.amin(mags2)))
#print('numSpecs: ' + str(len(mags2)))
print('bin1 : ' + str(g1))
print('bin2 : ' + str(g2))
print('bin3 : ' + str(g3))

#print(bin1[0])
bins = [bin1, bin2, bin3]

bin1Dir = '/Users/plira/india_temp/variable_spectra/auto/stacking/luminosity_bin1'
bin2Dir = '/Users/plira/india_temp/variable_spectra/auto/stacking/luminosity_bin2'
bin3Dir = '/Users/plira/india_temp/variable_spectra/auto/stacking/luminosity_bin3'
binDirs = [bin1Dir, bin2Dir, bin3Dir]
specDir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/spectral_classifications_final_renamed_internal_ID/updated_spectra'

i = 0
for current_bin in bins:
    binDir = binDirs[i]
    if os.path.exists(binDir):
        shutil.rmtree(binDir)  
    os.makedirs(binDir)
    
    for spec in current_bin[0]:
        copy_tree(specDir + '/' + spec, binDir + '/' + spec)

    i = i + 1
        


#n, bins, patches = plt.hist(mags2)
#plt.show()

#exactly 1/2 are less than 20 


#greater than 20
#19-20
#16-19






