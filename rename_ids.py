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
#import xlrd


#id_key = open('/Users/plira/india_temp/variable_spectra/auto/ids_key.txt', 'r').read()
#lines = id_key.split('\n')

id_key = '/Users/plira/india_temp/variable_spectra/auto/id_key_final.txt'
renameDir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/spectral_classifications_final_renamed_internal_ID/SII'
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
    numSpecs = 0
    for i in xrange(0,len(temp),4):
        #thisKey = [temp[i], temp[i+1]]
        #print(temp[i])
        #print(temp[i+1])
        spec = 'spec-' + temp[i] + '-' + temp[i+1] + '-' + temp[i+2]
        key[0].append(spec)
        key[1].append(temp[i+3])
        numSpecs = numSpecs + 1
        print('numSpecs: ' + str(numSpecs))
print(key)

i = 0
for filename in os.listdir(renameDir):
    if filename=='.DS_Store':
        continue
    else:
        #print(filename)
        if filename in key[0]:
            ind = key[0].index(filename)
            print("changing " + filename + " to " + key[1][ind] + '-' + filename)
            os.rename(renameDir + '/' + filename, renameDir + '/' + key[1][ind] + '-' + filename)
            #print('number specs in dir: ' + str(i))
            i = i + 1

print('num ids: ' + str(len(key[0])))

