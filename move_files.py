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


rootdir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications/NII/AGN'
outdir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications/updated_spectra/'

#this works! DO NOT DELETE
"""for subdir in os.listdir(rootdir):
    if subdir=='.DS_Store':
        continue
    elif subdir[-8:]!='redacted':
        print(subdir)
        copy_tree(rootdir + '/' + subdir, outdir + '/' + subdir + '/')"""

totalOut = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_final.txt', 'w')
out = '/Users/plira/india_temp/variable_spectra/auto/cumulative_final.txt'

for subdir in os.listdir(outdir):
    #print(subdir)
    if subdir=='.DS_Store':
        continue
    else:
        #for filename in os.listdir(outdir + subdir):
        #print(outdir + subdir + '/analysis_output')
        analysis = open(outdir + subdir + '/analysis_output', 'r').read()
        file = outdir + subdir + '/analysis_output'
        #os.rename(file, file.replace(" ", "_"))
        lines = analysis.split('\n')
        #print(lines)
        #print(subdir[-6:])
        if subdir[-6:]=='update':
            end = 11
        else:
            end = 12
        i = 0
        readin = []
        for line in lines:
            if i==3:
                readin = line.strip().split()
                #print(readin)
            elif i==end:
                readin.extend(line.strip().split())
                #print(readin)
            i = i + 1
        print(len(readin))
        if len(readin) > 0:
            temp = ",".join(readin)
            totalOut.write(temp + '\n')
            #totalOut.write('\n')
                
                
            

                
                


    


    
