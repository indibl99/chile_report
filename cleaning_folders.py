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

rootdir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/'
finaldir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/NII/AGN_final'
noisedir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/tooNoisy/'
#outdir = '/Users/plira/india_temp/variable_spectra/auto/output_spectra/'
outdir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/updated_spectra/'
indir = '/Users/plira/india_temp/variable_spectra/auto/input_spectra/'
cumulativeFile = '/Users/plira/india_temp/variable_spectra/auto/cumulative_final.txt'


def listAll(subname, finalArray):
    finalArray = []
    for subdir in os.listdir(rootdir + subname):
        if subdir=='.DS_Store':
            continue
        else:
            finalArray.append(subdir)

    print("total " + subname + ": " + str(len(finalArray)))
    return finalArray

def parseAnalysis(loc):
    file = outdir + loc + '/analysis_output'
    output = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    
    with open(file) as analysis:
        end = 11
        for i, line in enumerate(analysis):
            temp = line.strip().split()
            #print(temp)
            if i==0:
                if len(temp)==3: #might be an updated file
                    output[0] = temp[2]
                elif len(temp)==2:
                    output[0] = temp[1]
            if i==3:
                output[1] = temp[0] #hA
                output[2] = temp[1] #NII
                output[3] = temp[2] #hB
                output[4] = temp[3] #OIII
                output[5] = temp[4] #SII
                output[6] = temp[5] #OI
                output[7] = temp[6] #noise
            if i==9:
                output[14] = temp[0] #NII BPT
                output[15] = temp[1] #SII BPT
                output[16] = temp[2] #OI BPT
            if i==10 and temp==[]: #is there an extra line
                #print('extra line')
                end = 12
            if i==end:
                output[8] = temp[0] #hAClear
                output[9] = temp[1] #NIIClear
                output[10] = temp[2] #hBClear
                output[11] = temp[3] #OIIIClear
                output[12] = temp[4] #SIIClear
                output[13] = temp[5] #OIClear
                
    #print(output)
    #base, emLines.hAlpha, emLines.NII, emLines.hBeta, emLines.OIII, emLines.SII, emLines.OI, emLines.noise,
    #emLines.hAisClear, emLines.NIIisClear, emLines.hBisClear, emLines.OIIIisClear, emLines.SIIisClear, emLines.OIisClear, emLines.tooNoisy
     
    update = open(cumulativeFile, 'r')
    lines = update.readlines()
    update.close()
    update = open(cumulativeFile, 'w')
    for line in lines:
        params = line.split(',')
        if params[0]!=output[0]:
            update.write(line)
    update.write(str(output[0]) + ',' + str(output[1]) + ',' + str(output[2]) + ',' + str(output[3]) + ',' + str(output[4]) + ',' + str(output[5]) + ',' + str(output[6]) + ',' + str(output[7]) + ',' + str(output[8]) +',' + str(output[9]) +',' + str(output[10]) + ',' + str(output[11]) + ',' + str(output[12]) + ',' + str(output[13]) + ',' + str(output[14]) + ',' + str(output[15]) + ',' + str(output[16]) + '\n')
    update.close()

    if output[14]=='AGN' and output[15]=='Seyfert' and output[16]=='Seyfert':
        copy_tree(outdir + loc, outdir + 'likely_AGN_final/' + loc)
        
    

#318 total input spectra 
updateDir = []
noisySpec = []  
NII_AGN = []
NII_HII = []
NII_Composite = []
SII_Seyfert = []
SII_HII = []
SII_LINER = []
OI_Seyfert = []
OI_HII = []
OI_LINER = []


NII_AGN = listAll('NII/AGN_final', NII_AGN)
NII_HII = listAll('NII/HII_final', NII_HII)
NII_Composite = listAll('NII/Composite_final', NII_Composite)
NII = NII_AGN
print('\n')

SII_Seyfert = listAll('SII/Seyfert_final', SII_Seyfert)
SII_HII = listAll('SII/HII_final', SII_HII)
SII_LINER = listAll('SII/LINER_final', SII_LINER)
SII = SII_Seyfert
print('\n')

OI_Seyfert = listAll('OI/Seyfert_final', OI_Seyfert)
OI_HII = listAll('OI/HII_final', OI_HII)
OI_LINER = listAll('OI/LINER_final', OI_LINER)
OI = OI_Seyfert


updateDir = listAll('updated_spectra', updateDir)
noisySpec = listAll('tooNoisy', noisySpec)

for spec in updateDir:
    parseAnalysis(spec)

NII.extend(NII_HII)
NII.extend(NII_Composite)
SII.extend(SII_HII)
SII.extend(SII_LINER)
OI.extend(OI_HII)
OI.extend(OI_LINER)

print('\n')
print("total out: " + str(len(updateDir) + len(noisySpec)))
print("total NII: " + str(len(NII)))
print("total SII: " + str(len(SII)))
print("total OI: " + str(len(OI)))

"""for item in OI:
    if item not in updateDir:
        print("look!")
        print(item)"""


#uncomment and edit to find duplicates

uniqueDir = np.copy(OI)
copyUpdate = np.copy(updateDir)
new = np.append(uniqueDir, copyUpdate)
i = 1
print(len(new))

for spec in uniqueDir:
    seen = 0
    for checkAgain in new:
        if checkAgain==spec:
            seen = seen + 1
    #print(seen)
    if seen>2:
        print("found something!")
        print(spec)
        #print(seen)
    i = i + 1

#temp = updateDir
#totalDir = temp
#totalDir.extend(noisySpec)

#parseAnalysis('spec-0267-51608-0493')
