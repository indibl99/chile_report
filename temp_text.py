#!/usr/bin/python 
import os
from os import path
import sys
import numpy as np
import scipy.optimize as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import shutil
from shutil import copyfile
from collections import namedtuple
from peakutils import indexes
from subprocess import call
#import in_place
import fileinput
import csv
import warnings
import pickle as pl

totalOut = open('/Users/plira/india_temp/variable_spectra/auto/temp_output.txt','w')
output = ['13', '2304j', 'sfbgi3', 'sdg']
totalOut.write("{: <15} {: <15} {: <15} {: <15}".format('one', 'two', 'three', 'four'))
totalOut.write("\n{: <15} {: <15} {: <15} {: <15}".format(output[0], output[1], output[2], output[3]))
totalOut.close()


"""
def cumulative_BPT(totalFile):
    total = open(totalFile, 'r')
    for i, line in enumerate(total):
        if i==0:
            continue
        else:
            temp1 = line.strip().split()
            cBPT = namedtuple('cBPT', 'hAlpha NII hBeta OIII SII OI noise hAisClear NIIisClear hBisClear OIIIisClear SIIisClear OIisClear tooNoisy')
            BPTparams = cBPT(float(temp1[1]), float(temp1[2]), float(temp1[3]), float(temp1[4]), float(temp1[5]), float(temp1[6]), float(temp1[7]), float(temp1[8]), float(temp1[9]), float(temp1[10]), float(temp1[11]), float(temp1[12]), float(temp1[13]), float(temp1[14]))
            if BPTparams.tooNoisy==True:
                print("tooNoisy, skipping: " + str(i))
            elif BPTparams.tooNoisy==False:
                BPT(BPTparams, '/Users/plira/india_temp/variable_spectra/auto/dummy_out.txt', True)
                c, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
                ax1.plot( np.log(BPTparams.NII/BPTparams.hAlpha), np.log(BPTparams.OIII/BPTparams.hBeta) )
                ax1.set_title('log(NII/H-alpha) vs. log(OIII/H-beta'))
                ax2.plot( np.log(BPTparams.SII/BPTparams.hAlpha), np.log(BPTparams.OIII/BPTparams.hBeta) )
                ax2.set_title('log(SII/H-alpha) vs. log(OIII/H-beta'))                
                ax3.plot( np.log(BPTparams.OI/BPTparams.hAlpha), np.log(BPTparams.OIII/BPTparams.hBeta) )
                ax3.set_title('log(OI/H-alpha) vs. log(OIII/H-beta'))"""
