#!/usr/bin/python 

import os
from os import path
import sys

import numpy as np
import scipy.optimize as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys
import shutil
from shutil import copyfile
from collections import namedtuple
from peakutils import indexes
from subprocess import call
import fileinput
import warnings
import pickle as pl
from distutils.dir_util import copy_tree

#fig = plt.figure()
c, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 10))#3 subplots
for ax in [ax1, ax2, ax3]:
    ax.set_ylim(-3, 3)
    ax.set_xlim(-5, 3)
    ax.grid()

#NII BPT

    
xVals2 = np.arange(-10, 0.05, 0.01)
xVals = np.arange(-10, 0.47, 0.01)
nDelim2 = 0.61 / (xVals2 - 0.05) + 1.3
nDelim = 0.61 / (xVals - 0.47) + 1.19
ax1.set_title("BPT Diagram Using NII")
ax1.set_xlabel("[NII]/H-alpha")
ax1.set_ylabel("[OIII]/H-beta")
ax1.text(-0.6, -0.4, "HII", color='b')
ax1.text(-0.3, 0, "Composite", color='b')
ax1.text(0, 0.6, "AGN", color='b')
ax1.plot(xVals, nDelim, c='k', dashes=[6,2])
ax1.plot(xVals2, nDelim2, c='k')

#SII BPT

xS = np.arange(-6, 0.32, 0.03)
xS2 = np.arange(-0.315, 4, 0.01)
sDelim = 0.72 / (xS - 0.32) + 1.3
sDelim2 = 1.89*xS2 + 0.76 
ax2.set_title("BPT Diagram Using SII")
ax2.set_xlabel("[SII]/H-alpha")
ax2.set_ylabel("[OIII]/H-beta")
ax2.text(-0.4, -0.4, "HII", color='b')
ax2.text(-0.5, 0.6, "Seyfert", color='b')
ax2.text(0.2, -0.2, "LINER", color='b')
ax2.plot(xS, sDelim, c='k', dashes=[6,2])
ax2.plot(xS2, sDelim2, c='k')

#OI BPT 

xO = np.arange(-6, -0.59, 0.01)
xO2 = np.arange(-1.13, 4, 0.01)
oDelim = 0.73 / (xO + 0.59) + 1.33
oDelim2 = 1.18*xO2 + 1.3
ax3.set_title("BPT Diagram Using OI")
ax3.set_xlabel("[OI]/H-alpha")
ax3.set_ylabel("[OIII]/H-beta")
ax3.text(-1.6, -0.4, "HII", color='b')
ax3.text(-1, 0.6, "Seyfert", color='b')
ax3.text(-0.5, -0.2, "LINER", color='b')
ax3.plot(xO, oDelim, c='k', dashes=[6,2])
ax3.plot(xO2, oDelim2, c='k')

graph_data = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_final.txt', 'r').read()
#lines = graph_data.strip()
lines = graph_data.split('\n')
xNIIc = []
yNIIc = []
xSIIc = []
ySIIc = []
xOIc = []
yOIc = []

output = []
noise = []

hAc = []
NIIc = []
hBc = []
OIIIc = []
SIIc = []
OIc = []


NIIcolors = []
SIIcolors = []
OIcolors = []

for line in lines:
    if len(line) > 1:
        #'Filename','hAlpha','NII','hBeta','OIII','SII','OI','noise','hAisClear','NIIisClear','hBisClear','OIIIisClear','SIIisClear','OIisClear','tooNoisy', 'BPT NII', 'BPT SII', 'BPT OI'
        #name, hA, NII, hB, OIII, SII, OI, noise, hAc, NIIc, hBc, OIIIc, SIIc, OIc, tooNoisy, BPTNII, BPTSII, BPTOI = line.split(',')
        readin = line.split(',')
        #print(readin)
        end = len(readin)
        if end==15:
            print('tooNoisy')
            continue
        else:
            outputRow = []

            """hA = readin[1]
            NII = readin[2]
            hB = readin[3]
            OIII = readin[4]
            SII = readin[5]
            OI = readin[6]
            noise = readin[7]"""

            outputRow.append(readin[0]) #hA
            outputRow.append(readin[1]) #NII
            outputRow.append(readin[2]) #hB
            outputRow.append(readin[3]) #OIII
            outputRow.append(readin[4]) #SII
            outputRow.append(readin[5]) #OI
            outputRow.append(readin[6]) #noise
            outputRow.append(readin[7]) #hAClear
            outputRow.append(readin[8]) #NIIClear
            outputRow.append(readin[9]) #hBClear
            outputRow.append(readin[10]) #OIIIClear
            outputRow.append(readin[11]) #SIIClear
            outputRow.append(readin[12]) #OIClear
            

            """xNII = np.log(float(NII)/float(hA))
            yNII = np.log(float(OIII)/float(hB))
            outputRow.append(xNII)
            outputRow.append(yNII)
            #xNIIc.append(xNII)
            #yNIIc.append(yNII)
            #print(xNIIc)
            
            xSII = np.log(float(SII)/float(hA))
            ySII = np.log(float(OIII)/float(hB))
            outputRow.append(xSII)
            outputRow.append(ySII)
            #xSIIc.append(xSII)
            #ySIIc.append(ySII)
            #print(xSIIc)
            
            xOI = np.log(float(OI)/float(hA))
            yOI = np.log(float(OIII)/float(hB))
            outputRow.append(xOI)
            outputRow.append(yOI)
            #xOIc.append(xOI)
            #yOIc.append(yOI)
            #print(xOIc)"""

            """for i in range(7,end):
                if readin[i]=='True':
                    outputRow.append(1)
                elif readin[i]=='False':
                    outputRow.append(0)
            #print(outputRow)"""
            output.append(outputRow)

            """hAclear = output[8]
            NIIclear = output[9]
            hBclear = output[10]
            OIIIclear = output[11]
            SIIclear = output[12]
            OIclear = output[13]
            tooNoisy = output[14"""
            #BPTNII = output[15]
            #BPTSII = output[16]
            #BPTOI = output[17]
        
    
#print(len(output))
numAGN = 0
for j in range(len(output)):
    #print(j)
    print(output[j])
    #-1 is down, 1 is up
    flagXn = 0 #-1 is pointing to the left, 1 is pointing to the right
    flagXs = 0
    flagXo = 0
    flagY = 0
    flagX = 0
    allClear = True
    allClearNII = True
    allClearSII = True
    allClearOI = True

    Halpha = float(output[j][0])
    NII = float(output[j][1])
    Hbeta = float(output[j][2])
    OIII = float(output[j][3])
    SII = float(output[j][4])
    OI = float(output[j][5])
    noise = float(output[j][6])
    
    if output[j][7]==0: 
        Halpha = noise   
        flagX = 1
        allClear = False
    if output[j][8]==0:
        NII = noise
        flagXn = -1
        allClearNII = False
    if output[j][9]==0:
        Hbeta = noise
        flagY = 1
        allClear = False
    if output[j][10]==0:
        OIII = noise
        flagY = -1
        allClear = False
    if output[j][11]==0:
        SII = noise*2
        flagXs = -1
        allClearSII = False
    if output[j][12]==0:
        OI = noise
        flagXo = -1
        OIClear = False
    
    #color-code dot based on region
    colorVarNII = 'k'
    #print(Halpha)
    delim1NII = 0.61 / (np.log(NII/Halpha) - 0.47) + 1.19
    delim2NII = 0.61 / (np.log(NII/Halpha) - 0.05) + 1.3
    if np.log(OIII/Hbeta) >= delim1NII or np.log(NII/Halpha) >= 0.47:
        colorVarNII = '#66ff33'
        bptNII = 'AGN'
        print('AGN!')
        numAGN = numAGN + 1
        print(numAGN)
        
    elif np.log(OIII/Hbeta) >= delim2NII and np.log(OIII/Hbeta) < delim1NII:
        colorVarNII = '#ffff00'
        bptNII = 'Composite'
    elif np.log(OIII/Hbeta) < delim2NII:
        colorVarNII = 'r'
        bptNII = 'HII'

    ax1.plot(np.log(NII/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarNII, markeredgecolor=colorVarNII, alpha=0.5)
    if allClear==False or allClearNII==False: ax1.quiver(np.log(NII/Halpha), np.log(OIII/Hbeta), flagX, flagY,color=colorVarNII) #adds an arrow if uncertain


    #SII BPT

    colorVarSII = 'k'
    delim1SII = 0.72 / (np.log(SII/Halpha) - 0.32) + 1.3
    delim2SII = 1.89*np.log(SII/Halpha) + 0.76
    if np.log(OIII/Hbeta) >= delim1SII and np.log(OIII/Hbeta) >= delim2SII:
        colorVarSII = '#66ff33'
        bptSII = 'Seyfert'
    elif np.log(SII/Halpha) >= 0.32 and np.log(OIII/Hbeta) >= delim2SII:
        colorVarSII = '#66ff33'
        bptSII = 'Seyfert'
    elif np.log(OIII/Hbeta) < delim2SII and np.log(OIII/Hbeta) >= delim1SII:
        colorVarSII = '#ffff00'
        bptSII = 'LINER'
    elif np.log(SII/Halpha) >= 0.32 and np.log(OIII/Hbeta) < delim2SII:
        colorVarSII = '#ffff00'
        bptSII = 'LINER'
    elif np.log(OIII/Hbeta) < delim1SII: #check this
        colorVarSII = 'r'
        bptSII = 'HII'
        
    ax2.plot(np.log(SII/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarSII, markeredgecolor=colorVarSII, alpha=0.5)
    if allClear==False or allClearSII==False: ax2.quiver(np.log(SII/Halpha), np.log(OIII/Hbeta), flagX, flagY, color=colorVarSII)

    #OI BPT

    colorVarOI = 'k'
    delim1OI = 0.73 / (np.log(OI/Halpha) + 0.59) + 1.33
    delim2OI = 1.18*np.log(OI/Halpha) + 1.3

    if np.log(OIII/Hbeta) >= delim1OI and np.log(OIII/Hbeta) >= delim2OI:
        colorVarOI = '#66ff33'
        bptOI = 'Seyfert'
    elif np.log(OI/Halpha) >= -0.59 and np.log(OIII/Hbeta) >= delim2OI:
        colorVarOI = '#66ff33'
        bptOI = 'Seyfert'
    elif np.log(OIII/Hbeta) < delim2OI and np.log(OIII/Hbeta) >= delim1OI:
        colorVarOI = '#ffff00'
        bptOI = 'LINER'
    elif np.log(OI/Halpha) >= -0.59 and np.log(OIII/Hbeta) < delim2SII:
        colorVarOI = '#ffff00'
        bptOI = 'LINER'
    elif np.log(OIII/Hbeta) < delim1OI: #check this
        colorVarOI = 'r'
        bptOI = 'HII'

    ax3.plot(np.log(OI/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarOI, markeredgecolor=colorVarOI, alpha=0.5)
    if allClear==False or allClearOI==False: ax3.quiver(np.log(OI/Halpha), np.log(OIII/Hbeta), flagX, flagY, color=colorVarOI)


plt.tight_layout()
plt.show()
    
