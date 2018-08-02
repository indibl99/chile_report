#!/usr/bin/python 


import os
from os import path
import sys
#print '\n'.join(sys.path)

import numpy as np
import scipy.optimize as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.pyplot import figure

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

#has to be in the loop if you want it to update in real time

total = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_output.txt', 'r')

#creating figure with 2 subplots
c, (ax1, ax2, ax3) = plt.subplots(3, 1)

#plot delimeter equations here and initializing point objects
point1 = ax1.plot([], [], marker='o', c='k')
point2 = ax2.plot([], [], marker='o', c='k')
point3 = ax3.plot([], [], marker='o', c='k')
point = [point1, point2, point3]

xVals = np.arange(-10, 0.05, 0.01)
xVals2 = np.arange(-10, 0.47, 0.01)
nDelim = 0.61 / (xVals - 0.05) + 1.3
nDelim2 = 0.61 / (xVals2 - 0.47) + 1.19
ax1.set_title("BPT Diagram Using NII")
ax1.set_xlabel("[NII]/H-alpha")
ax1.set_ylabel("[OIII]/H-beta")
ax1.text(-0.6, -0.4, "HII", color='b')
ax1.text(-0.3, 0, "Composite", color='b')
ax1.text(0, 0.6, "AGN", color='b')
ax1.plot(xVals, nDelim, c='k', dashes=[6,2])
ax1.plot(xVals2, nDelim2, c='k')

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

plt.tight_layout()

#setting limits
for ax in [ax1, ax2, ax3]:
    ax.set_ylim(-3, 3)
    ax.set_xlim(-3, 3)
    ax.grid()

#initializing data arrays
xNdata, yNdata, xSdata, ySdata, xOdata, yOdata = [], [], [], [], [], []

def update(data):
    print("in update!")
    xN, yN, xS, yS, xO, yO = data #the named tuple in your version
    xNdata.append(xN)
    yNdata.append(yN)
    xSdata.append(xS)
    ySdata.append(yS)
    xOdata.append(xO)
    yOdata.append(yO)
    
    plt.draw()
    #return point
    

for i, line in enumerate(total):
    print("iteration " + str(i))
    #print(line)
    if i==0:
        continue
    else:
        data = [], [], [], [], [], []
        temp1 = line.strip().split()
        cBPT = namedtuple('cBPT', 'hAlpha NII hBeta OIII SII OI noise hAisClear NIIisClear hBisClear OIIIisClear SIIisClear OIisClear tooNoisy')
        BPT = cBPT(float(temp1[1]), float(temp1[2]), float(temp1[3]), float(temp1[4]), float(temp1[5]), float(temp1[6]), float(temp1[7]), float(temp1[8]), float(temp1[9]), float(temp1[10]), float(temp1[11]), float(temp1[12]), float(temp1[13]), float(temp1[14]))
        if BPT.tooNoisy==True:
            print("tooNoisy, skipping: " + str(i))
            continue
        elif BPT.tooNoisy==False:
            #xN, yN, xS, yS, xO, yO = data
            data = [np.log(BPT.NII/BPT.hAlpha), np.log(BPT.OIII/BPT.hBeta), np.log(BPT.SII/BPT.hAlpha), np.log(BPT.OIII/BPT.hBeta), np.log(BPT.OI/BPT.hAlpha), np.log(BPT.OIII/BPT.hBeta)]
            #update(data)
            print(data)
            update(data)
            ax1.plot(xNdata, yNdata)
            ax2.plot(xSdata, ySdata)
            ax3.plot(xOdata, yOdata)
            print("after ax.plot called!")
            c.show()
            #ani = animation.FuncAnimation(c, update, data, blit=True, interval=2, repeat=False)
            
    c.show()
