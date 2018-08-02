import numpy as np
import scipy.optimize as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys
import shutil
from shutil import copyfile
from collections import namedtuple
from pathlib import Path
import os

totalFile = '/Users/plira/india_temp/variable_spectra/auto/cumulative_output.txt'
with open(totalFile) as total:
    #c = plt.figure(7)
    for i, line in enumerate(total):
        if i==0:
            continue
        else:
            temp1 = line.strip().split()
            print(temp1)
            cBPT = namedtuple('cBPT', 'hAlpha NII hBeta OIII SII OI noise hAisClear NIIisClear hBisClear OIIIisClear SIIisClear OIisClear tooNoisy')
            BPTparams = cBPT(float(temp1[1]), float(temp1[2]), float(temp1[3]), float(temp1[4]), float(temp1[5]), float(temp1[6]), float(temp1[7]), float(temp1[8]), float(temp1[9]), float(temp1[10]), float(temp1[11]), float(temp1[12]), float(temp1[13]), float(temp1[14]))
            if BPTparams.tooNoisy==True:
                print("tooNoisy, skipping: " + str(i))
            elif cBPT.tooNoisy==False:
                BPT(BPTparams)
#key=1 is NII, =2 is SII, =3 is OI
def BPT(hAlpha, NII, hBeta, OIII, SII, OI, noise, hAisClear, NIIisClear, hBisClear, OIIIisClear, SIIisClear, OIisClear, tooNoisy): 
    for root, dirs, filenames in os.walk('/users/indiabhalla-ladd/desktop/variable_spectra/output_files/'):
        for f in filenames:
            if f.endswith(".txt"):
                #print(f)
                path = '/users/indiabhalla-ladd/desktop/variable_spectra/output_files/' + f
                #print(path)
                output = np.genfromtxt(path)

                Halpha = output[0]
                NII = output[1]
                Hbeta = output[2]
                OIII = output[3]
                SII = output[4]
                OI = output[5]
                noise = output[6]

                if key==1:
                    plt.plot(np.log(NII/Halpha), np.log(OIII/Hbeta), marker='o', c='k')
                elif key==2:
                    plt.plot(np.log(SII/Halpha), np.log(OIII/Hbeta), marker='o', c='k')
                elif key==3:
                    plt.plot(np.log(OI/Halpha), np.log(OIII/Hbeta), marker='o', c='k')
                    
        if key==1:
            xVals = np.arange(-10, 0.05, 0.01)
            xVals2 = np.arange(-10, 0.47, 0.01)
            nDelim = 0.61 / (xVals - 0.05) + 1.3
            nDelim2 = 0.61 / (xVals2 - 0.47) + 1.19
            plt.title("BPT Diagram Using NII")
            plt.xlabel("[NII]/H-alpha")
            plt.ylabel("[OIII]/H-beta")
            plt.text(-0.6, -0.4, "HII", color='b')
            plt.text(-0.3, 0, "Composite", color='b')
            plt.text(0, 0.6, "AGN", color='b')
            plt.plot(xVals, nDelim, c='k', dashes=[6,2])
            plt.plot(xVals2, nDelim2, c='k')

        elif key==2:
            xS = np.arange(-6, 0.32, 0.03)
            xS2 = np.arange(-0.315, 4, 0.01)
            sDelim = 0.72 / (xS - 0.32) + 1.3
            sDelim2 = 1.89*xS2 + 0.76 
            plt.title("BPT Diagram Using SII")
            plt.xlabel("[SII]/H-alpha")
            plt.ylabel("[OIII]/H-beta")
            plt.text(-0.4, -0.4, "HII", color='b')
            plt.text(-0.5, 0.6, "Seyfert", color='b')
            plt.text(0.2, -0.2, "LINER", color='b')
            plt.plot(xS, sDelim, c='k', dashes=[6,2])
            plt.plot(xS2, sDelim2, c='k')
            
        elif key==3:
            xO = np.arange(-6, -0.59, 0.01)
            xO2 = np.arange(-1.13, 4, 0.01)
            oDelim = 0.73 / (xO + 0.59) + 1.33
            oDelim2 = 1.18*xO2 + 1.3
            plt.title("BPT Diagram Using OI")
            plt.xlabel("[OI]/H-alpha")
            plt.ylabel("[OIII]/H-beta")
            plt.text(-1.6, -0.4, "HII", color='b')
            plt.text(-1, 0.6, "Seyfert", color='b')
            plt.text(-0.5, -0.2, "LINER", color='b')
            plt.plot(xO, oDelim, c='k', dashes=[6,2])
            plt.plot(xO2, oDelim2, c='k')

    plt.show()

plotBPT(1)
plotBPT(2)
plotBPT(3)

