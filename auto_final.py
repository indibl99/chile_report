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

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def double_gaussian(x, a1, x01, sigma1, a2, x02, sigma2):
    return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2)

def triple_gaussian(x, a1, x01, sigma1, a2, x02, sigma2, a3, x03, sigma3):
    return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2) + func(x, a3, x03, sigma3)

def errfunc(guess, xData, yData):
    return (func(xData, guess[0], guess[1], guess[2]) - yData)

def doubleErrFunc(guess, xData, yData):
    yFit = double_gaussian(xData, guess[0], guess[1], guess[2], guess[3], guess[4], guess[5])
    err2 = yData - yFit
    return err2

def tripleErrFunc(guess, xData, yData):
    yFit = triple_gaussian(xData, guess[0], guess[1], guess[2], guess[3], guess[4], guess[5], guess[6], guess[7], guess[8])
    err3 = yData - yFit
    return err3

def parse_datafile(base, filename):
    #parsing STARLIGHT output
    datapathIn = '/Users/plira/india_temp/variable_spectra/auto/starlight_models/' + base + '_output.txt'
    datapathOut = '/Users/plira/india_temp/variable_spectra/auto/output_spectra/' + base + '/'
    dataOut = datapathOut + base + '_parsed.txt'

    if os.path.exists(datapathOut):
        shutil.rmtree(datapathOut)
    os.makedirs(datapathOut)
    starModel = open(dataOut, 'w')

    with open(datapathIn) as star:
        for i, line in enumerate(star):
            if i == 25:
                temp = line.split(' ')
                fobs = float(temp[0])
                print("fobs: " + str(fobs))
            if i == 31:
                temp2 = line.strip().split()
                print(temp2)
                SN = float(temp2[0])
                print("SN: " + str(SN))
            elif i >= 216:
                starModel.writelines(line)    
    starData = [fobs, SN]
    return starData

def loadData(filename, fobs, outLoc):
    loc = filename
    spec = np.genfromtxt(filename)
    
    wave = spec[:,0]
    inSpec = fobs*spec[:,1]
    outSpec = fobs*spec[:,2]
    subSpec = inSpec - outSpec

   #plot observed spectrum, STARLIGHT model, subtracted model
    m = plt.figure(1, figsize=(15,5))
    plt.clf()
    plt.xlim(3700,7700)

    peaks = indexes(inSpec, thres=0.3, min_dist=4.0)
    peaksY = inSpec[peaks]
    peaksY = peaksY[:2000]
    yMax = np.amax(peaksY) + 30 #greater than 5500
    print('yMax: ' + str(yMax))
    
    plt.ylim(-20, yMax)
    plt.grid()
    plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
    plt.xlabel("wavelength (angstroms)")
    plt.ylabel("flux")
    title = "Spectrum and Subtracted STARLIGHT Model of \n" + filename
    plt.title(title)
    plt.plot(wave, inSpec, c='k', linewidth=0.4)
    plt.plot(wave, outSpec, c='r', linewidth=0.4)
    plt.plot(wave, subSpec, c='b', linewidth=0.4)

    #plt.autoscale(enable=True, axis='y', tight=True)

    plt.vlines(4861, -20, yMax, colors='r', linewidth=3.0, alpha=0.5)
    plt.text(4861, 40, "H-beta (4861)")
    plt.vlines(5007, -20, yMax, colors='r', linewidth=3.0, alpha=0.5)
    plt.text(5007, 30, "OIII (5007)")

    plt.vlines(6300, -20, yMax, colors='r', linewidth=3.0, alpha=0.5)
    plt.text(6300, 30, "OI (6300)")

    plt.vlines(6563, -20, yMax, colors='r', linewidth=2.0, alpha=0.5)
    plt.text(6563, 50, "H-alpha (6563)")
    plt.vlines(6583, -20, yMax, colors='r', linewidth=3.0, alpha=0.5)
    plt.text(6583, 42, "NII (6583)")

    plt.vlines(6716, -20, yMax, colors='r', linewidth=3.0, alpha=0.5)
    plt.text(6716, 35, "SII (6716)")
    plt.vlines(6736, -20, yMax, colors='r', linewidth=3.0, alpha=0.5)
    plt.text(6736, 28, "SII (6736)")

    m.savefig(outLoc + 'starlight_model.png',bbox_inches='tight')

    OutData = namedtuple('outData', 'wave inSpec outSpec subSpec')
    data = OutData(wave, inSpec, outSpec, subSpec)
    return data

def gaussianFits(fileName, fobs, outLoc, SN):

    data = loadData(fileName, fobs, outLoc)
    x = data.wave #could also write as hAregion[:,0]
    y = data.subSpec
    width = 3.0 #know this from experience 
    amphB = 0
    ampOIII = 0
    ampOI = 0
    amphA = 0
    ampNII = 0
    ampNII0 = 0
    ampSIIa = 0
    ampSIIb = 0
    xGuessSIIb = 6735
    
    #may need to make noise smarter soon - this is a test!!
    #NEED TO MAKE NOISE MUCH SMARTER
    noiseX = x[1450:1460]
    noiseY = y[1450:1460]
    
    noise = np.abs(3 * integrate.simps(noiseY, noiseX))
    
    #trying to fix the problem with the gaussian fits using peak utils indexes
    #automating guessing heights
    peaks = indexes(y, thres=0.3, min_dist=4.0) #these values are somewhat arbitrary right now
    numPeaks1 = 0
    peak1 = [0]
    peakIndices1 = [0]
    numPeaks2 = 0
    peak2 = [0]
    peakIndices2 = [0]
    peakX = x[peaks]
    peakY = y[peaks]
    #buffer = 0
    newRange = len(peakY)
    negRange = len(peakY)
    for height in peakY:
        if height > newRange:
            break
        if peakY[height] < 0:
            np.delete(peakY, height)
            np.delete(peakX, height)
            newRange = newRange - 1
    print(peakX)

    hBisClear = False
    OIIIisClear = False
    OIisClear = False
    hAisClear = False
    NIIisClear = False
    SIIisClear = False
    tooNoisy = False

    Halpha = 0
    NII = 0
    Hbeta = 0
    OIII = 0
    SII = 0
    OI = 0
    doublePeak = 0
    hBnum = 0
    hBsens = []
    
    for i in range(len(peakX)):
        """if peakX[i] > 6550 and peakX[i] < 6780: #H-alpha, NII, and SII
            numPeaks2 = numPeaks2 + 1
            peak2.append(peakX[i])
            peakIndices2.append(i)"""
        if peakX[i] > 6557 and peakX[i] < 6567:
            numPeaks2 = numPeaks2 + 1
            peak2.append(peakX[i])
            peakIndices2.append(i)
            hAisClear = True
            amphA = peakY[i]
        if peakX[i] > 6577 and peakX[i] < 6587:
            numPeaks2 = numPeaks2 + 1
            peak2.append(peakX[i])
            peakIndices2.append(i)
            NIIisClear = True
            ampNII = peakY[i]
        if peakX[i] > 6710 and peakX[i] < 6740:
            doublePeak = doublePeak + 1
            numPeaks2 = numPeaks2 + 1
            peak2.append(peakX[i])
            peakIndices2.append(i)
            SIIisClear = True
            """if doublePeak==1:
                ampSIIa = peakY[i]
            elif doublePeak==2:
                ampSIIb = peakY[i]
                if np.abs(peakX[i] - xGuessSIIb) >= 2:
                    xGuessSIIb = peakX[i]""" 
        if peakX[i] > 6713 and peakX[i] < 6722:
            ampSIIa = peakY[i]
            print('SIIa (' + str(peakX[i]) + ', ' + str(peakY[i]) + ')')
        if peakX[i] > 6727 and peakX[i] < 6733:
            ampSIIb = peakY[i]
            print('SIIb (' + str(peakX[i]) + ', ' + str(peakY[i]) + ')')
        if peakX[i] > 6577 and peakX[i] < 6587:
            numPeaks2 = numPeaks2 + 1
            peak2.append(peakX[i])
            peakIndices2.append(i)
            NIIisClear = True
            ampNII = peakY[i]
        if peakX[i] > 4850 and peakX[i] < 4870: #H-beta
            numPeaks1 = numPeaks1 + 1
            hBnum = hBnum + 1
            peak1.append(peakX[i])
            peakIndices1.append(i)
            hBisClear = True
            if hBnum >= 1:
                hBsens.append(peakY[i])
                amphB = np.amax(hBsens)
            #amphB = peakY[i]
        if peakX[i] > 5000 and peakX[i] < 5015: #OIII
            numPeaks1 = numPeaks1 + 1
            peak1.append(peakX[i])
            peakIndices1.append(i)
            OIIIisClear = True
            ampOIII = peakY[i]
        if peakX[i] > 6295 and peakX[i] < 6305: #OI
            numPeaks1 = numPeaks1 + 1
            peak1.append(peakX[i])
            peakIndices1.append(i)
            OIisClear = True
            ampOI = peakY[i]
    peak2.remove(0)
    peak1.remove(0)
    peakIndices2.remove(0)
    peakIndices1.remove(0)
    print("numPeaks2 = " + str(numPeaks2))
    print("numPeaks1 = " + str(numPeaks1))

    print(peak2)
    print(peak1)

    #note: OI is usually so small, indexes probably won't pick it up
    #tooNoisy = False
  
    """if numPeaks2==5:
        ampNII0 = peakY[peakIndices2[0]]
        amphA = peakY[peakIndices2[1]]
        ampNII = peakY[peakIndices2[2]]
        ampSIIa = peakY[peakIndices2[3]]
        ampSIIb = peakY[peakIndices2[4]]
    elif numPeaks2==4:
        ampNII0 = noise #placeholder because you need to do a triple gaussian fit 
        amphA = peakY[peakIndices2[0]]
        ampNII = peakY[peakIndices2[1]]
        ampSIIa = peakY[peakIndices2[2]]
        ampSIIb = peakY[peakIndices2[3]]
    if numPeaks2 < 3 or numPeaks2 > 9:
        tooNoisy = True
        print("too few peaks detected, check this spectra by hand")
    elif numPeaks2 > 25:
        tooNoisy = True
        print("too many peaks detected, check this spectra by hand")"""

    """if numPeaks1==2:
        amphB = peakY[peakIndices1[0]]
        ampOIII = peakY[peakIndices1[1]]
        ampOI = noise #this is just a place holder for now, means that .indexes couldn't find a peak at OI
    elif numPeaks1==3: #note! need to take care of uncertainties 
        amphB = peakY[peakIndices1[0]]
        ampOIII = peakY[peakIndices1[1]]
        ampOI = peakY[peakIndices1[2]]"""
    if numPeaks1 < 1: #NOT LESS THAN TWO (gets rid of cases where H-b is visible but not O3, which happens a lot)
        tooNoisy = True
        print("too few peaks detected, check this spectra by hand")
    if hAisClear==False: tooNoisy = True
    #if SN < 13: tooNoisy = True
    #if hBisClear==False and OIIIisClear==False: tooNoisy = True
    #print(tooNoisy)
    #print("LOOK HERE")
    tooNoisy = False
    if tooNoisy==False: 
        guesshB = [amphB, 4865, width]
        guessOIII = [ampOIII, 5007, width]
        guessOI = [ampOI, 6300, width]
        guesshA = [amphA, 6563, width]
        guessNII = [ampNII, 6583, width]
        guessSIIa = [ampSIIa, 6716, width] #6716
        guessSIIb = [ampSIIb, 6731, width] #6731

        #using the least square function to optimize the paramters for the gaussian fit(the params for the func() function)
        
        #singular-peak gaussian fits
        optimhB, flag = sp.leastsq(errfunc, guesshB, args=(x, y))
        optimOIII, flag = sp.leastsq(errfunc, guessOIII, args=(x, y))
        optimOI, flag = sp.leastsq(errfunc, guessOI, args=(x, y))

        #multi-peak gaussian fits
        guessSII = [guessSIIa[0], guessSIIa[1], guessSIIa[2], guessSIIb[0], guessSIIb[1], guessSIIb[2]]
        optimSII, flag = sp.leastsq(doubleErrFunc, guessSII, args=(x, y))
        ySII = double_gaussian(x, optimSII[0], optimSII[1], optimSII[2], optimSII[3], optimSII[4], optimSII[5])

        guesshANII = [ampNII0, 6548, width, guesshA[0], guesshA[1], guesshA[2], guessNII[0], guessNII[1], guessNII[2]]
        optimhANII, flag = sp.leastsq(tripleErrFunc, guesshANII, args=(x, y))
        yhANII = triple_gaussian(x, optimhANII[0], optimhANII[1], optimhANII[2], optimhANII[3], optimhANII[4], optimhANII[5], optimhANII[6], optimhANII[7], optimhANII[8])

        #now, calculating y-values and extracting each individual peak from multi-peak fit
        yNII0 = func(x, optimhANII[0], optimhANII[1], optimhANII[2])
        yhA = func(x, optimhANII[3], optimhANII[4], optimhANII[5])
        yNII = func(x, optimhANII[6], optimhANII[7], optimhANII[8])

        ySIIa = func(x, optimSII[0], optimSII[1], optimSII[2])
        ySIIb = func(x, optimSII[3], optimSII[4], optimSII[5])
        
        #calculating y values for a gaussian fit using new input paramters optimized above
        yhB = func(x, optimhB[0], optimhB[1], optimhB[2])
        yOIII = func(x, optimOIII[0], optimOIII[1], optimOIII[2])
        yOI = func(x, optimOI[0], optimOI[1], optimOI[2])

        
        Hbeta = integrate.simps(yhB, x)
        OIII = integrate.simps(yOIII, x)
        OI = integrate.simps(yOI, x)
        Halpha = integrate.simps(yhA, x)
        NII = integrate.simps(yNII, x)
        #SII = integrate.simps(ySIIa, x) + integrate.simps(ySIIb, x) #not sure which of these two it is (are they different?)
        SII = integrate.simps(ySII, x)

        if Hbeta < 0: hBisClear = False
        if OIII < 0: OIIIisClear = False
        if OI < 0: OIisClear = False
        if Halpha < 0: hAisClear = False
        if NII < 0: NIIisClear = False
        if SII < 0: SIIisClear = False

        f = plt.figure(2, figsize=(16,5))
        plt.clf()
        plt.xlim(4700, 7000)
        plt.grid()
        plt.tight_layout()
        plt.title('Gaussian Fit')

        plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
        plt.plot(x, y, c='k', linewidth=0.5)
        
        plt.plot(x, yhB, c='c', linewidth=0.75)
        plt.plot(x, yOIII, c='y', linewidth=0.75)
        plt.plot(x, yOI, c='m', linewidth=0.75)
        plt.plot(x, yhANII, c='b', lw=0.75)
        plt.plot(x, yhA, c='g', lw =0.75)
        plt.plot(x, yNII0, c='y', lw =0.75)
        plt.plot(x, yNII, c='g', lw =0.75)
        plt.plot(x, ySII, c='b', lw=0.75)
        plt.plot(x, ySIIa, c='r', lw=0.75)
        plt.plot(x, ySIIb, c='r', lw=0.75)
        
        f.savefig(outLoc + 'gaussian_fit.png', bbox_inches='tight')
        
    elif tooNoisy==True:
        print("spectrum is too noisy to analyze automatically, please check by hand")

    elineFluxes = namedtuple('EmLines', 'hAlpha NII hBeta OIII SII OI noise hAisClear NIIisClear hBisClear OIIIisClear SIIisClear OIisClear tooNoisy')
    emLines = elineFluxes(Halpha, NII, Hbeta, OIII, SII, OI, noise, hAisClear, NIIisClear, hBisClear, OIIIisClear, SIIisClear, OIisClear, tooNoisy)
    return emLines #potential problem here! maybe it has to return something no matter what!
        

def BPT(emLines, outLoc):

    bptNII=  ''
    bptSII = ''
    bptOI = ''

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

    Halpha = emLines.hAlpha
    NII = emLines.NII
    Hbeta = emLines.hBeta
    OIII = emLines.OIII
    SII = emLines.SII
    OI = emLines.OI
    noise = emLines.noise

    hAisClear = emLines.hAisClear
    hAClear = emLines.hAisClear
    NIIisClear = emLines.NIIisClear
    NIIClear = emLines.NIIisClear
    hBisClear = emLines.hBisClear
    hBClear = emLines.hBisClear
    OIIIisClear = emLines.OIIIisClear
    OIIIClear = emLines.OIIIisClear
    SIIisClear = emLines.SIIisClear
    SIIClear = emLines.SIIisClear
    OIisClear = emLines.OIisClear
    OIClear = emLines.OIisClear

    flagY = 0 #-1 is down, 1 is up
    flagX = 0 #-1 is pointing to the left, 1 is pointing to the right
    allClear = True
    allClearNII = True
    allClearSII = True
    allClearOI = True
    
    
    if hAisClear==False and Halpha < noise: #right now there is a problem if both H-alpha and NII/SII/OI are both uncertain, same for OIII and HB
        Halpha = noise   #but if H-alpha is uncertain, the spectrum is too noisy anyways 
        flagX = 1
        print("hA is uncertain")
        hAClear = False
        allClear = False
    if NIIisClear==False and NII < noise:
        NII = noise
        flagX = -1
        print("NII is uncertain")
        NIIClear = False
        allClearNII = False
    if hBisClear==False and Hbeta < noise:
        Hbeta = noise
        flagY = 1
        print("hB is uncertain")
        hBClear = False
        allClear = False
    if OIIIisClear==False and OIII < noise:
        OIII = noise
        flagY = -1
        print("OIII is uncertain")
        OIIIClear = False
        allClear = False
    if SIIisClear==False and SII < noise:
        SII = noise
        flagX = -1
        print("SII is uncertain")
        SIIClear = False
        allClearSII = False
    if OIisClear==False or OI < noise: #this one is different because OI is very small and often unclear
        OI = noise
        flagX = -1
        print("OI is uncertain")
        allClearOI = False
        OIClear = False
        
    print("H-beta: ", Hbeta)
    print("OIII: ", OIII)
    print("OI: ", OI)
    print("H-alpha: ", Halpha)
    print("NII: ", NII)
    print("SII: ", SII)
    print("noise: ", noise)

    plt.clf()
    c, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 10))#3 subplots
    for ax in [ax1, ax2, ax3]:
        ax.set_ylim(-3, 3)
        ax.set_xlim(-5, 3)
        ax.grid()

    #format and plot the delimeter equations
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

    #BPT NII
    
    #color-code dot based on region
    colorVarNII = 'k'
    delim1NII = 0.61 / (np.log(NII/Halpha) - 0.47) + 1.19
    delim2NII = 0.61 / (np.log(NII/Halpha) - 0.05) + 1.3
    if np.log(OIII/Hbeta) >= delim1NII or np.log(NII/Halpha) >= 0.47:
        colorVarNII = '#66ff33'
        bptNII = 'AGN'
    elif np.log(OIII/Hbeta) >= delim2NII and np.log(OIII/Hbeta) < delim1NII:
        colorVarNII = '#ffff00'
        bptNII = 'Composite'
    elif np.log(OIII/Hbeta) < delim2NII:
        colorVarNII = 'r'
        bptNII = 'HII'

    ax1.plot(np.log(NII/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarNII, markeredgecolor=colorVarNII)
    if allClear==False or allClearNII==False: ax1.quiver(np.log(NII/Halpha), np.log(OIII/Hbeta), flagX, flagY,color=colorVarNII) #adds an arrow if uncertain


    #SII BPT

    #format and plot delimiter equations
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
        
    ax2.plot(np.log(SII/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarSII, markeredgecolor=colorVarSII)
    if allClear==False or allClearSII==False: ax2.quiver(np.log(SII/Halpha), np.log(OIII/Hbeta), flagX, flagY, color=colorVarSII)

    #OI BPT 
    #k = plt.figure(5)
    
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
    elif np.log(OI/Halpha) >= -0.59 and np.log(OIII/Hbeta) < delim2OI:
        colorfVarOI = '#ffff00'
        bptOI = 'LINER'
    elif np.log(OIII/Hbeta) < delim1OI: #check this
        colorVarOI = 'r'
        bptOI = 'HII'

    ax3.plot(np.log(OI/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarOI, markeredgecolor=colorVarOI)
    if allClear==False or allClearOI==False: ax3.quiver(np.log(OI/Halpha), np.log(OIII/Hbeta), flagX, flagY, color=colorVarOI)
    plt.tight_layout()
    c.savefig(outLoc + 'BPT Diagrams', bbox_inches='tight')

    bptR = [bptNII, bptSII, bptOI, hAClear, NIIClear, hBClear, OIIIClear, SIIClear, OIClear]
    
    return bptR

def funcToDeleteItems(fullPathToDir):
    shutil.rmtree(fullPathToDir)
                
def run_code():
    #for a first run, fill the directories with like 5 spectra of different kinds before you do them all at once, and then when you do, make sure you have a couple hours
    #and save as you go, so that you can quit in the middle of the loop if need be 
    
    outpath = '/Users/plira/india_temp/variable_spectra/auto/output_spectra_2/'
    starlight = '/Users/plira/software/STARLIGHTv04/'
    oldFile = 'old'
    totalFile = '/Users/plira/india_temp/variable_spectra/auto/cumulative_output_2.txt'
    aniBPT = '/Users/plira/india_temp/variable_spectra/auto/cumulative_BPT_2.txt'
    #totalBPT = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_BPT.txt', 'w')
    inputDir = '/Users/plira/india_temp/variable_spectra/auto/input_spectra_2/'
    classDir = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_2/'

    cumulative = open(totalFile, 'w')
    cumulative.write("{: <20} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <15} {: <15} {: <15}".format('Filename','hAlpha','NII','hBeta','OIII','SII','OI','noise','hAisClear','NIIisClear','hBisClear','OIIIisClear','SIIisClear','OIisClear','tooNoisy', 'BPT NII', 'BPT SII', 'BPT OI'))
    
    if os.path.exists(classDir + 'likely_AGN/'): funcToDeleteItems(classDir + 'likely_AGN/')
    if os.path.exists(classDir + 'NII/'): funcToDeleteItems(classDir + 'NII/')
    if os.path.exists(classDir + 'OI/'): funcToDeleteItems(classDir + 'OI/')
    if os.path.exists(classDir + 'SII/'): funcToDeleteItems(classDir + 'SII/')
    if os.path.exists(classDir + 'tooNoisy/'): funcToDeleteItems(classDir + 'tooNoisy/')
    
    for filename in os.listdir('/Users/plira/india_temp/variable_spectra/auto/input_spectra_2'): 
        if filename=='.DS_Store':
            continue
        else:
            totalBPT = open(aniBPT, 'a')
            print('old filename: ' + filename)
            os.rename(inputDir + filename, inputDir+filename[-24:])
            filename2 = filename[-24:]
            print('new filename: ' + filename2)
            base = filename2[:-4]
            print("base: " + base)
            #temporarily commenting out starlight
            """grid = '/Users/plira/software/STARLIGHTv04/grid_auto.in'
            for line in fileinput.input([grid], inplace=True):
                line = line.replace(oldFile, base)
                sys.stdout.write(line)

            os.chdir(starlight)
            call(["pwd"])
            os.system("./Starlight_v04_Mac.exe < grid_auto.in") #:)"""
            #the way this is set up, if you just comment out the STARLIGHT section, it should run without running starlight? check that to be sure
            starData = parse_datafile(base, filename2)
            fobs = starData[0]
            SN = starData[1]
                
            locI = outpath + base + '/'
            locO = '/Users/plira/india_temp/variable_spectra/auto/input_spectra_2/' + filename2
            fileOut = open(locI + 'analysis output', 'w')
            print("locI: " + locI)
            shutil.copy2(locO, locI)
            shutil.copy2('/Users/plira/india_temp/variable_spectra/auto/starlight_models_2/'+base+'_output.txt', locI)
            loc = outpath + base + '/' + base + '_parsed.txt' #this is the PARSED model that gets passed
            emLines = gaussianFits(loc, fobs, locI, SN);
            output = [base, emLines.hAlpha, emLines.NII, emLines.hBeta, emLines.OIII, emLines.SII, emLines.OI, emLines.noise, emLines.hAisClear, emLines.NIIisClear, emLines.hBisClear, emLines.OIIIisClear, emLines.SIIisClear, emLines.OIisClear, emLines.tooNoisy]
            cumulative.write("\n{: <20} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10}".format(output[0], output[1], output[2], output[3], output[4], output[5], output[6], output[7], output[8], output[9], output[10], output[11], output[12], output[13], output[14]))

            fileOut.write('Filename: '+ output[0])
            fileOut.write('\n ')
            fileOut.write("\n{: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15}".format('hAlpha','NII','hBeta','OIII','SII', 'OI', 'noise'))
            fileOut.write("\n{: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15}".format(output[1], output[2], output[3], output[4], output[5], output[6], output[7]))
            fileOut.write('\n ')
            fileOut.write("\n{: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10}".format('hAisClear','NIIisClear','hBisClear','OIIIisClear','SIIisClear','OIisClear','tooNoisy'))
            fileOut.write("\n{: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10}".format(output[8], output[9], output[10], output[11], output[12], output[13], output[14]))        

            if emLines.tooNoisy==False:
                bptO = BPT(emLines, locI)
                fileOut.write('\n ')
                fileOut.write("\n{: <15} {: <15} {: <15}".format('BPT NII', 'BPT SII', 'BPT OI'))
                fileOut.write("\n{: <15} {: <15} {: <15}".format(bptO[0], bptO[1], bptO[2]))
                cumulative.write("{: <15} {: <15} {: <15}".format(bptO[0], bptO[1], bptO[2]))
                totalBPT.write('\n'+ str(output[0]) + ',' + str(output[1]) + ',' + str(output[2]) + ',' + str(output[3]) + ',' + str(output[4]) + ',' + str(output[5]) + ',' + str(output[6]) + ',' + str(output[7]) + ',' + str(output[8]) + ',' + str(output[9]) + ',' + str(output[10]) + ',' + str(output[11]) + ',' + str(output[12]) + ',' + str(output[13]) + ',' + str(output[14]))
                fileOut.write('\n ')
                fileOut.write("\n{: <10} {: <10} {: <10} {: <10} {: <10} {: <10}".format('hAisClear2','NIIisClear2','hBisClear2','OIIIisClear2','SIIisClear2','OIisClear2'))
                fileOut.write("\n{: <13} {: <13} {: <13} {: <13} {: <13} {: <13}".format(bptO[3], bptO[4], bptO[5], bptO[6], bptO[7], bptO[8]))        
                fileOut.close()
                             
                copy_tree(outpath + base, classDir + 'NII/' + bptO[0] + '/' + base + '/')
                copy_tree(outpath + base, classDir + 'SII/' + bptO[1] + '/' + base + '/')
                copy_tree(outpath + base, classDir + 'OI/' + bptO[2] + '/' + base + '/')
                if(bptO[0]=='AGN' and bptO[1]=='Seyfert' and bptO[2]=='Seyfert'):
                    copy_tree(outpath + base, classDir + 'likely_AGN/' + '/' + base + '/')
                
            elif emLines.tooNoisy==True:
                notWorthIt = open(locI + 'ErrorMessage_noisy_spectrum', 'w')
                notWorthIt.write('spectrum was too noisy, please inspect file ' + base + '.txt by hand')
                notWorthIt.close()
                copy_tree(outpath + base, classDir + 'tooNoisy/' + base + '/')
            totalBPT.close()

        #fileOut.close()
        oldFile = base
        print("oldFile: " + oldFile)
        
    cumulative.close()

if __name__ == '__main__':
    run_code()

