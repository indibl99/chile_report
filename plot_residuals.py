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

def parse_datafile(base):
    #parsing STARLIGHT output
    datapathIn = '/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/load_individual/' + base + '_in/' + base + '_output.txt'
    datapathOut = '/Users/plira/india_temp/variable_spectra/auto/output_spectra/spectral_classifications_final/load_individual/' + base + '/'

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

def loadData(filenames):
    stack = np.zeros(4201)
    waveStack = []
    for spectrum in filenames:
        spec = np.genfromtxt(spectrum)
        
        wave = spec[:,0] #no fobs
        inSpec = spec[:,1]
        outSpec = spec[:,2]
        subSpec = inSpec - outSpec

        stack = stack + subSpec #needs to be smarter than this, and use gaussian to align the centers
        waveStack = wave

        plt.xlim(3700,7700)
        plt.grid()
        plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
        plt.xlabel("wavelength (angstroms)")
        plt.ylabel("flux")
        
        #plt.plot(wave, subSpec, c='r', linewidth=0.4)

    plt.plot(waveStack, stack, c = 'k', linewidth=0.4)
    plt.show()
    raw_input()

    data = [waveStack, stack]

    OutData = namedtuple('data', 'wave spec')
    data = OutData(waveStack, stack)
    return data
    


def gaussianFits(fileName):

    data = loadData(fileName)
    x = data.wave
    y = data.spec
    width = 3.0
    amphB = 0
    ampOIII = 0
    ampOI = 0
    amphA = 0
    ampNII = 0
    ampNII0 = 0
    ampSIIa = 0
    ampSIIb = 0
    
    peaks = indexes(y, thres=0.3, min_dist=4.0) #these values are somewhat arbitrary right now
    numPeaks1 = 0
    peak1 = [0]
    peakIndices1 = [0]
    numPeaks2 = 0
    peak2 = [0]
    peakIndices2 = [0]
    peakX = x[peaks]
    peakY = y[peaks]

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
    tooNoisy = False
    
    if tooNoisy==False: 
        guesshB = [amphB, 4865, width]
        guessOIII = [ampOIII, 5007, width]
        guessOI = [ampOI, 6300, width]
        guesshA = [amphA, 6563, width]
        guessNII = [ampNII, 6583, width]
        guessSIIa = [ampSIIa, 6716, width] #6716
        guessSIIb = [ampSIIb, 6731, width] #6731

        """guesshB = [48, 4865, 3]
        guessOIII = [15, 5007, 3]
        guessOI = [5, 6300, 3]
        guesshA = [10, 6563, 3]
        guessNII = [5, 6583, 3]
        ampNII0 = 4
        width = 3"""
        #guessSIIa = [3.1, 6716, 3] #6716
        #guessSIIb = [2.1, 6735, 1]
        #4.4 and 3.6

        #using the least square function to optimize the paramters for the gaussian fit(the params for the func() function)
        
        #singular-peak gaussian fits
        optimhB, flag = sp.leastsq(errfunc, guesshB, args=(x, y))
        optimOIII, flag = sp.leastsq(errfunc, guessOIII, args=(x, y))
        optimOI, flag = sp.leastsq(errfunc, guessOI, args=(x, y))

        #multi-peak gaussian fits
        guessSII = [guessSIIa[0], guessSIIa[1], guessSIIa[2], guessSIIb[0], guessSIIb[1], guessSIIb[2]]
        #guessSII = [35, 6716, 3, 30, 6736, 3]
        optimSII, flag = sp.leastsq(doubleErrFunc, guessSII, args=(x, y))
        ySII = double_gaussian(x, optimSII[0], optimSII[1], optimSII[2], optimSII[3], optimSII[4], optimSII[5])

        #writing H-alpha and NII as a double gaussian rn, un comment when done
        guesshANII = [ampNII0, 6548, width, guesshA[0], guesshA[1], guesshA[2], guessNII[0], guessNII[1], guessNII[2]]
        optimhANII, flag = sp.leastsq(tripleErrFunc, guesshANII, args=(x, y))
        yhANII = triple_gaussian(x, optimhANII[0], optimhANII[1], optimhANII[2], optimhANII[3], optimhANII[4], optimhANII[5], optimhANII[6], optimhANII[7], optimhANII[8])
        
        
        """guesshANII2 = [guesshA[0], guesshA[1], guesshA[2], guessNII[0], guessNII[1], guessNII[2]]
        guesshANII2 = [38.9, 6565, 7, 13.7, 6586, 3]
        optimhANII, flag = sp.leastsq(doubleErrFunc, guesshANII2, args=(x, y))
        yhANII2 = double_gaussian(x, optimhANII[0], optimhANII[1], optimhANII[2], optimhANII[3], optimhANII[4], optimhANII[5])"""
        
        optimhA, flag = sp.leastsq(errfunc, guesshA, args=(x,y))
        optimNII, flag = sp.leastsq(errfunc, guessNII, args=(x,y))
        
        #now, calculating y-values and extracting each individual peak from multi-peak fit
        yNII0 = func(x, optimhANII[0], optimhANII[1], optimhANII[2])
        yhA = func(x, optimhANII[3], optimhANII[4], optimhANII[5])
        yNII = func(x, optimhANII[6], optimhANII[7], optimhANII[8])
                                 
        #yhA = func(x, optimhANII[0], optimhANII[1], optimhANII[2])
        #yNII = func(x, optimhANII[3], optimhANII[4], optimhANII[5])

        ySIIa = func(x, optimSII[0], optimSII[1], optimSII[2])
        ySIIb = func(x, optimSII[3], optimSII[4], optimSII[5])
        
        #calculating y values for a gaussian fit using new input paramters optimized above
        yhB = func(x, optimhB[0], optimhB[1], optimhB[2])
        yOIII = func(x, optimOIII[0], optimOIII[1], optimOIII[2])
        yOI = func(x, optimOI[0], optimOI[1], optimOI[2])

        """#comment below out after this run
        yhA = func(x, optimhA[0], optimhA[1], optimhA[2])
        yNII = func(x, optimNII[0], optimNII[1], optimNII[2])"""      

        
        Hbeta = integrate.simps(yhB, x)
        OIII = integrate.simps(yOIII, x)
        OI = integrate.simps(yOI, x)
        Halpha = integrate.simps(yhA, x)
        NII = integrate.simps(yNII, x)
        #SII = integrate.simps(ySIIa, x) + integrate.simps(ySIIb, x) #not sure which of these two it is (are they different?)
        SII = integrate.simps(ySII, x)

        f = plt.figure(2, figsize=(16,5))
        plt.clf()
        plt.xlim(4700, 7000)
        plt.grid()
        plt.tight_layout()
        plt.title('Gaussian Fit')

        plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
        plt.plot(x, y, c='k', linewidth=0.5)
        
        plt.plot(x, yhB, c='b', linewidth=0.75)
        plt.plot(x, yOIII, c='b', linewidth=0.75)
        plt.plot(x, yOI, c='b', linewidth=0.75)
        plt.plot(x, yhANII, c='b', lw=0.75)
        plt.plot(x, yhA, c='r', lw =0.75)
        plt.plot(x, yNII0, c='b', lw =0.75)
        plt.plot(x, yNII, c='b', lw =0.75)
        plt.plot(x, ySII, c='b', lw=0.75)
        plt.plot(x, ySIIa, c='b', lw=0.75)
        plt.plot(x, ySIIb, c='b', lw=0.75)

        #plot residuals

        plt.plot(x, ( (yhA+yNII0+yNII+ySIIa+ySIIb+yhB+yOIII+yOI) - y) -5, c='k', lw=0.75)
        
        f.show()
        raw_input()
        
    elif tooNoisy==True:
        print("spectrum is too noisy to analyze automatically, please check by hand")

    #elineFluxes = namedtuple('EmLines', 'hAlpha NII hBeta OIII SII OI noise hAisClear NIIisClear hBisClear OIIIisClear SIIisClear OIisClear tooNoisy')
    #emLines = elineFluxes(Halpha, NII, Hbeta, OIII, SII, OI, noise, hAisClear, NIIisClear, hBisClear, OIIIisClear, SIIisClear, OIisClear, tooNoisy)
    #return emLines 
        

def funcToDeleteItems(fullPathToDir):
    shutil.rmtree(fullPathToDir)
                
def run_code():
    #game plan: plot two different test specs without multiplying by fobs and see if they align, then try adding them
    stackDir = '/Users/plira/india_temp/variable_spectra/auto/stacking/stack_dir'
    stackArray = []
    """for filename in os.listdir(stackDir):
        if filename=='.DS_Store':
            continue
        else:
            stackArray.append(stackDir + '/' + filename + '/' + filename + '_parsed.txt')
            print(stackDir + '/' + filename + '/' + filename + '_parsed.txt')

    stackArray = ['/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/table_cross_match/BROADLINE/spec-0275-51910-0438/spec-0275-51910-0438_parsed.txt']
    print(str(len(stackArray)))"""
    stackArray = ['/Users/plira/india_temp/variable_spectra/auto/spectral_classifications_final/table_cross_match/BROADLINE/spec-0275-51910-0438/spec-0275-51910-0438_parsed.txt']


    gaussianFits(stackArray)

    plt.show()

if __name__ == '__main__':
    run_code()


"""PURGATORY


        
    #plt.vlines(6564, -20, 60, color='b')

    #trying to fit BLR
    guessBF = [0.9, 6564, 8]
    yBF = func(overallX, guessBF[0], guessBF[1], guessBF[2])
    yTest = originalStack - yBF
    
    guessTest = [7.5, 6564, 3]
    optimyTest, flag = sp.leastsq(errfunc, guessTest, args=(overallX,yTest))
    yNew = func(overallX, optimyTest[0], optimyTest[1], optimyTest[2])
    
    guess2 = [0.9, 6564, 8, 7.5, 6564, 3]
    optim2, flag = sp.leastsq(compositeErrFunc, guess2, args=(overallX, originalStack))
    print(optim2)
    
    #y2 = double_gaussian(overallX, optim2[0], optim2[1], optim2[2], optim2[3], optim2[4], optim2[5])
    y2broad = func(overallX, optim2[0], optim2[1], optim2[2])
    y2narrow = func(overallX, optim2[3], optim2[4], optim2[5])
    #plt.plot(overallX, y2broad, c='g', lw=1.0)
    #plt.plot(overallX, y2narrow, c='g', lw=1.0)

    #plt.plot(overallX, yBF, c='g', lw=1.5)
    #plt.plot(overallX, yNew, c='g', lw=1.5)
    #plt.plot(overallX, yBF + yNew, c='r', lw=1.5)"""
