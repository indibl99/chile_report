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
import math

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def double_gaussian(x, a1, x01, sigma1, a2, x02, sigma2):
    return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2)

def triple_gaussian(x, a1, x01, sigma1, a2, x02, sigma2, a3, x03, sigma3):
    return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2) + func(x, a3, x03, sigma3)

"""def quad_gaussian(x, a1, x01, sigma1, a2, x02, sigma2, a3, x03, sigma3, a4, x04, sigma4): #need fixed narrow line comp widths
    #return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2) + func(x, a3, x03, sigma3) + func(x, a4, x04, sigma4)
    return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2) + func(x, a3, x03, sigma1) + func(x, a4, x04, sigma1)"""

def quad_gaussian(x, a1, x01, sigma1, a2, x02, sigma2, a3, x03, a4, x04): #need fixed narrow line comp widths
    return func(x, a1, x01, sigma1) + func(x, a2, 6564, sigma2) + func(x, a3, 6564, sigma1) + func(x, a4, x04, sigma1) #try fiddling with which sigma?

def composite_gaussian(x, a1, x01, sigma1, a2, x02, sigma2):
    return func(x, a1, x01, sigma1) + func(x, a2, x01, sigma2) #same x01

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

"""def quadErrFunc(guess, xData, yData):#need fixed narrow line comp widths
    yFit = quad_gaussian(xData, guess[0], guess[1], guess[2], guess[3], guess[4], guess[5], guess[6], guess[7], guess[8], guess[9], guess[10], guess[11])
    err4 = yData - yFit
    return err4"""

def quadErrFunc_fix_width(guess, xData, yData):#need fixed narrow line comp widths
    yFit = quad_gaussian(xData, guess[0], guess[1], guess[2], guess[3], guess[4], guess[5], guess[6], guess[7], guess[8], guess[9])
    err4 = yData - yFit
    return err4


def compositeErrFunc(guess, xData, yData):
    blrFit = func(xData, guess[0], guess[1], guess[2])
    tempY = yData - blrFit #remaining spec after guess BLR component is taken away
    nlrFit = func(xData, guess[3], guess[4], guess[5])
    errOverall = yData - (blrFit + nlrFit)
    return errOverall

#def blr_gaussian(x, a1, x01, sigma1, a2, x02, sigma2):


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
    thisWave = []
    thisSpec = []
    for spectrum in filenames:
        """if spectrum[0:4]=='Blind':
            #spec = np.genfromtxt(spectrum)
            print(spectrum)
        else:
            spec = np.genfromtxt(spectrum)"""
        spec = np.genfromtxt(spectrum)
        
        wave = spec[:,0] #no fobs
        inSpec = spec[:,1]
        outSpec = spec[:,2]
        subSpec = inSpec - outSpec

        thisWave.append(wave)
        thisSpec.append(subSpec)

        stack = stack + subSpec #needs to be smarter than this, and use gaussian to align the centers
        waveStack = wave

        plt.xlim(3700,7700)
        plt.grid()
        plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
        plt.xlabel("wavelength (angstroms)")
        plt.ylabel("flux")
        
        plt.plot(wave, subSpec, c='r', linewidth=0.4)
        

    plt.plot(waveStack, stack, c = 'k', linewidth=0.4)
    plt.show()
    #raw_input()

    data = [waveStack, stack]
    np.savetxt('/Users/plira/india_temp/variable_spectra/auto/stacking/save_stack.txt', np.c_[waveStack, stack])

    OutData = namedtuple('data', 'wave spec thisWave thisSpec')
    data = OutData(waveStack, stack, thisWave, thisSpec)
    return data
    


def gaussianFits(fileName):
    
    data = loadData(fileName)
    stack = data.thisSpec

    hAGauss = []
    NIIGauss = []
    NII0Gauss = []
    original = []
    Wave = data.wave
    
    for spectrum in stack:
        #x = data.wave
        #y = data.spec
        x = data.wave
        y = spectrum
        original.append(y)
        
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
            if peakX[i] > 6713 and peakX[i] < 6722:
                ampSIIa = peakY[i]
                #print('SIIa (' + str(peakX[i]) + ', ' + str(peakY[i]) + ')')
            if peakX[i] > 6727 and peakX[i] < 6733:
                ampSIIb = peakY[i]
                #print('SIIb (' + str(peakX[i]) + ', ' + str(peakY[i]) + ')')
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
        #print("numPeaks2 = " + str(numPeaks2))
        #print("numPeaks1 = " + str(numPeaks1))

        #print(peak2)
        #print(peak1)

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

        #centers needed for stacking!!
        print("optimal value for center of H-alpha gaussian: " + str(optimhANII[4]))
        centerhA.append(optimhANII[4])
        
        hAGauss.append(yhA)
        NIIGauss.append(yNII)
        NII0Gauss.append(yNII0)
                                 
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

        #f = plt.figure(2, figsize=(16,5))
        #plt.clf()
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

        #plt.plot(x, yBF, c='g', lw=0.75)
        #plt.plot(x, yNew, c='g', lw=0.75)
        #plt.plot(x, yBF + yNew, c='r', lw=0.75)
        #plt.plot(x, y - yBF - yNew - 5, c='#999999', lw=0.75)

        #plot residuals

        plt.plot(x, ( (yhA+yNII0+yNII+ySIIa+ySIIb+yhB+yOIII+yOI) - y) -20, c='k', lw=0.75)
    
    plt.show()
    #raw_input()
    return hAGauss, NIIGauss, NII0Gauss, original

def funcToDeleteItems(fullPathToDir):
    shutil.rmtree(fullPathToDir)


def resample(xarr1, data1):
    xarr = xarr1.copy()
    data = data1.copy()

    xarr[0] = math.ceil(xarr1[0])
    data[0] = ( (xarr[0]* (data1[1]-data1[0])) + (xarr1[1]*data1[0]) - (xarr1[0]*data1[1])) / (xarr1[1] - xarr1[0])

    endIndex = len(xarr1)-1

    for i, wave in enumerate(xarr1, 1): #need to start at index = 1
        if i >= endIndex:
            break
        else:
            x1 = xarr1[i]
            x2 = xarr1[i+1] #store these two values for later calculation

            y1 = data1[i]
            y2 = data1[i+1] #stored surrounding y values

            xarr[i] = xarr[i-1] + 1.0 #now xarr1 is evenly spaced, xarr1[wave
            interpFlux = ( (xarr[i]* (y2-y1)) + (x2*y1) - (x1*y2) ) / (x2 - x1)
            data[i] = interpFlux
    #3502.0
    #3499.0
    #print('first element before filling in: ' + str(xarr[0]))
    if xarr[0] != 3500.0:
        ogStart = xarr[0]
        #print('ogStart: ' + str(ogStart))
        offset = 3500.0 - xarr[0]
        if offset < 0:
            tempInd = ogStart - 1
            while tempInd != 3499.0:
                #print('tempInd: ' + str(tempInd))
                #xarr.insert(0, tempInd) #add another index
                xarr = np.insert(xarr, 0, tempInd)
                #xarr.pop(len(xarr)-1) #remove the last element
                xarr = np.delete(xarr,len(xarr)-1)
                tempInd = tempInd - 1

                #data.insert(0, 0.0)
                #data.pop((len(xarr)-1))
                data = np.insert(data, 0, 0.0)
                data = np.delete(data,len(data)-1)
                
        if offset > 0:
            tempInd = ogStart
            tempEnd = xarr[len(xarr)-1]
            while tempInd != 3500.0:
                #xarr.pop(0) #remove first element, getting closer to 3500
                xarr = np.delete(xarr,0)
                #xarr.insert(len(xarr)-1, tempEnd + 1)
                xarr = np.insert(xarr, len(xarr)-1, tempEnd + 1)
                tempInd = tempInd + 1
                tempEnd = tempEnd + 1

                #data.pop(0)
                data = np.delete(data,0)
                #data.insert(len(xarr)-1, 0.0)
                data = np.insert(data, len(data)-1, 0.0)
    #print('first element after filling in: ' + str(xarr[0]))
                
            
    #del xarr[len(xarr)-1]
    #del data[len(data)-1]
    #ignore the last value (should still just be 0)

    return xarr, data

def stack(hAGauss, NIIGauss, NII0Gauss, original):
    print(centerhA)
    originalStack = np.zeros(len(hAGauss[0]))
    alignedStack = np.zeros(len(hAGauss[0]))
    overallX = np.zeros(len(hAGauss[0]))
    NIIStack =np.zeros(len(hAGauss[0]))
    hAStack =np.zeros(len(hAGauss[0]))

    saveOriginal = []
    
    #print("wave " + str(Wave[0]))
    i = 0
    for hAspec in hAGauss:
        print('hAGauss: ' + str(len(hAGauss)))
        peaks = indexes(hAspec)
        peakX = Wave[peaks]
        peakY = hAspec[peaks]
        xSpec = Wave.copy()
        NIIspec = NIIGauss[i].copy()
        NII0spec = NII0Gauss[i].copy()
        ogSpec = original[i].copy()

        #first need to resample each spec so that they have an index where the gaussian function's peak is...
        """peaks0 = indexes(hAspec)
        peaksX0 = xSpec[peaks0]
        peaksY0 = hAspec[peaks0]
        print('look here: ' + str(peaksX0[0]))
        preShift = peaksX0[0] - centerhA[i] #how much the peaks are offset from their indices

        xSpec = xSpec + preShift

        xSpechA, hAspec = resample(xSpec, hAspec)
        xSpecNII, NIIspec = resample(xSpec, NIIspec)
        xSpecNII0, NII0spec = resample(xSpec, NII0spec)
        xSpecOg, ogSpec = resample(xSpec, ogSpec)

        peaks0 = indexes(hAspec)
        peaksX0 = xSpechA[peaks0]
        peaksY0 = hAspec[peaks0]
        print('now, centerhA is ' + str(centerhA[i]) + ' and the index of the peak of hAspec is ' + str(peaksX0[0]))"""

        if centerhA[i]==6564.0: #target center
            print('already centered')
        elif centerhA[i]!=6564.0:
            #peaksTest = indexes(hAspec)
            #peaksTestX = xSpec[peaksTest]
            #peaksTestY = hAspec[peaksTest]
            #print('not shifted hA index: ' + str(peaksTestX[0]))
            print('current target center of hA: ' + str(centerhA[i]))
            
            shift = 6564.0 - centerhA[i]
            print('shift: ' + str(shift))
            xSpec = xSpec + shift
            #peaksTest = indexes(hAspec)
            #peaksTestX = xSpec[peaksTest]
            #peaksTestY = hAspec[peaksTest]
            #print('shifted hA index: ' + str(peaksTestX[0]))
            
        #need to resample the entire function so that the peak of h-alpha rests on 6564, and then can add them all with that aligned, filling in the ends

        print('first val shifted, not resampled: ' + str(xSpec[0]))
        xSpechA, hAspec = resample(xSpec, hAspec) #need to return these!!
        print('first val shifted, resampled: ' + str(xSpechA[0]))
        xSpecNII, NIIspec = resample(xSpec, NIIspec)
        xSpecNII0, NII0spec = resample(xSpec, NII0spec)
        xSpecOg, ogSpec = resample(xSpec, ogSpec)
        
        overallX = xSpecOg
        #need to chop off ends of all to account for uneven sampling -- done
        
        peaks = indexes(NII0spec)
        peakX = xSpecNII0[peaks]
        peakY = NII0spec[peaks]
        print('after both shifts and resampling, the peak of the NII0 gaussian is at index ' + str(peakX[0]) + ' of xSpecNII0')
        #print('first val: ' + str(xSpec[0]))
        
        #plt.clf()
        plt.xlim(4700, 7000)
        plt.grid()
        #print('len(xSpec): ' + str(len(xSpec)))
        #print('len(hAspec): ' + str(len(hAspec)))
        #plt.plot(xSpechA, hAspec, color='m')
        #plt.plot(xSpecNII, NIIspec, color='m')
        #plt.plot(xSpecNII0, NII0spec, color='m')
        plt.plot(xSpecOg, ogSpec, color ='#999999')
        saveOriginal.append(ogSpec)
        
        alignedStack = alignedStack + hAspec + NIIspec + NII0spec
        originalStack = originalStack + ogSpec
        NIIStack = NIIStack + NIIspec + NII0spec
        hAStack = hAStack + hAspec
        #also need to deal with residuals!!
        
        i = i + 1

    #plotting, was r and k and r
    
    #plt.plot(overallX, alignedStack, color='b')
    plt.plot(overallX, originalStack, color='k')
    #plt.plot(overallX, NIIStack, color='b')
    #plt.plot(overallX, hAStack, color='b')

    #amp NII0, center NII0, width NII0 (and all others), amp broad, center broad, width broad, amp Ha, center hA, amp NII, center NII
    guessQuad = [33.9, 6549, 4, 7.0, 6564, 12, 273.5, 6564, 105.5, 6584]
    guessQuad = [9.9, 6549, 4, 7.0, 6564, 12, 130.4, 6564, 27.3, 6584]
    #guessQuad = [1.2, 6549, 2.8, 0.9, 6564, 8, 7.7, 6564, 3.5, 6584]

    #I want all the widths to be the same, even if it means they get bigger as needed and the BLR shrinks
    optim4, flag = sp.leastsq(quadErrFunc_fix_width, guessQuad, args=(overallX,originalStack))
    #lambda guessQuad[0], guessQuad[1], guessQuad[3], guessQuad[4], guessQuad[5], guessQuad[6], guessQuad[7], guessQuad[9], guessQuad[10], xData, yData: quadErrFunc(guessQuad[0], guessQuad[1], 2.8, guessQuad[3], guessQuad[4], guessQuad[5], guessQuad[6], guessQuad[7], 2.8, guessQuad[9], guessQuad[10], 2.8, xData, yData)
    
    print('look here for widths!')
    print(optim4)
    """yQuad = quad_gaussian(overallX, optim4[0], optim4[1], optim4[2], optim4[3], optim4[4], optim4[5], optim4[6], optim4[7], optim4[8], optim4[9], optim4[10], optim4[11])    
    yQuadBroad = func(overallX, optim4[3], optim4[4], optim4[5])
    yQuadNarrow = func(overallX, optim4[6], optim4[7], optim4[8])
    yQuadNIIa = func(overallX, optim4[0], optim4[1], optim4[2])
    yQuadNIIb = func(overallX, optim4[9], optim4[10], optim4[11])"""

    yQuad = quad_gaussian(overallX, optim4[0], optim4[1], optim4[2], optim4[3], optim4[4], optim4[5], optim4[6], optim4[7], optim4[8], optim4[9])    
    yQuadBroad = func(overallX, optim4[3], optim4[4], optim4[5])
    yQuadNarrow = func(overallX, optim4[6], optim4[7], optim4[2])
    yQuadNIIa = func(overallX, optim4[0], optim4[1], optim4[2])
    yQuadNIIb = func(overallX, optim4[8], optim4[9], optim4[2])


    plt.plot(overallX, yQuadBroad, c='r', lw=1.5, linestyle='--')
    plt.plot(overallX, yQuadNarrow, c='r', lw=1.5)
    plt.plot(overallX, yQuadNIIa, c='b', lw=1.5)
    plt.plot(overallX, yQuadNIIb, c='b', lw=1.5)
    
    
    #plt.plot(overallX, alignedStack - originalStack - 5, color='#999999')
    #plt.plot(overallX, alignedStack - yBF - yNew - 5, color='#999999')
    plt.show()
    raw_input()

    return overallX, saveOriginal, originalStack, yQuadBroad, yQuadNarrow, yQuadNIIa, yQuadNIIb, guessQuad
                
def run_code():
    #game plan: plot two different test specs without multiplying by fobs and see if they align, then try adding them
    ask = raw_input("Do you want to run a new sample (y/n)? ")
    internalID = raw_input('Please enter the name for this run: ')
    cmfont = {'fontname':'Georgia'}
    #plt.rcParams.update(cmfont)

    if ask=='y':
        stackDir = '/Users/plira/india_temp/variable_spectra/auto/stacking/' + internalID
        saveDir = '/Users/plira/india_temp/variable_spectra/auto/stacking/saved_runs/' + internalID + '/'
        if os.path.exists(saveDir):
            shutil.rmtree(saveDir)
        os.makedirs(saveDir)
        copy_tree(stackDir + '/', saveDir)
        
        stackArray = []
        for filename in os.listdir(stackDir):
            if filename=='.DS_Store':
                continue
            else:
                print(filename[0:5])
                if filename[0:5]=='Blind':
                    shortFilename = filename[-20:]
                    print(shortFilename)
                    stackArray.append(stackDir + '/' + filename + '/' + shortFilename + '_parsed.txt')
                else:
                    stackArray.append(stackDir + '/' + filename + '/' + filename + '_parsed.txt')
                print(stackDir + '/' + filename + '/' + filename + '_parsed.txt')
        
        print(str(len(stackArray)))
        #stackArray = ['/Users/plira/india_temp/variable_spectra/auto/stacking/test_dir/spec-0267-51608-0246_parsed.txt', '/Users/plira/india_temp/variable_spectra/auto/stacking/test_dir/spec-0267-51608-0493_parsed.txt']

        #loadData(stackArray)
        #np.savetxt('/Users/plira/india_temp/variable_spectra/auto/stacking/save_stack.txt', np.c_[waveStack, stack])
        
        hAGauss, NIIGauss, NII0Gauss, original = gaussianFits(stackArray)
        overallX, saveOriginal, originalStack, yQuadBroad, yQuadNarrow, yQuadNIIa, yQuadNIIb, guessQuad = stack(hAGauss, NIIGauss, NII0Gauss, original)
        np.savetxt(saveDir + internalID + '_save_stack.txt', np.c_[overallX, originalStack, yQuadBroad, yQuadNarrow, yQuadNIIa, yQuadNIIb])
        
        
        tempArray = []
        #tempArray.append(overallX)
        for individual in saveOriginal:
            tempArray.append(individual)        
        np.savetxt(saveDir + internalID + '_save_individuals.txt', np.c_[tempArray])

        paramFile = open(saveDir + 'blr_params.txt', 'w')
        for guess in guessQuad:
            paramFile.write(str(guess) + ', ')
        paramFile.close()
        
    elif ask=='n':
        #print("work in progress")
        loc = raw_input("Please enter filepath of saved run: ")
        saved = np.genfromtxt(loc)
        individuals = np.genfromtxt(loc[:-9] + 'individuals.txt')
        
        x = saved[:,0]
        originalStack = saved[:,1]
        yQuadBroad = saved[:,2]
        yQuadNarrow = saved[:,3]
        yQuadNIIa = saved[:,4]
        yQuadNIIb = saved[:,5]

        for spec in individuals:
            plt.plot(x, spec, c='#999999')

        plt.title('stacked ' + internalID, **cmfont)
        plt.grid()
        plt.plot(x, originalStack, c='k', lw=0.75)
        plt.plot(x, yQuadBroad + yQuadNarrow, c='b', lw=1.5, alpha=0.75)
        plt.plot(x, yQuadBroad, c='r', lw=1.5, linestyle='--')
        plt.plot(x, yQuadNarrow, c='r', lw=1.5)
        plt.plot(x, yQuadNIIa, c='b', lw=1.5)
        plt.plot(x, yQuadNIIb, c='b', lw=1.5)
        
        plt.show()

        raw_input()
        
        

centerhA = []
Wave = np.arange(3500, 7701, 1)
cmfont = {'fontname':'Georgia'}
if __name__ == '__main__':
    run_code()

