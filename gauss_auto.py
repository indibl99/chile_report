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
    datapathIn = '/Users/plira/india_temp/variable_spectra/auto/temp_out/' + base + '_output.txt'
    datapathOut = '/Users/plira/india_temp/variable_spectra/auto/output_models/' + base + '/'
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
            elif i >= 216:
                starModel.writelines(line)    
    return fobs

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
    plt.xlim(3500,7700)
    plt.grid()
    plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
    plt.xlabel("wavelength (angstroms)")
    plt.ylabel("flux")
    title = "Spectrum and Subtracted STARLIGHT Model of \n" + filename
    plt.title(title)
    plt.plot(wave, inSpec, c='k', linewidth=0.4)
    plt.plot(wave, outSpec, c='r', linewidth=0.4)
    plt.plot(wave, subSpec, c='b', linewidth=0.4)
    m.show()
    m.savefig(outLoc + 'starlight_model.png',bbox_inches='tight')
    #pickleloc = outLoc + 'starlight_model.pickle'
    #pl.dump(m,file(pickleloc, 'wb'))

    OutData = namedtuple('outData', 'wave inSpec outSpec subSpec')
    data = OutData(wave, inSpec, outSpec, subSpec)
    return data

def gaussianFits(fileName, fobs, outLoc):

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
    
    #may need to make noise smarter soon - this is a test!!
    #NEED TO MAKE NOISE MUCH SMARTER
    noiseX = x[1450:1460]
    noiseY = y[1450:1460]
    
    noise = 3 * integrate.simps(noiseY, noiseX)
    
    #trying to fix the problem with the gaussian fits using peak utils indexes
    #automating guessing heights
    peaks = indexes(y, thres=0.3, min_dist=8.0) #these values are somewhat arbitrary right now
    numPeaks1 = 0
    peak1 = [0]
    peakIndices1 = [0]
    numPeaks2 = 0
    peak2 = [0]
    peakIndices2 = [0]
    peakX = x[peaks]
    peakY = y[peaks]

    hBisClear = False
    OIIIisClear = False
    OIisClear = False
    hAisClear = False
    NIIisClear = False
    SIIisClear = False

    Halpha = 0
    NII = 0
    Hbeta = 0
    OIII = 0
    SII = 0
    OI = 0
    
    
    for i in range(len(peakX)):
        if peakX[i] > 6500 and peakX[i] < 6780: #H-alpha, NII, and SII
            numPeaks2 = numPeaks2 + 1
            peak2.append(peakX[i])
            peakIndices2.append(i)
            if peakX[i] > 6560 and peakX[i] < 6570:
                hAisClear = True
            if peakX[i] > 6575 and peakX[i] < 6590:
                NIIisClear = True
            if peakX[i] > 6710 and peakX[i] < 6780:
                SIIisClear = True
        if peakX[i] > 4850 and peakX[i] < 4880: #H-beta
            numPeaks1 = numPeaks1 + 1
            peak1.append(peakX[i])
            peakIndices1.append(i)
            hBisClear = True
        if peakX[i] > 5000 and peakX[i] < 5015: #OIII
            numPeaks1 = numPeaks1 + 1
            peak1.append(peakX[i])
            peakIndices1.append(i)
            OIIIisClear = True
        if peakX[i] > 6290 and peakX[i] < 6310: #OI
            numPeaks1 = numPeaks1 + 1
            peak1.append(peakX[i])
            peakIndices1.append(i)
            OIisClear = True
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
  
    if numPeaks2==5:
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
    else:
        tooNoisy = True
        print("too few peaks detected, check this spectra by hand")

    if numPeaks1==2:
        amphB = peakY[peakIndices1[0]]
        ampOIII = peakY[peakIndices1[1]]
        ampOI = noise #this is just a place holder for now, means that .indexes couldn't find a peak at OI
    elif numPeaks1==3: #note! need to take care of uncertainties 
        amphB = peakY[peakIndices1[0]]
        ampOIII = peakY[peakIndices1[1]]
        ampOI = peakY[peakIndices1[2]]

    if hAisClear==False: tooNoisy = True
    if hBisClear==False and OIIIisClear==False: tooNoisy = True
    #ADD MORE PARAMETERS TO THIS FLAG TO ACCOUNT FOR SPECTRA THAT ARE JUST TOO NOISY AND NOT WORTH IT
    
    if tooNoisy==False: 
        guesshB = [amphB, 4865, width]
        guessOIII = [ampOIII, 5007, width]
        guessOI = [ampOI, 6300, width]
        guesshA = [amphA, 6563, width]
        guessNII = [ampNII, 6583, width]
        guessSIIa = [ampSIIa, 6716, width] #6716
        guessSIIb = [ampSIIb, 6731, width] #6731

        #using the least square function to optimize the paramters for the gaussin fit(the params for the func() function)
        
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
        

def BPT(emLines, outLoc, isCumulative):

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
    NIIisClear = emLines.NIIisClear
    hBisClear = emLines.hBisClear
    OIIIisClear = emLines.OIIIisClear
    SIIisClear = emLines.SIIisClear
    OIisClear = emLines.OIisClear
    
    flagY = 0 #-1 is down, 1 is up
    flagX = 0 #-1 is pointing to the left, 1 is pointing to the right
    allClear = True
    allClearNII = True
    allClearSII = True
    allClearOI = True
    
    
    if hAisClear==False: #right now there is a problem if both H-alpha and NII/SII/OI are both uncertain, same for OIII and HB
        Halpha = noise   #but if H-alpha is uncertain, the spectrum is too noisy anyways 
        flagX = 1
        print("hA is uncertain")
        allClear = False
    if NIIisClear==False:
        NII = noise
        flagX = -1
        print("NII is uncertain")
        allClearNII = False
    if hBisClear==False:
        Hbeta = noise
        flagY = 1
        print("hB is uncertain")
        allClear = False
    if OIIIisClear==False:
        OIII = noise
        flagY = -1
        print("OIII is uncertain")
        allClear = False
    if SIIisClear==False:
        SII = noise
        flagX = -1
        print("SII is uncertain")
        allClearSII = False
    if OIisClear==False:
        OI = noise
        flagX = -1
        print("OI is uncertain")
        allClearOI = False
        
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
    print('debugging BPTs...')
    print('x = ' + str(np.log(OI/Halpha)) + ' y = ' + str(np.log(OIII/Hbeta)))
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

    ax3.plot(np.log(OI/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarOI, markeredgecolor=colorVarOI)
    if allClear==False or allClearOI==False: ax3.quiver(np.log(OI/Halpha), np.log(OIII/Hbeta), flagX, flagY, color=colorVarOI)
    #k.savefig(outLoc + 'BPT_OI.png', bbox_inches='tight')

    if isCumulative==True:
        plt.show()
    else:
        plt.tight_layout()
        c.savefig(outLoc + 'BPT Diagrams', bbox_inches='tight')

    #raw_input()

    bptR = [bptNII, bptSII, bptOI]

    return bptR
                
    
def run_code():
    #grid file is going to already set up but for the future, we need to edit lines in the grid file in the loop
    #for a first run, fill the directories with like 5 spectra of different kinds before you do them all at once, and then when you do, make sure you have a couple hours
    #and save as you go, so that you can quit in the middle of the loop if need be 
    
    outpath = '/Users/plira/india_temp/variable_spectra/auto/output_models/'
    starlight = '/Users/plira/software/STARLIGHTv04/'
    oldFile = '5339079575614451712'
    totalFile = '/Users/plira/india_temp/variable_spectra/auto/cumulative_output.txt'
    totalBPT = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_BPT.txt', 'w')

    cumulative = open(totalFile, 'w')
    cumulative.write("{: <20} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <15} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <15} {: <15} {: <15}".format('Filename','hAlpha','NII','hBeta','OIII','SII','OI','noise','hAisClear','NIIisClear','hBisClear','OIIIisClear','SIIisClear','OIisClear','tooNoisy', 'BPT NII', 'BPT SII', 'BPT OI'))
        
    for filename in os.listdir('/Users/plira/india_temp/variable_spectra/auto/temp_input'): #right now, this directory has the models, but eventually it will just have raw data
        print('filename: ' + filename)
        base = filename[:-4]
        print("base: " + base)
        #you are so fucking good
        #temporarily commenting out starlight
        grid = '/Users/plira/software/STARLIGHTv04/grid_temp.in'
        for line in fileinput.input([grid], inplace=True):
            line = line.replace(oldFile, base)
            sys.stdout.write(line)

        os.chdir(starlight)
        call(["pwd"])
        os.system("./Starlight_v04_Mac.exe < grid_auto.in") #:)"""
        #the way this is set up, if you just comment out the STARLIGHT section, it should run without running starlight? check that to be sure
        
        fobs = parse_datafile(base, filename)
        locI = outpath + base + '/'
        locO = '/Users/plira/india_temp/variable_spectra/auto/temp_input/' + filename
        fileOut = open(locI + 'analysis output', 'w')
        print("locI: " + locI)
        shutil.copy2(locO, locI)
        shutil.copy2('/Users/plira/india_temp/variable_spectra/auto/temp_out/'+base+'_output.txt', locI)
        loc = outpath + base + '/' + base + '_parsed.txt' #this is the PARSED model that gets passed
        emLines = gaussianFits(loc, fobs, locI);
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
            bptO = BPT(emLines, locI, False)
            fileOut.write('\n ')
            fileOut.write("\n{: <15} {: <15} {: <15}".format('BPT NII', 'BPT SII', 'BPT OI'))
            fileOut.write("\n{: <15} {: <15} {: <15}".format(bptO[0], bptO[1], bptO[2]))
            cumulative.write("{: <15} {: <15} {: <15}".format(bptO[0], bptO[1], bptO[2]))
            totalBPT.write('\n'+ str(output[0]) + ',' + str(output[1]) + ',' + str(output[2]) + ',' + str(output[3]) + ',' + str(output[4]) + ',' + str(output[5]) + ',' + str(output[6]) + ',' + str(output[7]) + ',' + str(output[8]) + ',' + str(output[9]) + ',' + str(output[10]) + ',' + str(output[11]) + ',' + str(output[12]) + ',' + str(output[13]) + ',' + str(output[14]))

            #put dynamic plot animation here!
        elif emLines.tooNoisy==True:
            notWorthIt = open(locI + 'ErrorMessage_noisy_spectrum', 'w')
            notWorthIt.write('spectrum was too noisy, please inspect file ' + base + '.txt by hand')
            notWorthIt.close()

        fileOut.close()
        oldFile = base
        print("oldFile: " + oldFile)

    cumulative.close()
    totalBPT.close()

    #holy shit this works


if __name__ == '__main__':
    run_code()

#ok game plan - first of all, you did a fucking amazing job today
#right now as is, this works
#tomorrow, you need to try and get the executable file with starlight working and figure out how to parse into and edit the grid file
#also, right now the arrows in BPT diagrams don't mean much and uncertainty in general needs work
#specifically and importantly, you need to have a cop-out flag for those spectra that are just so noisy you're just going to get crazy data if you try to analyze
#maybe consider implementing a small S/N calculation - otherwise just use the uncertainty flags you already have 
#you also need to make sure you have a file keeping track of the spectra as you go
#and forget the spobjid, you decided to jst change your filename system

#ALSO you need to find a way to save the graphs as you go, or run them individualy (could just change the path) but would be easier to save the graphs as you go 



#game plan pt 2
#ANOTHER crazy productive day, so proud of you
#almost everything is ready, now you're working on graphing the cumulative BPT diagram results (need to have a way for the function to skip over
#the files we said were too noisy
#you're kick ass I love it
#also maybe try your automated fit against one of the handmade ones and see how they compare!
#test a file of like 5 before all of them at once
#also maybe print out a lot of analysis and put it in a text file so you have a little more detailed info for each? 
#also want to be able to return some statistics on your runs, so be able to see automatically where a point falls

#need to catch runtime warnings and figure out why it wants you to press enter in btw valid runs - rawinput?

#MAKE BPTS SUBPLOTS AND STILL FIXING CUMULATIVE THING
#ideas from last night - color code with if greater than statements, save properly fitted BPT subplots 
