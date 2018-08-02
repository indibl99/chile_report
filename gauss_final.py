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

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def triple_gaussian(x, a1, x01, sigma1, a2, x02, sigma2, a3, x03, sigma3):
    return func(x, a1, x01, sigma1) + func(x, a2, x02, sigma2) + func(x, a3, x03, sigma3)

def errfunc(guess, xData, yData):
    return (func(xData, guess[0], guess[1], guess[2]) - yData)

def loadData(filename, fobs):
    loc = filename
    spec = np.genfromtxt(filename)
    
    wave = spec[:,0]
    inSpec = fobs*spec[:,1]
    outSpec = fobs*spec[:,2]
    subSpec = inSpec - outSpec

    ####    EMISSION LINES WE ARE LOOKING FOR      ###
    
    eNII = np.empty(100)
    ehA = np.empty(100)
    ehB = np.empty(100)
    eOIII = np.empty(100)
    eSIIa = np.empty(100)
    eSIIb = np.empty(100)
    eOI = np.empty(100)
    for i in range(100):
        eNII[i] = 6584
        ehA[i] = 6563
        ehB[i] = 4860
        eOIII[i] = 5007
        eSIIa[i] = 6716
        eSIIb[i] = 6731
        eOI[i] = 6300

    """eNII = np.full( (100, 1), 6584)
    ehA = np.full( (100, 1), 6563)
    ehB = np.full( (100, 1), 4860)
    eOIII = np.full( (100, 1), 5007) #the equations you use only care about OIII 5007
    eSIIa = np.full( (100, 1), 6716)
    eSIIb = np.full( (100, 1), 6731)
    eOI = np.full( (100,1), 6300)"""

    eY = np.arange(-10, 90, 1)

    ##plot emission lines
    j = plt.figure(1)
    
    plt.xlabel("wavelength (angstroms)")
    plt.ylabel("flux")
    title = "Subtracted STARLIGHT Model of \n" + filename
    plt.title(title)
    
    plt.plot(eNII, eY, c='r', linewidth=3.0, alpha=0.5)
    plt.text(6586, 30, "NII")
    plt.plot(ehA, eY, c='r', linewidth=3.0, alpha=0.5)
    plt.text(6566, 30, "H-alpha")
    plt.plot(ehB, eY, c='r', linewidth=3.0, alpha=0.5)
    plt.text(4865, 30, "H-beta")
    plt.plot(eOIII, eY, c='r', linewidth=3.0, alpha=0.5)
    plt.text(5010, 30, "OIII")
    plt.plot(eSIIa, eY, c='r', linewidth=3.0, alpha=0.5)
    plt.text(6720, 30, "SII")
    plt.plot(eSIIb, eY, c='r', linewidth=3.0, alpha=0.5)
    plt.text(6735, 25, "SII")
    plt.plot(eOI, eY, c='r', linewidth=1.3, alpha=0.5)
    plt.text(6305, 30, "OI")
    
    plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
    plt.plot(wave, subSpec, c='k', linewidth=0.5) #plot subSpec and ask for user input for guesses

    j.show()
    
    #plot observed spectrum, STARLIGHT model, subtracted model

    m = plt.figure(2)
    plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
    plt.xlabel("wavelength (angstroms)")
    plt.ylabel("flux")
    title = "Spectrum and Subtracted STARLIGHT Model of \n" + filename
    plt.title(title)
    plt.plot(wave, inSpec, c='k', linewidth=0.4)
    plt.plot(wave, outSpec, c='r', linewidth=0.4)
    plt.plot(wave, subSpec, c='b', linewidth=0.4)
    m.show()

    OutData = namedtuple('outData', 'wave inSpec outSpec subSpec')
    data = OutData(wave, inSpec, outSpec, subSpec)
    return data

def gaussianFits(reRun, oldInput, newInput, output):
    
    if reRun=='yes':
        fileName = raw_input("Please enter the file you're trying to analyze: ")
        #fileName = oldInput[0]
        fobs = float(oldInput[1])
        data = loadData(fileName, fobs)
        width = float(oldInput[2])
        
        amphA = float(oldInput[3])
        ampNII = float(oldInput[4])
        ampSIIa = float(oldInput[5])
        ampSIIb = float(oldInput[6])
        
        ampOI = float(oldInput[7])

        amphB = float(oldInput[8])
        ampOIII = float(oldInput[9])

        xi = int(oldInput[10])
        xf = int(oldInput[11])

    else:
        fileName = raw_input("Please enter the file you're trying to analyze: ")
        fobs = raw_input("And now, the corresponding fobs value: ")
        fobs = float(fobs)
        data = loadData(fileName, fobs)

        print("Looking at the data, please enter your best guess for the following prompts...")
        width = float(raw_input("approximate width(sigma) of all lines: "))
        amphA = float(raw_input("amplitude of H-alpha? "))
        ampNII = float(raw_input("amplitude of NII? "))
        ampSIIa = float(raw_input("amplitude of the first SII? "))
        ampSIIb = float(raw_input("amplitude of the second SII? "))
        
        ampOI = float(raw_input("amplitude of OI? "))

        amphB = float(raw_input("amplitude of H-beta? "))
        ampOIII = float(raw_input("amplitude of OIII? "))

        #find a good, tall piece of noise, *3, and enter bounds here to be used in case of upperLimit
        print("In case of uncertainity, please fine a nice tall piece of noise and enter the limits when prompted.")
        #data file goes 3500 to 7700
        xi = int(raw_input("from x: ")) - 3500
        xf = int(raw_input("to x: ")) - 3500

        newInput.write(fileName)
        newInput.write('\n' + str(fobs) + '\n')
        newInput.write(str(width) + '\n')
        newInput.write(str(amphA) + '\n')
        newInput.write(str(ampNII) + '\n')
        newInput.write(str(ampSIIa) + '\n')
        newInput.write(str(ampSIIb) + '\n')
        newInput.write(str(ampOI) + '\n')
        newInput.write(str(amphB) + '\n')
        newInput.write(str(ampOIII) + '\n')
        newInput.write(str(xi) + '\n')
        newInput.write(str(xf) + '\n')

    #trying to fix the problem with the gaussian fits using peak utils indexes
    peaks = indexes(y, thres=0.3, min_dist=8.0) #these values are somewhat arbitrary right now
    print(peaks)
    print(x[peaks], y[peaks]) #want to only look at 6500-6800

    numPeaks = 0
    peak2 = [0]
    peakX = x[peaks]
    for i in range(len(peakX)):
        if peakX[i] > 6500 and peakX[i] < 6780:
            numPeaks = numPeaks + 1
            peak2.append(peakX[i])
    peak2.remove(0)
    print("numPeaks = " + str(numPeaks))
    print(peak2)
    
    """guesshA = [63.5, 6563, 3]
    guessNII = [32.8, 6583, 3]

    guesshB = [16.9, 4865, 3]
    guessOIII = [4.1, 5009, 3]

    guessSIIa = [9.0, 6716, 3]
    guessSIIb = [7.1, 6736, 3]

    guessOI = [2.5, 6300, 3]"""

    guesshA = [amphA, 6563, width]
    guessNII = [ampNII, 6583, width]

    guesshB = [amphB, 4865, width]
    guessOIII = [ampOIII, 5007, width]

    guessSIIa = [ampSIIa, 6716, width] #6716
    guessSIIb = [ampSIIb, 6731, width] #6731

    guessOI = [ampOI, 6300, width]

    x = data.wave #could also write as hAregion[:,0]
    y = data.subSpec

    #using the least square function to optimize the paramters for the gaussin fit(the params for the func() function)
 
    optimhA, flag = sp.leastsq(errfunc, guesshA, args=(x, y))
    optimNII, flag = sp.leastsq(errfunc, guessNII, args=(x, y))

    optimhB, flag = sp.leastsq(errfunc, guesshB, args=(x, y))
    optimOIII, flag = sp.leastsq(errfunc, guessOIII, args=(x, y))

    optimSIIa, flag = sp.leastsq(errfunc, guessSIIa, args=(x, y))
    optimSIIb, flag = sp.leastsq(errfunc, guessSIIb, args=(x, y))

    optimOI, flag = sp.leastsq(errfunc, guessOI, args=(x, y))

    #calculating y values for a gaussian fit using new input paramters optimized above
    y2 = func(x, optimhA[0], optimhA[1], optimhA[2])
    y3 = func(x, optimNII[0], optimNII[1], optimNII[2])

    y4 = func(x, optimhB[0], optimhB[1], optimhB[2])
    y5 = func(x, optimOIII[0], optimOIII[1], optimOIII[2])

    y6 = func(x, optimSIIa[0], optimSIIa[1], optimSIIa[2])
    y7 = func(x, optimSIIb[0], optimSIIb[1], optimSIIb[2])

    y8 = func(x, optimOI[0], optimOI[1], optimOI[2])

        
    
    #TOMORROW MORNING HERE'S WHERE YOU'RE AT: KIND OF GAVE UP ON VIRTUAL ENVIRONMENTS FOR NOW (FOR NOW)
    #FINALLY GOT PEAKUTILS INSTALLED AND USING THE INDEXES FUNCTION TO TRY AND FIX THE GAUSSIAN FIT PROBLEMS
    #RIGHT NOW THE THRESHOLD IS A GOOD NUMBER, BUT TRYING TO CONSTRAIN TO THE H-ALPHA/NII - SII REGION
    #AND THEN MOVE THIS ABOVE THE GUESSES TO REPLACE USER INPUT WITH AUTOMATION (WILL EVENTUALLY HAPPEN TO ALL)
    #YOU DO HAVE ONE WORKING VIRTUAL ENV (A VENV) CALLED ENV
    #I'M ALSO PRETTY SURE YOU'RE STILL USING PYTHON 2.7 HERE
    #SEE PAGE IN JOURNAL FOR GAME PLAN, ONCE YOU'VE FIXED THE GAUSSIANS, THE IDEA IS TO MAKE A PIPELINE IN A PYTHON SCRIPT
    #WHICH WILL AUTOMATE WHAT YOU'VE BEEN DOING USING INDEXES TO GUESS PEAKS
    #MAYBE SET A FAIL SAFE FLAG THAT WILL ALERT YOU IF THE S/N RATIO IS TOO BAD TO AUTOMATICALLY ANALYZE
    #MOTIVATION, PART OF RATIONALE FOR DOING THIS IS CAUSE THE EXCEL FILE GOT CORRUPTED SO YOU'D HAVE TO RE DO
    #ALL OF THOSE ANYWAYS
    #IF YOU'RE STILL HAVING TROUBLE AFTER LUNCH TOMORROW, PLEASE GO TO TECH GUYS FOR HELP
    #breathe, you're doing really well, you still have a month, you deserve the breaks you take,
    #you're gonna be ok
    #screw that, you're going to be amazing

    """for i in range(20): #may need to change max value 
        if(x[peaks][i] <= 6800 and x[peaks][i] >= 6800):
            peakNum = peakNum + 1
            np.append(peak2, x[peaks][i])

    print("where are the peaks???")
    print(peak2)
    print("peakNum: " + str(peakNum))"""
        

    

    Halpha = integrate.simps(y2, x)
    NII = integrate.simps(y3, x)

    Hbeta = integrate.simps(y4, x)
    OIII = integrate.simps(y5, x)

    SII = integrate.simps(y6, x) + integrate.simps(y7, x)

    OI = integrate.simps(y8, x)
    
    f = plt.figure(3)
    plt.grid(c='k', linestyle='-', linewidth=1.0, alpha=0.25)
    
    plt.plot(x, y, c='k', linewidth=0.5)

    plt.plot(x, y2, c='r', linewidth=0.75)
    plt.plot(x, y3, c='g', linewidth=0.75)
    plt.plot(x, y4, c='c', linewidth=0.75)
    plt.plot(x, y5, c='y', linewidth=0.75)
    plt.plot(x, y6, c='b', linewidth=0.75)
    plt.plot(x, y7, c='r', linewidth=0.75)
    plt.plot(x, y8, c='m', linewidth=0.75)
    
    
    f.show()

    print("H-alpha: ", Halpha)
    print("NII: ", NII)
    print("H-beta: ", Hbeta)
    print("OIII: ", OIII)
    print("SII: ", SII)
    print("OI: ", OI)

    output.write(str(Halpha) + '\n')
    output.write(str(NII) + '\n')
    output.write(str(Hbeta) + '\n')
    output.write(str(OIII) + '\n')
    output.write(str(SII) + '\n')
    output.write(str(OI) + '\n')

    #was 1460 t0 1464
    noiseX = x[xi:xf]
    noiseY = y[xi:xf]
    noise = 3 * integrate.simps(noiseY, noiseX)
    print("noise: ", noise)
    output.write(str(noise) + '\n')

    elineFluxes = namedtuple('EmLines', 'hAlpha NII hBeta OIII SII OI noise')
    emLines = elineFluxes(Halpha, NII, Hbeta, OIII, SII, OI, noise)
    return emLines

def BPT(emLines, unclearOIII, unclearHA, unclearHB, unclearNII, unclearSII, unclearOI): 
    #emLines = gaussianFits()
    Halpha = emLines.hAlpha
    NII = emLines.NII
    Hbeta = emLines.hBeta
    OIII = emLines.OIII
    SII = emLines.SII
    OI = emLines.OI
    noise = emLines.noise
    flagY = 0 #-1 is down, 1 is up
    flagX = 0 #-1 is pointing to the left, 1 is pointing to the right

    if unclearOIII==1:
        OIII = noise
    if unclearHA==1:
        Halpha = noise
    if unclearHB==1:
        Hbeta = noise
    if unclearNII==1:
        NII = noise
    if unclearSII==1:
        SII = noise
    if unclearOI==1:
        OI = noise

    if unclearOIII==1:
        flagY = -1
    if unclearNII==1 or unclearSII==1 or unclearOI==1:
        flagX = -1
    if unclearHB==1:
        flagY = 1
    if unclearHA==1:
        flagX = 1
    #uncomment if logs are right, log(OIII/Hb) and log(NII/SII/OI / Ha)
    """if unclearOIII==1:
        flagY = -1
    if unclearNII==1 or unclearSII==1 or unclearOI==1: #yes
        flagX = -1
    if unclearHB==1:
        flagY = 1
    if unclearHA==1: #yes
        flagX = 1"""

    #NII BPT 
    g = plt.figure(4)

    #format and plot the delimeter equations
    xVals = np.arange(-3, 0.05, 0.01)
    xVals2 = np.arange(-3, 0.47, 0.01)
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

    #testing BPTs
    #plt.plot(xVals2, log_OIII_Hb_NII(xVals2), '-k')

    plt.plot(np.log(NII/Halpha), np.log(OIII/Hbeta), marker='o', c='k')
    if unclearNII==1 or unclearOIII==1 or unclearHA==1 or unclearHB==1: plt.quiver(np.log(NII/Halpha), np.log(OIII/Hbeta), flagX, flagY) #adds an arrow if uncertain
    g.show()

    #SII BPT
    h = plt.figure(5)

    #format and plot delimiter equations
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

    #plt.plot(xS, log_OIII_Hb_SII(xS), '-k')

    plt.plot(np.log(SII/Halpha), np.log(OIII/Hbeta), marker='o', c='k')
    if unclearSII==1 or unclearOIII==1 or unclearHA==1 or unclearHB==1: plt.quiver(np.log(SII/Halpha), np.log(OIII/Hbeta), flagX, flagY)
    h.show()

    #OI BPT 
    k = plt.figure(6)

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

    plt.plot(np.log(OI/Halpha), np.log(OIII/Hbeta), marker='o', c='k')
    if unclearOI==1 or unclearOIII==1 or unclearHA==1 or unclearHB==1: plt.quiver(np.log(OI/Halpha), np.log(OIII/Hbeta), flagX, flagY)
    k.show()
    
    raw_input()

def run_code():
    #saveFile = '/users/plira/india_temp/variable_spectra/starlight_model_graphs'
    reRun = raw_input("Enter yes if you would like to re-do the previous run with the same user input: ")
    previousRun = open('/users/plira/india_temp/previous_run.txt', 'r')
    #previousRun = open('/Users/plira/india_temp/variable_spectra/starlight_model_graphs/533803308112963584/saved_run_533803308112963584.txt', 'r') #enter file path for saved path and uncomment if necessary

    spObjId = raw_input("Please enter the SpObjId: ")
    outFile = '/users/plira/india_temp/variable_spectra/starlight_model_graphs/' + spObjId + '/output_' + spObjId + '.txt'
    output = open(outFile, 'w')

    oldInput = previousRun.readlines()
    print('reRun: ' + reRun)

    if reRun != 'yes': #if not doing a re run, open a file to write to, if you open this regardless, will overwrite
        newInput = open('/users/plira/india_temp/previous_run.txt', 'w')
    else:
        newInput = open('/users/plira/india_temp/dummy.txt', 'w')
        #output = open('/users/plira/india_temp/dummy_2.txt', 'w')
        #outFile = '/users/plira/india_temp/variable_spectra/starlight_model_graphs/' + spObjId + '/output_' + spObjId + '.txt'
        #output = open(outFile, 'w')

    emLines = gaussianFits(reRun, oldInput, newInput, output);

    Halpha = 0
    NII = 0
    Hbeta = 0
    OIII = 0
    SII = 0
    OI = 0
    if reRun=='yes':
        OIII = int(oldInput[12])
        output.write(oldInput[12])
        Halpha = int(oldInput[13])
        Hbeta = int(oldInput[14])
        NII = int(oldInput[15])
        SII = int(oldInput[16])
        OI = int(oldInput[17])

        output.write(oldInput[13]) #Ha
        output.write(oldInput[15]) #NII
        output.write(oldInput[14]) #Hb
        output.write(oldInput[12]) #OIII
        output.write(oldInput[16]) #SII
        output.write(oldInput[17]) #OI

    else:
        print("if any lines were uncertain, please answer yes when prompted...")
        unclearOIII = raw_input("OIII: ")
        unclearHA = raw_input("H-alpha: ")
        unclearHB = raw_input("H-beta: ")
        unclearNII = raw_input("NII: ")
        unclearSII = raw_input("SII: ")
        unclearOI = raw_input("OI: ")

        if unclearOIII == 'yes':
            OIII = 1
            newInput.write(str(OIII) + '\n')
            output.write(str(OIII) + '\n')
        else:
            newInput.write('0\n')
            output.write('0\n')

        if unclearHA == 'yes':
            Halpha = 1
            newInput.write(str(Halpha)+ '\n')
            output.write(str(Halpha)+ '\n')
        else:
            newInput.write('0\n')
            output.write('0\n')

        if unclearHB == 'yes':
            Hbeta = 1
            newInput.write(str(Hbeta)+ '\n')
            output.write(str(Hbeta)+ '\n')
        else:
            newInput.write('0\n')
            output.write('0\n')

        if unclearNII == 'yes':
            NII = 1
            newInput.write(str(NII) + '\n')
            output.write(str(NII) + '\n')
        else:
            newInput.write('0\n')
            output.write('0\n')

        if unclearSII == 'yes':
            SII = 1
            newInput.write(str(SII)+'\n')
            output.write(str(SII)+'\n')
        else:
            newInput.write('0\n')
            output.write('0\n')

        if unclearOI == 'yes':
            OI = 1
            newInput.write(str(OI)+ '\n')
            output.write(str(OI)+ '\n')
        else:
            newInput.write('0\n')
            output.write('0\n')
        newInput.close()

    BPT(emLines, OIII, Halpha, Hbeta, NII, SII, OI);

    output.close()
    saveFile = '/users/plira/india_temp/variable_spectra/starlight_model_graphs/' + spObjId + '/saved_run_' + spObjId + '.txt'
    print(saveFile)

    if reRun != 'yes':
        copyfile('/users/plira/india_temp/previous_run.txt', saveFile)

    shutil.copy2(outFile, '/users/plira/india_temp/variable_spectra/output_files/')

    previousRun.close()

if __name__ == '__main__':
    run_code()


