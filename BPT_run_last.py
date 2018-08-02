#!/usr/bin/python 

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

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

graph_data = open('/Users/plira/india_temp/variable_spectra/auto/cumulative_BPT.txt', 'r').read()
#lines = graph_data.strip()
lines = graph_data.split('\n')
xNIIc = []
yNIIc = []
xSIIc = []
ySIIc = []
xOIc = []
yOIc = []

NIIcolors = []
SIIcolors = []
OIcolors = []

for line in lines:
    if len(line) > 1:
        #'Filename','hAlpha','NII','hBeta','OIII','SII','OI','noise','hAcl','NIIisClear','hBisClear','OIIIisClear','SIIisClear','OIisClear','tooNoisy', 'BPT NII', 'BPT SII', 'BPT OI'
        #name, hA, NII, hB, OIII, SII, OI, noise, hAc, NIIc, hBc, OIIIc, SIIc, OIc, tooNoisy, BPTNII, BPTSII, BPTOI = line.split(',')
        output = line.split(',')
        #print(output)
        hA = output[1]
        NII = output[2]
        hB = output[3]
        OIII = output[4]
        SII = output[5]
        OI = output[6]
        noise = output[7]
        hAcl = output[8] #1 for clear, 0 for unclear
        NIIcl = output[9]
        hBcl = output[10]
        OIIIcl = output[11]
        SIIcl = output[12]
        OIcl = output[13]
        tooNoisy = output[14]
        #BPTNII = output[15]
        #BPTSII = output[16]
        #BPTOI = output[17]
        
        xNII = np.log(float(NII)/float(hA))
        yNII = np.log(float(OIII)/float(hB))
        xNIIc.append(xNII)
        yNIIc.append(yNII)
        #print(xNIIc)
        
        xSII = np.log(float(SII)/float(hA))
        ySII = np.log(float(OIII)/float(hB))
        xSIIc.append(xSII)
        ySIIc.append(ySII)
        #print(xSIIc)
        
        xOI = np.log(float(OI)/float(hA))
        yOI = np.log(float(OIII)/float(hB))
        xOIc.append(xOI)
        yOIc.append(yOI)
        #print(xOIc)

        NIIcolors.append(NIIcl, hAcl, OIIIcl, hBcl)
        SIIcolors.append(SIIcl, hAcl, OIIIcl, hBcl)
        OIcolors.append(OIcl, hAcl, OIIIcl, hBcl)

    flagY = 0 #-1 is down, 1 is up
    flagX = 0 #-1 is pointing to the left, 1 is pointing to the right
    allClear = True
    allClearNII = True
    allClearSII = True
    allClearOI = True
    
    if hAcl==False and Halpha < noise: #right now there is a problem if both H-alpha and NII/SII/OI are both uncertain, same for OIII and HB
        Halpha = noise   #but if H-alpha is uncertain, the spectrum is too noisy anyways 
        flagX = 1
        print("hA is uncertain")
        hAClear = False
        allClear = False
    if NIIcl==False and NII < noise:
        NII = noise
        flagX = -1
        print("NII is uncertain")
        NIIClear = False
        allClearNII = False
    if hBcl==False and Hbeta < noise:
        Hbeta = noise
        flagY = 1
        print("hB is uncertain")
        hBClear = False
        allClear = False
    if OIIIcl==False and OIII < noise:
        OIII = noise
        flagY = -1
        print("OIII is uncertain")
        OIIIClear = False
        allClear = False
    if SIIcl==False and SII < noise:
        SII = noise
        flagX = -1
        print("SII is uncertain")
        SIIClear = False
        allClearSII = False
    if OIcl==False or OI < noise: #this one is different because OI is very small and often unclear
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
    elif np.log(OI/Halpha) >= -0.59 and np.log(OIII/Hbeta) < delim2SII:
        colorVarOI = '#ffff00'
        bptOI = 'LINER'
    elif np.log(OIII/Hbeta) < delim1OI: #check this
        colorVarOI = 'r'
        bptOI = 'HII'

    ax3.plot(np.log(OI/Halpha), np.log(OIII/Hbeta), marker='o', c=colorVarOI, markeredgecolor=colorVarOI)
    if allClear==False or allClearOI==False: ax3.quiver(np.log(OI/Halpha), np.log(OIII/Hbeta), flagX, flagY, color=colorVarOI)
    plt.tight_layout()
    c.savefig(outLoc + 'BPT Diagrams', bbox_inches='tight')
    else:
        continue
    
for j in range(len(xNIIc)):
    ax1.plot(xNIIc[j], yNIIc[j], marker='o', c='k')
    ax2.plot(xSIIc[j], ySIIc[j], marker='o', c='k')
    ax3.plot(xOIc[j], yOIc[j], marker='o', c='k')


    

plt.tight_layout()
plt.show()
    
