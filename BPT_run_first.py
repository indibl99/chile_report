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

def animate(i):
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
            #'Filename','hAlpha','NII','hBeta','OIII','SII','OI','noise','hAisClear','NIIisClear','hBisClear','OIIIisClear','SIIisClear','OIisClear','tooNoisy', 'BPT NII', 'BPT SII', 'BPT OI'
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
            hAc = output[8] #1 for clear, 0 for unclear
            NIIc = output[9]
            hBc = output[10]
            OIIIc = output[11]
            SIIc = output[12]
            OIc = output[13]
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
        else:
            continue
        
    for j in range(len(xNIIc)):
        ax1.plot(xNIIc[j], yNIIc[j], marker='o', c='k')
        ax2.plot(xSIIc[j], ySIIc[j], marker='o', c='k')
        ax3.plot(xOIc[j], yOIc[j], marker='o', c='k')

ani = animation.FuncAnimation(c, animate, interval=1000)
plt.tight_layout()
plt.show()
    
