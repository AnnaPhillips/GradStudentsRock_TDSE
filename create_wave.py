#! /usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt

#This file create various initialwave functions to be pushed through TDSE.py and defines the xMesh function that potentials.py calls on. 

#xMesh is array of points (xmin, xmin+L/gridpoints, ...., xmax) that is gridpoints+1 long to capture both beginning and end points. 
def xMesh(xMin, xMax, gridpoints):
    return np.array([xMin + i*(xMax-xMin)/gridpoints for i in range(gridpoints+1)])


#This defines a complex gaussian wavepacket. x0 is the initial position, k0 is initial wave number, and var is the width. xAR is the array spat out of xMesh
def gaussian(x0, k0, var, xAR):
    initialPos = (var**2/(2*np.pi))**0.25*np.exp(0.25*(xAR - x0)*(4.j*k0 + (-xAR + x0)*var**2))
    return initialPos

#This is a real cosine wave. It takes in a harmonic number n, xAR, and the width of the well. 
def cosine(n, xAR, barWidth):
    initialPos=np.cos(xAR*np.pi*n/barWidth)*(abs(xAR) < barWidth/2)
    return initialPos


def expo(amp, xAR):
    initialPos=np.exp(-1*math.sqrt(amp/2)*xAR**2)
    return initialPos


#Testing code is commented out. 
#xAR = xMesh(-10., 10., 200) #make sure the bounds are float types
#wave = gaussian(0, 0, 1)

#print xAR
#print wave

#plt.scatter(xAR, wave)
#plt.show()
