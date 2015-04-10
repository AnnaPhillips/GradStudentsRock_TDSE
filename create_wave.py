#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

#This file creates a real gaussian wave packet and defines the xMesh function that potentials.py calls on. 

def xMesh(xMin, xMax, gridpoints):
    return np.array([xMin + i*(xMax-xMin)/gridpoints for i in range(gridpoints+1)])


def gaussian(x0, k0, var):
    initialPos = (var**2/(2*np.pi))**0.25*np.exp(0.25*(xAR - x0)*(4.j*k0 + (-xAR + x0)*var**2))
    return initialPos

xAR = xMesh(-10., 10., 200) #make sure the bounds are float types
wave = gaussian(0, 0, 1)

#print xAR
#print wave

#plt.scatter(xAR, wave)
#plt.show()
