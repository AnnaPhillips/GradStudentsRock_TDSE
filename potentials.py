#! /usr/bin/python
import numpy as np
import matplotlib.pyplot as plt


#This file creates arrays for each type of potential we might want to use. We call on these in our RunTests.py file and the results get pushed through the calculations done by TDSE.py. 

class Potentials:
    #The inputs for this class of objects are the minimum and maximum spacial points, the number of spacial grid points and an overall amplitude factor for the potential. 
    def __init__(self, xMin, xMax, gridpoints, amplitude):
        self.xMin = xMin
        self.xMax = xMax
        self.gridpoints = gridpoints
        self.A = amplitude
        self.xMesh = np.array([self.xMin + i*(self.xMax-self.xMin)/float(self.gridpoints) for i in range(gridpoints+1)])


    #The free particle potential simply creates an array of zeros that is the length of our grid points. 
    def freeParticle(self):
        return self.xMesh*0
        #plt.plot(self.xMesh,self.xMesh**2)
        #plt.show()

    #The Harmonic Oscillator potential creates a grid of values A*x**2 for all of our grid points. Note that this means the center of our oscillator is fixed at x=0. This should be used with non-periodic boundary conditions. 
    def harmonicOscillator(self):
        return self.A*self.xMesh**2
        #plt.plot(self.xMesh,self.xMesh**2)
        #plt.show()

    #The triangle potential creates a grid of values A*abs(x) for all of our grid points. Note that this means the center of our triangle is fixed at x=0. This should be used with non-periodic boundary conditions.
    def triangle(self):
        return self.A*np.abs(self.xMesh)
        #plt.plot(self.xMesh,np.abs(self.xMesh))
        #plt.show()

    #This creates a barrier with a defined width centered on X=0. This should be used with non-periodic boundary conditions. 
    def barrier(self, barWidth):
        #enter width as float
        list = np.zeros([self.gridpoints+1])
        for i in range(self.gridpoints+1):
            if self.xMesh[i] > (self.xMin/(barWidth/2)) and self.xMesh[i] < (self.xMax/(barWidth/2)):
                list[i] = self.A
        return list
        #print list
        #plt.plot(self.xMesh,list)
        #plt.show()

    #This creates a semi-infinite square well of a defined barrier width centered on X=0. This should be used with our non-periodic boundary conditions.
    def squareWell(self, barWidth):
        #enter width as float
        list = np.zeros([self.gridpoints+1])
        for i in range(self.gridpoints+1):
            if self.xMesh[i] < (self.xMin/(barWidth/2)) or self.xMesh[i] > (self.xMax/(barWidth/2)):
                list[i] = self.A*1000.
        return list
        #print list
        #plt.plot(self.xMesh,list)
        #plt.show()

    #This defines our "saw tooth" potential. It should be used with the periodic boundary conditions. 
    def kronigPenney(self, gapWidth, numberWells):
        list = np.zeros([self.gridpoints+1])
        tooth=(self.xMax - self.xMin)/numberWells
        for i in range(self.gridpoints+1):
            if (self.xMesh[i] - self.xMin)%tooth > gapWidth:
                list[i]=self.A
        return list
        #print list
        #plt.plot(self.xMesh,list)
        #plt.show()
                
    #This defines an imaginary potential. It should be used with the non-periodic boundary conditions.            
    def imagPotential(self):
        # V=ix
        return self.A*1.j*self.xMesh
        
    #This defines another potential. It should be used with the non-periodic boundary conditions.            
    def complexPotential(self):
        # V=x+ix
        return self.xMesh + 1.j*self.xMesh
    
