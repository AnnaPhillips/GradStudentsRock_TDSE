#! /usr/bin/python
import numpy as np
import matplotlib.pyplot as plt


class Potentials:
    def __init__(self, xMin, xMax, gridpoints, amplitude):
        self.xMin = xMin
        self.xMax = xMax
        self.gridpoints = gridpoints
        self.A = amplitude
        self.xMesh = np.array([self.xMin + i*(self.xMax-self.xMin)/float(self.gridpoints) for i in range(gridpoints+1)])


    def freeParticle(self):
        return self.xMesh*0
        #plt.plot(self.xMesh,self.xMesh**2)
        #plt.show()


    def harmonicOscillator(self):
        return self.xMesh**2
        #plt.plot(self.xMesh,self.xMesh**2)
        #plt.show()


    def triangle(self):
        return np.abs(self.xMesh)
        #plt.plot(self.xMesh,np.abs(self.xMesh))
        #plt.show()


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
                
                
    def imagPotential(self):
        # V=ix
        return 1.j*self.xMesh
        
        
    def complexPotential(self):
        # V=x+ix
        return self.xMesh + 1.j*self.xMesh
    
