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


    def harmonicOscillator(self):
        return self.xMesh**2
        #plt.plot(self.xMesh,self.xMesh**2)
        #plt.show()


    def triangle(self):
        return np.abs(self.xMesh)
        #plt.plot(self.xMesh,np.abs(self.xMesh))
        #plt.show()
    

    def barrier(self, width):
        #enter width as float
        list = np.zeros([self.gridpoints+1])
        for i in range(self.gridpoints+1):
            if self.xMesh[i] < (self.xMin/width):
               list[i] = 0
            elif self.xMesh[i] > (self.xMax/width):
                list[i] = 0
            else:
                list[i] = self.A
        print list
        plt.plot(self.xMesh,list)
        plt.show()


    def teeth(self):
        list = np.zeros([self.gridpoints+1])
        for i in range(self.gridpoints+1):
            if i < self.gridpoints/4.0:
                list[i]=self.A
            elif i > 3*self.gridpoints/4.0:
                list[i]=self.A
            else:
                list[i]=0
        print list
        #plt.plot(self.xMesh,list)
        #plt.show()