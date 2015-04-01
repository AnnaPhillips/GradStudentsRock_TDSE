#! /usr/bin/python
import numpy as np
from scipy import sparse

class TDSE:

    #input: an array of size 1xN with the initial position, delta x, delta t, "output
    #file name"
    def __init__(self, initialPos, delx, delt, outputFile):
        self.N = len(initialPos)
        self.delx = delx
        self.delt = delt
        self.zeroMatrix = np.zeros([self.N,self.N])
        #saving an N by N zero matrix whose first row is the real part of our
        #initialPos array
        self.matrixReal = np.vstack([initialPos.real,self.zeroMatrix])
        #saving an N by N zero matrix whose first row is the imaginary part of
        #our initialPos array (if the array is all real, it will be all zeros)
        self.matrixImag = np.vstack([initialPos.imag,self.zeroMatrix])
        self.outputFile = outputFile


    #Build NxN tridiagonal matrix to represent the lhs of equation.
    def getA(self):
        a = -2*np.ones(self.N)/self.delx**2-complex(0,1)*np.ones(self.N)/self.delt
        a_sub = np.ones(self.N)
        data = [a_sub,a,a_sub]
        diags = [-1,0,1]
        A = sparse.spdiags(data,diags,self.N,self.N)
        return A


    def run(self):
        for n in range(self.N):
            b = np.zeros(self.N)
            Adense = self.getA().todense()
            for i in range(self.N):
                #This is where it seems to have trouble:
                b[i-1] = -1.0/self.delt*complex(self.matrixImag[n-1,i-1],
                                           self.matrixReal[n-1,i-1])
            solution = np.linalg.solve(Adense,b)
            self.matrixReal[n,:] = solution.real
            self.matrixImag[n,:] = solution.imag

        np.savetxt(self.outputFile + "_Real", self.matrixReal, delimiter=",")
        np.savetxt(self.outputFile + "_Imag", self.matrixImag, delimiter=",")





