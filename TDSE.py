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
        a = -2*np.ones(self.N)/self.delx**2-1j*np.ones(self.N)/self.delt
        a_sub = np.ones(self.N)
        data = [a_sub,a,a_sub]
        diags = [-1,0,1]
        A = sparse.spdiags(data,diags,self.N,self.N)
        print A
        return A



    def run(self):
        for n in range(self.N):
            print("start of a time step")
            print(self.matrixReal)
            print(self.matrixImag)
            b = np.zeros(self.N,dtype=complex)
            Adense = self.getA().todense()
            print(Adense)
            for i in range(self.N):
                b[i] = -1.0/self.delt*(self.matrixImag[n,i]+1j*self.matrixReal[n,i])
            print("b=")
            print(b)
            solution = np.linalg.solve(Adense,b)
            print(solution)
            self.matrixReal[n+1,:] = solution.real
            self.matrixImag[n+1,:] = solution.imag
            print(self.matrixReal)
            print(self.matrixImag)
            

        np.savetxt(self.outputFile + "_Real", self.matrixReal, delimiter=",")
        np.savetxt(self.outputFile + "_Imag", self.matrixImag, delimiter=",")