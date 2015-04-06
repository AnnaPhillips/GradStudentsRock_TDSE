#! /usr/bin/python
import numpy as np
from scipy import sparse

class TDSE:

    #input: an array of size 1xN with the initial position, potential, delta x, delta t, "output
    #file name"
    def __init__(self, initWaveFunc, potential , delx, delt, timesteps, outputFile):
        self.N = len(initWaveFunc)
        self.V=potential
        self.timesteps=timesteps
        self.delx = delx
        self.delt = delt
        self.zeroMatrix = np.zeros([timesteps,self.N])
        #saving an N by N zero matrix whose first row is the real part of our
        #initWaveFunc array
        self.matrixReal = np.vstack([initWaveFunc.real,self.zeroMatrix])
        #saving an N by N zero matrix whose first row is the imaginary part of
        #our initWaveFunc array (if the array is all real, it will be all zeros)
        self.matrixImag = np.vstack([initWaveFunc.imag,self.zeroMatrix])
        self.outputFile = outputFile


    #Build NxN tridiagonal matrix to represent the lhs of equation.
    def getA(self):
        a = 2*np.ones(self.N)/self.delx**2-1j*np.ones(self.N)/self.delt
        a_sub = np.ones(self.N)
        data = [-1*a_sub,a,-1*a_sub]
        diags = [-1,0,1]
        A = sparse.spdiags(data,diags,self.N,self.N)
        print "original A"
        print A
        return A

    def getAwithVnonperiodic(self):
        A=np.zeros([self.N,self.N],dtype=complex)
        A[0,0]=-2/self.delx**2-1j +self.V[0]
        A[0,1]=1
        A[self.N-1,self.N-2]= 1
        A[self.N-1,self.N-1]=-2/self.delx**2-1j +self.V[self.N-1]
        for i in range(self.N-2):
            A[i+1,i]=-1
            A[i+1,i+1]=+2/self.delx**2-1j +self.V[i+1]
            A[i+1,i+2]=-1
        print "A with V non periodic"
        print A
        return A

    def getAwithVperiodic(self):
        A=np.zeros([self.N,self.N],dtype=complex)
        for i in range(self.N):
            A[i,(i-1) % self.N]=-1
            A[i,i % self.N]=2/self.delx**2-1j +self.V[i]
            A[i,(i+1) % self.N]=-1
        print "A with V periodic"
        print A
        return A
            

    def energy(self):
        normalization=1./(self.N)
        Elist=[]
        V=np.zeros(self.N)
        for time in range(self.timesteps):
            E=0
            for i in range(self.N):
                EtempUnnormed = (self.matrixReal[time,i]-1j*self.matrixImag[time,i]
                                 )*(-0.5*(self.matrixReal[time -1,i]+1j*self.matrixImag[time-1,i]-2*(self.matrixReal[time,i]
                                    +1j*self.matrixImag[time,i])+self.matrixReal[time+1,i]+1j*self.matrixImag[time+1,i])/ float(self.delx**2)
                                    + V[i]*(self.matrixReal[time, i]+1j*self.matrixImag[time,i]))
                E=E+EtempUnnormed*normalization
            Elist.append(E)
        return Elist


    def runNonPeriodicNoV(self):
        for n in range(self.timesteps):
            print("start of a time step")
            print(self.matrixReal)
            print(self.matrixImag)
            b = np.zeros(self.N,dtype=complex)
            Adense = self.getA().todense()
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
        hopefullythisworks=self.energy()
        print(hopefullythisworks)
        np.savetxt(self.outputFile + "_RealFiniteDifference.csv", self.matrixReal, delimiter=",")
        np.savetxt(self.outputFile + "_ImagFiniteDifference.csv", self.matrixImag, delimiter=",")


    def runNonPeriodicWithV(self):
        for n in range(self.timesteps):
            print("start of a time step")
            print(self.matrixReal)
            print(self.matrixImag)
            b = np.zeros(self.N,dtype=complex)
            A = self.getAwithVnonperiodic()
            for i in range(self.N):
                b[i] = -1.0/self.delt*(self.matrixImag[n,i]+1j*self.matrixReal[n,i])
            print("b=")
            print(b)
            solution = np.linalg.solve(A,b)
            print(solution)
            self.matrixReal[n+1,:] = solution.real
            self.matrixImag[n+1,:] = solution.imag
            print(self.matrixReal)
            print(self.matrixImag)
        hopefullythisworks=self.energy()
        print(hopefullythisworks)
        np.savetxt(self.outputFile + "_RealFiniteDifference.csv", self.matrixReal, delimiter=",")
        np.savetxt(self.outputFile + "_ImagFiniteDifference.csv", self.matrixImag, delimiter=",")

    def runPeriodicWithV(self):
        for n in range(self.timesteps):
            print("start of a time step")
            print(self.matrixReal)
            print(self.matrixImag)
            b = np.zeros(self.N,dtype=complex)
            A = self.getAwithVperiodic()
            for i in range(self.N):
                b[i] = -1.0/self.delt*(self.matrixImag[n,i]+1j*self.matrixReal[n,i])
            print("b=")
            print(b)
            solution = np.linalg.solve(A,b)
            print(solution)
            self.matrixReal[n+1,:] = solution.real
            self.matrixImag[n+1,:] = solution.imag
            print(self.matrixReal)
            print(self.matrixImag)
        hopefullythisworks=self.energy()
        print(hopefullythisworks)
        np.savetxt(self.outputFile + "_RealFiniteDifference.csv", self.matrixReal, delimiter=",")
        np.savetxt(self.outputFile + "_ImagFiniteDifference.csv", self.matrixImag, delimiter=",")
