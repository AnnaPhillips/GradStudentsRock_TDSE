#! /usr/bin/python
import numpy as np
from scipy import sparse

class TDSE:

    #input: an array of size 1xN with the initial position, potential, delta x, delta t, time steps, periodic, "output
    #file name"
    def __init__(self, initWaveFunc, potential, delx, delt, timesteps, periodic, useCN, outputFile):
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
        self.periodic=periodic
        self.useCN=useCN


    def getAnonPeriodic(self):
        A = np.zeros([self.N,self.N],dtype=complex)
        bPrime = np.zeros([self.N,self.N],dtype=complex)

        A[0,0] = 2/self.delx**2 - 1j/self.delt + self.V[0]
        A[0,1] = -1/self.delx**2
        A[self.N-1,self.N-2] = -1/self.delx**2
        A[self.N-1,self.N-1] = 2/self.delx**2 - 1j/self.delt + self.V[self.N-1]

        bPrime[0,0] = -1j/self.delt
        bPrime[self.N-1,self.N-1] = -1j/self.delt

        for i in range(self.N-2):
            A[i+1,i] = -1/self.delx**2
            A[i+1,i+1] = 2/self.delx**2 - 1j/self.delt + self.V[i+1]
            A[i+1,i+2] = -1/self.delx**2

            bPrime[i+1,i+1] = -1j/self.delt
            
        print "A is non periodic"
        #print A
        #print bPrime
        return A, bPrime


    def getAperiodic(self):
        A = np.zeros([self.N,self.N],dtype=complex)
        bPrime = np.zeros([self.N,self.N],dtype=complex)
        
        for i in range(self.N):
            #negative array references count backwards, so we only need to mod the last one of these
            A[i,(i-1)] = -1/self.delx**2
            A[i,i] = 2/self.delx**2 - 1j/self.delt + self.V[i]
            A[i,(i+1) % self.N] = -1/self.delx**2
           
            bPrime[i,i] = -1j/self.delt
            
        print "A is periodic"
        #print A
        #print bPrime
        return A, bPrime


    def getAnonPeriodicCN(self):#CN is for Crank-Nicolson
        A = np.zeros([self.N,self.N],dtype=complex)
        bPrime = np.zeros([self.N,self.N],dtype=complex)

        A[0,0] =  1 + 1j*self.delt/2*(2/self.delx**2 + self.V[0])
        A[0,1] = -1j*self.delt/(2*self.delx**2)
        A[self.N-1,self.N-2] = -1j*self.delt/(2*self.delx**2)
        A[self.N-1,self.N-1] =  1 + 1j*self.delt/2*(2/self.delx**2 + self.V[self.N-1])

        bPrime[0,0] =  1 - 1j*self.delt/2*(2/self.delx**2 + self.V[0])
        bPrime[0,1] = 1j*self.delt/(2*self.delx**2)
        bPrime[self.N-1,self.N-2] = 1j*self.delt/(2*self.delx**2)
        bPrime[self.N-1,self.N-1] =  1 - 1j*self.delt/2*(2/self.delx**2 + self.V[self.N-1])

        for i in range(self.N-2):
            A[i+1,i] = -1j*self.delt/(2*self.delx**2)
            A[i+1,i+1] =  1 + 1j*self.delt/2*(2/self.delx**2 + self.V[i+1])
            A[i+1,i+2] = -1j*self.delt/(2*self.delx**2) 
            
            bPrime[i+1,i] = 1j*self.delt/(2*self.delx**2)
            bPrime[i+1,i+1] = 1 - 1j*self.delt/2*(2/self.delx**2 + self.V[i+1])
            bPrime[i+1,(i+2) % self.N] = 1j*self.delt/(2*self.delx**2)

        print "A is non periodic"
        #print A
        #print bPrime
        return A, bPrime


    def getAperiodicCN(self):#CN is for Crank-Nicolson
        A = np.zeros([self.N,self.N],dtype=complex)
        bPrime = np.zeros([self.N,self.N],dtype=complex)
        
        for i in range(self.N):
            #negative array references count backwards, so we only need to mod the last one of these
            A[i,(i-1)] = -1j*self.delt/(2*self.delx**2)
            A[i,i] = 1 + 1j*self.delt/2*(2/self.delx**2 + self.V[i])
            A[i,(i+1) % self.N] = -1j*self.delt/(2*self.delx**2)

            bPrime[i,(i-1)] = 1j*self.delt/(2*self.delx**2)
            bPrime[i,i] = 1 - 1j*self.delt/2*(2/self.delx**2 + self.V[i])
            bPrime[i,(i+1) % self.N] = 1j*self.delt/(2*self.delx**2)
            
        print "A is periodic"
        #print A
        #print bPrime
        return A, bPrime


    def getEigenStates(self):
        if self.periodic:
            A=self.getAperiodic()
        else:
            A=self.getAnonPeriodic()
        Areal=np.real(A)
        Energies,Evectors=np.linalg.eig(Areal)
        print "energies"
        print Energies
        print "eigenvector"
        print Evectors
        #print A

    def energy(self):
        normalization=1./(self.N)
        Elist=[]
        V=np.zeros(self.N) #we need to change this to reflect the potentials
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


    def run(self):
        if self.useCN:
            if self.periodic:
                A, bPrime = self.getAperiodicCN()
            else:
                A, bPrime = self.getAnonPeriodicCN()
        else:
            if self.periodic:
                A, bPrime = self.getAperiodic()
            else:
                A, bPrime = self.getAnonPeriodic()
        self.getEigenStates()
        for n in range(self.timesteps):
            print("start of a time step")
            print(self.matrixReal)
            print(self.matrixImag)
            b = np.zeros(self.N,dtype=complex)
            for k in range(self.N):
                b[k] = bPrime[n,k]*(self.matrixReal[n,k] + 1j*self.matrixImag[n,k])
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

