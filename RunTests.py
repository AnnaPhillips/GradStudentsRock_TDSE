import os
import numpy as np
import TDSE
import potentials
import create_wave
import sys

#This file is what should be run to run the code. It runs our different cases one by one. 

# create function to make a directory only if it does not already exist
def makedir(name):
    if not os.path.isdir(name):
        os.mkdir(name)
    

# change directory to where the data will be stored
makedir('Runs')
os.chdir('Runs')

# set up testing conditions
timesteps = 800
xMin=-20.0
xMax=20.0
gridpoints=800
xAR = create_wave.xMesh(xMin, xMax, gridpoints)
delx=(xMax-xMin)/gridpoints
delt=delx/4.
outputFilePeriodic = "periodic"
outputFileNonPeriodic = "nonPeriodic"

# set up our potentials generator
runPotential = potentials.Potentials(xMin, xMax, gridpoints,xAR)



#----------------------------------------
# Free Particle
#----------------------------------------
caseName='freeParticle'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=-5.0
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)

# pull the free particle potential from the potentials generator
potential = runPotential.freeParticle()

# run free particle non-periodic conditions (simple finite difference method)
nonPeriodicSFD = TDSE.TDSE(initWaveFunc, potential, delx, delt,timesteps,
                           False, False, outputFileNonPeriodic + '_SFD')
nonPeriodicSFD.run()

# run free particle periodic conditions (simple finite difference method)
periodicSFD = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                        True, False, outputFilePeriodic + '_SFD')
periodicSFD.run()

# run free particle non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

# run free particle periodic conditions (Crank-Nicolson scheme)
periodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                       True, True, outputFilePeriodic + '_CN')
periodicCN.run()

os.chdir('..')

# just for now until we want to run all these
#sys.exit(0)



#----------------------------------------
# Square Well with initial wave packet
#----------------------------------------
caseName='squareWell'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=0.1
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)
# pull the potential from the potentials generator
potential = runPotential.squareWell(20.,1000.)

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')


#----------------------------------------
# Square Well with an eigenstate
#----------------------------------------
caseName='squareWellwithEigenState'
makedir(caseName)
os.chdir(caseName)
print(caseName)

harmonic=1. #Pick higher energy eigenstates if you want!
barWidth=20. #Our well will be 20 across. It needs to be a float. 

initWaveFunc = create_wave.cosine(harmonic,xAR,barWidth)
# pull the potential from the potentials generator
potential = runPotential.squareWell(barWidth,1000.)
#Basically, this is a semi-fininite well. The barriers are 1000 units high. Our initial wave functions all have intitial energies on the order of 1 to 10. 

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')


#----------------------------------------
# Harmonic Oscillator with wave packet. 
#----------------------------------------
caseName='harmonicOscillator'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=0.1
waveVar=1.0
waveK0=2.0
amp=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)
# pull the initial wave function from create_wave

# pull the potential from the potentials generator
potential = runPotential.harmonicOscillator(amp)

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')




#----------------------------------------
# Harmonic Oscillator with lowest energy state
#----------------------------------------
caseName='harmonicOscillatorEigenState'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=0.1
waveVar=1.0
waveK0=2.0
amp=2.0
initWaveFunc = create_wave.expo(amp, xAR)
# pull the initial wave function from create_wave

# pull the potential from the potentials generator
potential = runPotential.harmonicOscillator(amp)

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')



#----------------------------------------
# Triangle
#----------------------------------------
caseName='triangle'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=0.1
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)

# pull the potential from the potentials generator
potential = runPotential.triangle(2.)

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')



#----------------------------------------
# Kronig-Penney
#----------------------------------------
caseName='kronigPenney'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=-5.0
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)

# pull the potential from the potentials generator
potential = runPotential.kronigPenney(7.,5.,4.)

# run periodic conditions (Crank-Nicolson scheme)
PeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          True, True, outputFilePeriodic + '_CN')
PeriodicCN.run()

os.chdir('..')



#----------------------------------------
# Imaginary Potential (V=ix)
#----------------------------------------
caseName='imagPotential'
makedir(caseName)
os.chdir(caseName)
print(caseName)


waveX0=5.0
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)

# pull the potential from the potentials generator
potential = runPotential.imagPotential(5.)

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')



#----------------------------------------
# Complex Potential (V=x+ix)
#----------------------------------------
caseName='complexPotential'
makedir(caseName)
os.chdir(caseName)
print(caseName)

waveX0=-5.0
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)

# pull the potential from the potentials generator
potential = runPotential.complexPotential(4.)

# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()

os.chdir('..')



#----------------------------------------
#Barriers
#----------------------------------------

#We need to have a larger domain for the barriers

timesteps = 800
xMin=-60.0
xMax=60.0
gridpoints=800
xAR = create_wave.xMesh(xMin, xMax, gridpoints)
delx=(xMax-xMin)/gridpoints
delt=delx/4.

# barrier height=E

waveX0=-15.0
waveVar=1.0
waveK0=2.0
initWaveFunc = create_wave.gaussian(waveX0,waveK0,waveVar,xAR)
caseName='barrier1'
amplitude = 4.
makedir(caseName)
os.chdir(caseName)
print(caseName)

# pull the potential from the potentials generator
potential = runPotential.barrier(1.0,amplitude) #can set width of barrier
# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()
os.chdir('..')


# barrier height<E
caseName='barrier2'
amplitude = amplitude*.5
makedir(caseName)
os.chdir(caseName)
print(caseName)

# pull the potential from the potentials generator
potential = runPotential.barrier(1.0, amplitude) #can set width of barrier
# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()
os.chdir('..')


# barrier height>E
caseName='barrier3'
amplitude = amplitude*3
makedir(caseName)
os.chdir(caseName)
print(caseName)
# pull the potential from the potentials generator
potential = runPotential.barrier(1.0,amplitude) #can set width of barrier
# run non-periodic conditions (Crank-Nicolson scheme)
nonPeriodicCN = TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps,
                          False, True, outputFileNonPeriodic + '_CN')
nonPeriodicCN.run()
os.chdir('..')


# back out to main folder
os.chdir('..')
