import os
import shutil
import numpy as np
import TDSE
import potentials
import create_wave
import sys

def makedir(name):
    if os.path.isdir(name):
        shutil.rmtree(name)
    os.mkdir(name)
    

# change directory to where the data will be stored
makedir('Runs')
os.chdir('Runs')

# set up testing conditions
initWaveFunc = create_wave.gaussian(0,0,1)
timesteps = 1000
xMin=-20.0
xMax=20.0
gridpoints=400
delx=(xMax-xMin)/gridpoints
delt=delx/1.
outputFilePeriodic = "periodic"
outputFileNonPeriodic = "nonPeriodic"


#----------------------------------------
# Free Particle
#----------------------------------------
makedir('freeParticle')
os.chdir('freeParticle')
amplitude = 1
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.freeParticle()
#set up and run free particle non-periodic conditions (simple finite difference method)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, True, False, outputFilePeriodic)
funstuff.run() # run finite difference scheme
#set up and run free particle periodic conditions (simple finite difference method)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFileNonPeriodic)
funstuff.run()
#set up and run free particle non-periodic conditions (Crank-Nicolson scheme)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, True, True, outputFilePeriodic + '_CN')
funstuff.run() # run finite difference scheme
#set up and run free particle periodic conditions (Crank-Nicolson scheme)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, True, outputFileNonPeriodic + '_CN')
funstuff.run()
os.chdir('..')

#i'll make this work properly later
sys.exit(0)

#----------------------------------------
# Square Well
#----------------------------------------
makedir('squareWell')
os.chdir('squareWell')
amplitude = 1

---------------------------------------
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.freeParticle()
#set up and run free particle non-periodic conditions (simple finite difference method)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, True, False, outputFilePeriodic)
funstuff.run() # run finite difference scheme
#set up and run free particle periodic conditions (simple finite difference method)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFileNonPeriodic)
funstuff.run()
#set up and run free particle non-periodic conditions (Crank-Nicolson scheme)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, True, True, outputFilePeriodic + '_CN')
funstuff.run() # run finite difference scheme
#set up and run free particle periodic conditions (Crank-Nicolson scheme)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, True, outputFileNonPeriodic + '_CN')
funstuff.run()
os.chdir('..')




#----------------------------------------
# Harmonic Oscillator
#----------------------------------------
makedir('harmonicOscillator')
os.chdir('harmonicOscillator')
amplitude = 1
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.harmonicOscillator()
outputFile = "harmonicOscillator"
#set up and run harmonic oscillator with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme
os.chdir('..')



#----------------------------------------
# Triangle
#----------------------------------------
makedir('triangle')
os.chdir('triangle')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.triangle()
outputFile = "triangle"
#set up and run triangle with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme
os.chdir('..')



#----------------------------------------
# Kronig-Penney
#----------------------------------------





#----------------------------------------
# Barrier
#----------------------------------------
makedir('barrier')
os.chdir('barrier')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.barrier(20.0) #can set width of barrier
outputFile = "nonPeriodic"
#set up and run barrier with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme
os.chdir('..')


#----------------------------------------
# Teeth
#----------------------------------------
makedir('teeth')
os.chdir('teeth')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.teeth()
outputFile = "periodic"
#set up and run barrier with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme
os.chdir('..')



#----------------------------------------
# V=ix
#----------------------------------------




#----------------------------------------
# v=x+ix
#----------------------------------------





os.chdir('..')