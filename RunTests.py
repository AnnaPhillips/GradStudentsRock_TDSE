import os
import shutil
import numpy as np
import TDSE
import potentials
import create_wave

def makedir(name):
    if os.path.isdir(name):
        shutil.rmtree(name)
    os.mkdir(name)
    

# change directory to where the data will be stored
makedir('Runs')
os.chdir('Runs')

# set up testing conditions
initWaveFunc = create_wave.gaussian(0,0,1)
timesteps = 10
xMin=-10.0
xMax=10.0
gridpoints=200
delx=(xMax-xMin)/gridpoints
delt=delx/10


#----------------------------------------
# Free Particle
#----------------------------------------
makedir('freeParticle')
os.chdir('freeParticle')
amplitude = 1
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.freeParticle()
outputFilePeriodic = "free_periodic"
outputFileNonPeriodic = "free_nonPeriodic"
#set up and run free particle non-periodic conditions (simple finite difference method)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, True, False, outputFilePeriodic)
funstuff.run() # run finite difference scheme
#set up and run free particle periodic conditions (simple finite difference method)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, False, outputFileNonPeriodic)
funstuff.run()
#set up and run free particle non-periodic conditions (Crank-Nicolson scheme)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, True, True, outputFilePeriodic + 'CN')
funstuff.run() # run finite difference scheme
#set up and run free particle periodic conditions (Crank-Nicolson scheme)
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, True, outputFileNonPeriodic + 'CN')
funstuff.run()
#funstuff.runOS()
os.chdir('..')


#----------------------------------------
# Square Well
#----------------------------------------




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
makedir('Triangle')
os.chdir('Triangle')
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
makedir('Barrier')
os.chdir('Barrier')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.barrier(20.0) #can set width of barrier
outputFile = "barrier"
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