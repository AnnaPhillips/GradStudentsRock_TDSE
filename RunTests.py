import os
import numpy as np
import TDSE
import potentials
import create_wave

# change directory to where the data will be stored
os.mkdir('Runs')
os.chdir('Runs')

# set up testing conditions
initWaveFunc = create_wave.gaussian(0,0,1)
timesteps = 10
xMin=-10.0
xMax=10.0
gridpoints=200
delx=1
delt=1

#----------------------------------------
# Free Particle
#----------------------------------------
os.mkdir('freeParticle')
os.chdir('freeParticle')
potential = 0
outputFilePeriodic = "free_periodic"
outputFileNonPeriodic = "free_nonPeriodic"
#set up and run free particle non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, False, outputFilePeriodic)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme
#set up and run free particle periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, timesteps, True, outputFileNonPeriodic)
funstuff.run()
#funstuff.runOS()
os.chdir('..')



#----------------------------------------
# Square Well
#----------------------------------------




#----------------------------------------
# Harmonic Oscillator
#----------------------------------------
os.mkdir('harmonicOscillator')
os.chdir('harmonicOscillator')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.harmonicOscillator()
outputFile = "harmonicOscillator"
#set up and run harmonic oscillator with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, timesteps, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme




#----------------------------------------
# Triangle
#----------------------------------------
os.mkdir('Triangle')
os.chdir('Triangle')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.triangle()
outputFile = "triangle"
#set up and run triangle with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, timesteps, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme




#----------------------------------------
# Kronig-Penney
#----------------------------------------





#----------------------------------------
# Barrier
#----------------------------------------
os.mkdir('Triangle')
os.chdir('Triangle')
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)
potential = runPotential.barrier(20.0) #can set width of barrier
outputFile = "barrier"
#set up and run barrier with non-periodic conditions
funstuff=TDSE.TDSE(initWaveFunc, potential, timesteps, False, outputFile)
funstuff.run() # run finite difference scheme
#funstuff.runOS() # run other scheme



#----------------------------------------
# V=ix
#----------------------------------------




#----------------------------------------
# v=x+ix
#----------------------------------------





os.chdir('..')