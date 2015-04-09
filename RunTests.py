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
# set the amplitude for the potential
amplitude = 1

# set up our potentials generator
runPotential = potentials.Potentials(xMin, xMax, gridpoints, amplitude)

#----------------------------------------
# Free Particle
#----------------------------------------
makedir('freeParticle')
os.chdir('freeParticle')

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

#i'll make this work properly later
sys.exit(0)

#----------------------------------------
# Square Well
#----------------------------------------
makedir('squareWell')
os.chdir('squareWell')


#---------------------------------------
os.chdir('..')



#----------------------------------------
# Harmonic Oscillator
#----------------------------------------
makedir('harmonicOscillator')
os.chdir('harmonicOscillator')


os.chdir('..')


#----------------------------------------
# Triangle
#----------------------------------------
makedir('triangle')
os.chdir('triangle')


os.chdir('..')


#----------------------------------------
# Kronig-Penney
#----------------------------------------
makedir('kronigPenney')
os.chdir('kronigPenney')


os.chdir('..')



#----------------------------------------
# Teeth
#----------------------------------------
makedir('teeth')
os.chdir('teeth')


os.chdir('..')



#----------------------------------------
# V=ix
#----------------------------------------
makedir('imag1')
os.chdir('imag1')


os.chdir('..')



#----------------------------------------
# v=x+ix
#----------------------------------------
makedir('imag2')
os.chdir('imag2')


os.chdir('..')


#----------------------------------------
# Barrier * I moved this last because we want to set different amplitudes and I had a global amplitude set at the beginning
#----------------------------------------

# barrier height=E
amplitude = 1

makedir('barrier1')
os.chdir('barrier1')



os.chdir('..')

# barrier height<E
amplitude = 0.5

makedir('barrier3')
os.chdir('barrier3')



os.chdir('..')


# barrier height>E
amplitude = 2

makedir('barrier3')
os.chdir('barrier3')



os.chdir('..')





os.chdir('..')