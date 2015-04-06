import os
import numpy as np
import TSDE

# change directory to where the data will be stored
os.chdir('runs')

# set up testing conditions
initWaveFunc = ?
delx = np.array([testing different delx?])
delt = np.array([testing different delt?])
timesteps = timesteps

#----------------------------------------
# Non-periodic No V
#----------------------------------------
os.makedirs('nonperiodic_noV')
os.chdir('nonperiodic_noV')
potential = 0
outputFile = "noV"
funstuff=TDSE.TDSE(initWaveFunc, potential, delx,
                   delt, timesteps, outputFile)
funstuff.runNonPeriodicNoV()
os.chdir('..')



#----------------------------------------
# Non-periodic With V
#----------------------------------------
os.makedirs('nonPeriodic_V')
os.chdir('nonPeriodic_V')
potentials = np.array(['Andrews output of different potentials'])
for x in potentials:
    potential = x
    dirname = '%.d' % (x)
    os.mkdir(dirname)
    os.chdir(dirname)
    outputFile = "nonPeriodic"
    funstuff=TDSE.TDSE(initWaveFunc, potential, delx,
                       delt, timesteps, outputFile)
    funstuff.runNonPeriodicWithV()
    os.chdir('..')
os.chdir('..')



#----------------------------------------
# Periodic With V
#----------------------------------------
os.makedirs('periodic_V')
os.chdir('periodic_V')
potentials = np.array(['Andrews output of different potentials'])
for x in potential:
    potential = x
    dirname = '%.d' % (x)
    os.mkdir(dirname)
    os.chdir(dirname)
    outputFile = "periodic"
    funstuff=TDSE.TDSE(initWaveFunc, potential, delx,
                       delt, timesteps, outputFile)
    funstuff.runPeriodicWithV()
    os.chdir('..')
os.chdir('..')
