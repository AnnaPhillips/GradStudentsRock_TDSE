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
potential=0




#----------------------------------------
# Non-periodic With V
#----------------------------------------
os.makedirs('nonperiodic_V')
os.chdir('nonperiodic_V')
potentials = np.array([any potentials])
for x in potential:
    dirname = '%.d' % (x)
    os.mkdir(dirname)
    os.chdir(dirname)
    outputFile = "nonperiodic"
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
potential = np.array([any potentials])
potentials = np.array([any potentials])
for x in potential:
    dirname = '%.d' % (x)
    os.mkdir(dirname)
    os.chdir(dirname)
    outputFile = "periodic"
    funstuff=TDSE.TDSE(initWaveFunc, potential, delx,
                       delt, timesteps, outputFile)
    funstuff.runPeriodicWithV()
    os.chdir('..')
os.chdir('..')
