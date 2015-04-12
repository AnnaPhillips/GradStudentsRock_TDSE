# The Time Dependent Schroedinger Equation
### Computational Physics -- Tufts University, Sprint 2015
### Andrew DeBendedictis, Anna Phillips, & Jennifer Radoff

#### Overview
This program runs simulations of the time dependent schroedinger equation for 
different potentials. This program runs a simple finite difference scheme and the 
Crank-Nicolson scheme. It also runs for either periodic or non-periodic boundary 
conditions. It runs 8 different potential conditions: zero potential (free 
particle), square well, harmonic oscillator, triangle, Kronig-Penney, a single 
barrier, and two different imaginary potentials: V=ix and V=x+ix.

#### Usage & Options
This program was developed in Python 2.7.x, with Numpy 1.9. 

We created 3 different main programs, `TDSE.py`, `potentials.py`, and 
`create_wave.py`, and 1 program that coordinates the 3 main scripts and speficies 
the run conditions, `RunTests.py`.

________________________

`potentials.py` define our potentials. In Python, the user can import 
`potentials.py` and call up partiuclar potentials. The potentials can be run using 
this string of code once imported in Python:

`potentials.Potentials(xMin, xMax, gridpoints, xAR).POTENTIALFUNCTION(args)`

where `xMin` and `xMax` are the minimum and maximum spatial points, `gridpoints` is 
the number of spatial grid points and `xAR` is an overall amplitude factor for the 
potential. Every potential accepts an amplitude (except for the free particle), and 
potentials with special argument values are as follows:

*   `barrier` accepts a barrier width
*   `squareWell` accepts a barrier width
*   `kronigPenney` accepts a gap width and number of wells

Each potential outputs a 1xN array the size of our grid (which we use in `TDSE.py`
to create an A matrix for our linear solver).

________________________

`create_wave.py` creates an initial wave function to pass through our potential. The 
user can call either a mesh array or a potential in a script by inputting: 

`create_wave.xMesh(xMin, xMax, gridpoints)` or `create_wave.WAVETYPE()`

The script creates a mesh, which is an array of points 1 longer than our number of 
gridpoints to normalize the space for the wave. It uses the mesh function to create 
a gaussian wave packet, a consine wave, and an exponential wave.

`xMesh` accepts a minimum and maximum point, as well as number of gridpoints. Each 
initial `WAVETYPE` accepts the `xMesh` array in order to normalize it to the 
grid. In addition, the `gaussian` accepts an initial position, an initial wave 
number, and a width. The `cosine` wave accepts a harmonic number, and a well width 
(it is used for the well potential). The `exponential` accepts an amplitude.

_________________________

`TDSE.py` creates an A matrix for any given potential for both the simple finite 
difference scheme and the Crank-Nicolson scheme, for either periodic or non-periodic 
boundary conditions. It then solves a system of linear equations which updates a 
matrix that describes how the wave function evolves over time. The real and the 
imaginary parts of this matrix are output to seperate .csv files. It also calculates 
the eigenvalues and eigenvectors for each A matrix and outputs both the eigenvalues 
and real and imaginary parts of the eigenvectors in seperate .csv files. 

In Python, the user can import `TDSE.py` and run the script using this string of 
code:

`TDSE.TDSE(initWaveFunc, potential, delx, delt, timesteps, periodic, useCN, 
outputFile)`

where `initWaveFunc` is the intial wave function array created by `create_wave.py`, 
`potential` is the potential array created by `potentials.py`, `delx` and `delt` are 
the position-step and time-step, `timesteps` is the number of total time steps, 
`periodic` specifies periodic or non-periodic conditions (where the input `True` 
calls periodic conditions and the input `False` calls non-periodic conditions), 
`useCN` specifies the scheme (where the input `True` calls Crank-Nicolson scheme and 
the input `False` calls the simple finite difference scheme), and `outputFile` 
specifies the name of the output file and should be input as a string without 
specifying a file type (i.e. "filename").

__________________________

`RunTests.py` creates directories for the different potential cases and runs 
`TDSE.py` for each potential and at different conditions (with different types of 
wave functions and with different potential and wave function amplitudes, widths, 
etc.). The wave function matrix output is arranged into directories named as
follows:

`Runs/<POTENTIAL>/<periodic_condition>/_<scheme>/_<real_or_imag>.csv`

and the eigenstate matrices are arranged as follows:

`Runs/<POTENTIAL>/<periodic_condition>/_<scheme>/_<real_or_imag>/_<eValues_or_eVectors>.csv`





