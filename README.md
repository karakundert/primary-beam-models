Simulation
==========

This simulation is intended to be used to demonstrate the effects of
distortions in the primary beam on radio astronomical observations. By changing
the parameters, the user can determine the point at which these distortions
will cross the noise threshold, and become visible to the observer.

The user can make changes to the parameters of the beam distortion in run_test.py.

The CASA package must be installed in order to be able to run this simulation.

To run the program, type these commands into the CASA terminal.

execfile('makesources.py')
execfile('simsky.py')
execfile('run_test.py')
run_test()
