Simulation
==========

This simulation is intended to be used to demonstrate the effects of
distortions in the primary beam on radio astronomical observations. By changing
the parameters, the user can determine the point at which these distortions
will cross the noise threshold, and become visible to the observer.

The user can make changes to the parameters of the beam distortion in run_test.py.

The CASA package must be installed in order to be able to run this simulation.

To run the simulation, type these commands into the CASA terminal.

    execfile('init.py')
    run()

You may want to choose tests to run individually, as opposed to running all
tests together. The simulation takes several hours to run with all tests. To do
so, look at the source code for run_test.py and choose the tests desired.

The test is currently configured for ALMA parameters. The test can also be run
with VLA parameters. To do so, open run_test.py and change the execfile to
simsky.py instead of alma_simsky.py.

Make sure before running the simulation that the source directory has been
changed in all files to match the location of each file. The simulation can
currently be run from any Socorro NRAO computer, but will require
reconfiguration for any other machine.

Plotting the Results
====================

To plot the results of an ALMA simulation, type these commands into the CASA
terminal.

    execfile('make_cor_image.py')

This will make plots of the rms-levels and image fidelity of all runs. Make
sure that the directory list in the code matches the tests you want to plot,
and the pathnames of those tests.
