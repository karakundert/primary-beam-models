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

Plotting the Results
====================

Using the run_test() command, plots depicting trends of perturbations will
automatically be generated. However, it is possible - if one simply wants to
make plots on preexisting data sets. Look at the init.py source code to find
the necessary commands for plotting the right data sets.
