Simulation
==========

This simulation is intended to be used to demonstrate the effects of
distortions in the primary beam on radio astronomical observations. By changing
the parameters, the user can determine the point at which these distortions
will cross the noise threshold, and become visible to the observer.

The user can make changes to the parameters of the beam distortion in run_test.py.

The CASA package must be installed in order to be able to run this simulation.

To run the simulation, type these commands into the CASA terminal.

    execfile('makesources.py')
    execfile('simsky.py')
    execfile('run_test.py')
    run_test()

You may want to choose tests to run individually, as opposed to running all
tests together. The simulation takes several hours to run with all tests.

Plotting the Results
====================

By using the plot_data script, you can see the effects that these
perturbations have on the primary beam, and the data we collect. 

To run the data plotting script, type these commands into the CASA terminal.

    execfile('plot_data.py')

Be careful with how the run numbers combine in an individual test. If you
plot a group of points that spans from Data8-Data13, the plotting script will
place Data 10-13 on the plot before Data 8-9, making it difficult to see
trends.
