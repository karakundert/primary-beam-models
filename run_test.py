source_dir = "/lustre/kkundert/Code/"

#execfile('simsky.py')
execfile(source_dir+"alma_simsky.py")

from numpy import *
from scipy import *
from scipy import fftpack
from scipy import signal
from numpy import fft
from scipy import ndimage

import matplotlib.pylab as pl
import random
import os
import re
from glob import glob
from matplotlib import axes

# choose the comparable data set
numRegex = r"[+-]?([0-9]+\.)?[0-9]+([eE][+-][0-9]+)?"
requiredWhiteSpaceRegex = r"[ \t]+"
optionalWhiteSpaceRegex = r"[ \t]*"
rmsRegex = r"(?P<rms>%s)" % numRegex
newline = r"\n"
headerRegex = (
    optionalWhiteSpaceRegex
    + "#"
    + optionalWhiteSpaceRegex
    + "Npts"
    + requiredWhiteSpaceRegex
    + "Sum"
    + requiredWhiteSpaceRegex
    + "Mean"
    + requiredWhiteSpaceRegex
    + "Rms"
    + requiredWhiteSpaceRegex
    + "Std dev"
    + requiredWhiteSpaceRegex
    + "Minimum"
    + requiredWhiteSpaceRegex
    + "Maximum"
    + optionalWhiteSpaceRegex
)
regex = (
    headerRegex
    + newline
    + optionalWhiteSpaceRegex
    + 3 * (numRegex + requiredWhiteSpaceRegex)
    + rmsRegex
    + 3 * (requiredWhiteSpaceRegex + numRegex)
    + optionalWhiteSpaceRegex
)
pattern = re.compile(regex, re.IGNORECASE)

class Test:
    def __init__(self, name, desc, xaxis=None, settings=[], xlabel='', xtransform=None):
        self.name = name
        self.desc = desc
        self.settings = settings
        self.xaxis = xaxis
        self.xlabel = xlabel
        self.xtransform = xtransform
            # xtransform may be None,
            #     in which case the value plotted on the x axis is exactly the
            #     values used in the analysis,
            # or it may be a number,
            #     in which case the value used in the plot is the value used
            #     in the analysis scaled by xtransform,
            # otherwise xtransform can be a function,
            #     in which case the value used in the plot is the value returned
            #     by the function when the value used in the analysis is passed
            #     in as the only argument.
        self.runnum = 0

    def run(self):
        # calls makeMS function in simsky code to actually run the individual
        # tests
        # increments runnum to make tabulating data easy
        for val in self.settings:
            print self.desc
            print "Run number =", self.runnum
            print "%s: %s" % (self.name, val)
            makeMS(runnum="%03d" % self.runnum, **val)
            self.runnum += 1

        os.system('rm -rf %s' % self.name)
        os.system('mkdir %s' % self.name)
        os.system('mv Data* %s' % self.name)

    def plot(self):
        for level in ['single', 'multi']:
            # gemerates plots for all sources for given test run
            # code below specifies x-axis for plot
            if callable(self.xtransform):
                xfrm = self.xtransform
            elif self.xtransform:
                xfrm = lambda x: x*self.xtransform
            else:
                xfrm = lambda x: x
            x_axis = [xfrm(each[self.xaxis]) for each in self.settings]
            if len(x_axis) < 2:
                print "plotting of %s skipped" % self.name
                return
            print "\n###", level.upper()
            filenames = glob('%s/Data*-%s' % (self.name, level))
            rmsValues = []
            rmsValuesNear = []
            rmsValuesOff = []

            for filename in filenames:
                # generates plots for all statistical locations
                with open(filename+"/all_stats.txt") as f:
                    contents = f.read()
                    match = pattern.search(contents)
                    assert match
                    rmsValues.append(float(match.group('rms')))
                with open(filename+"/near_source_stats.txt") as f:
                    contents = f.read()
                    match = pattern.search(contents)
                    assert match
                    rmsValuesNear.append(float(match.group('rms')))
                with open(filename+"/off_source_stats.txt") as f:
                    contents = f.read()
                    match = pattern.search(contents)
                    assert match
                    rmsValuesOff.append(float(match.group('rms')))

            print rmsValues
            print rmsValuesNear
            print rmsValuesOff

            pl.clf()
            pl.title("%s - %s" % (self.desc, level))
            p1 = pl.plot(x_axis,rmsValues,'b')
            p2 = pl.plot(x_axis,rmsValuesNear,'g')
            p3 = pl.plot(x_axis,rmsValuesOff,'r')
            pl.xlabel(self.xlabel)
            pl.legend([p1,p2,p3],["Whole Sky", "Near Source", "Off Source"], loc = 2)
            pl.savefig("%s-%s.png" % (self.name, level))

available_tests = [
    # list of tests possible to be executed
    Test(
        name='no_perturbation', desc='No Perturbation',
        settings=[{'makeBeams': True}]),
    Test(
        name='no_perturbation_no_make', desc='No Perturbation',
        settings=[{'makeBeams': False}]),
    Test(
        name='noise', desc='Noise', xaxis='noise', xtransform=100,
        settings=[{'noise': val} for val in linspace(0.3,0.5,3)],
        xlabel="Percent Noise Across Aperture"),
    Test(
        name='supports', desc='Support Beams On/Off',
        settings=[{'supports': False}]),
    Test(
        name='rotate', desc='Rotation Angle', xaxis='rotate',
        settings=[{'rotate': val} for val in linspace(0,45,10)],
        xlabel="Degree of Beam Rotation"),
    Test(
        name='eccentricity', desc='Eccentricity', xaxis='ell_u', xtransform=lambda x: 100*(x-1),
        settings=[{'ell_u': val} for val in linspace(1.0,1.05,6)],
        xlabel="Percent Increase in Semimajor Axis"),
    Test(
        name='pointing', desc='Pointing Offset', xtransform=1/6.8,
        settings=[{'pointing': True}]),
    Test(
        name='pointing_no_make', desc='Pointing Offset',
        xtransform=1/6.8,
        settings=[{'pointing': False, 'makeBeams': False},
            {'pointing': True, 'makeBeams': False}]),
    Test(
        name='rot_aper', desc='With/Without Rotated Apertures',
        settings=[{'rot': False}, {'rot': True}]),
    Test(
        name='point_eccentricity', desc='Pointing Offset & Eccentricity',
        xaxis='offset_u',
        settings=[{'offset_u': offset_u, 'ell_u': ell_u}
            for offset_u, ell_u in zip(
                    [True, False],
                    [1.00, 1.01, 1.02, 1.03, 1.04, 1.05])])
]

def run_test(tests=[], run=True, plot=True):
    # function that runs the tests and makes the plots, as defined by the
    # function parameters
    # specify which tests to run by filling tests list, or default to running
    # all tests
    assert run or plot
    for test in available_tests:
        if not tests or test.name in tests:
            if run:
                test.run()
            if plot:
                test.plot()
