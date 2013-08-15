execfile('simsky.py')

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
        for level in ['centered', 'half-power,' 'low-power']:
            if callable(self.transform):
                xfrm = self.xtransform
            elif self.transform:
                xfrm = lambda x: x
            else:
                xfrm = lambda x: x*self.transform
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
    Test(
        name='noise', desc='Noise', xaxis='noise', xtransform=100,
        settings=[{'noise': val} for val in linspace(0.0,0.5,6)],
        xlabel="Percent Increase in Semimajor Axis"),
    Test(
        name='supports', desc='Support Beams On/Off',
        settings=[{'supports': False}]),
    Test(
        name='rotate', desc='Rotation Angle', xaxis='rotate',
        settings=[{'rotate': val} for val in linspace(0,45,10)],
        xlabel="Percent Increase in Semimajor Axis"),
    Test(
        name='eccentricity', desc='Eccentricity', xaxis='ell_u', xtransform=lambda x: 100*(x-1),
        settings=[{'ell_u': val} for val in linspace(1.0,1.05,6)],
        xlabel="Percent Increase in Semimajor Axis"),
    Test(
        name='pointing', desc='Pointing Offset', xaxis='offset_u', xtransform=1/6.8,
        settings=[{'offset_u': val} for val in linspace(0,2.0,10)],
        xlabel="Percent Increase in Semimajor Axis"),
    Test(
        name='point_eccentricity', desc='Pointing Offset & Eccentricity',
        xaxis='offset_u',
        settings=[{'offset_u': offset_u, 'ell_u': ell_u} 
            for offset_u, ell_u in zip(
                    linspace(0,2.0,10), 
                    [1.00, 1.01, 1.02, 1.03, 1.04, 1.05])],
        xlabel="Percent Increase in Semimajor Axis")
]

def run_test(tests=[], run=True, plot=True):
    assert run or (plot if self.axis else False)
    for test in available_tests:
        if not tests or test.name in tests:
            if run:
                test.run()
            if plot:
                test.plot()
