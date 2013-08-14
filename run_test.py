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

class TestSettings:
    def __init__(self, name, desc, settings=[]):
        self.name = name
        self.desc = desc
        self.settings = settings
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

available_tests = [
    TestSettings(
        name='noise', desc='Noise',
        settings=[{'noise': val} for val in linspace(0.0,0.5,6)]),
    TestSettings(
        name='supports', desc='Support Beams On/Off',
        settings=[{'supports': False}]),
    TestSettings(
        name='rotate', desc='Rotation Angle',
        settings=[{'rotate': val} for val in linspace(0,45,10)]),
    TestSettings(
        name='ell_u', desc='Eccentricity',
        settings=[{'ell_u': val} for val in linspace(1.0,1.05,6)]),
    TestSettings(
        name='offset_u', desc='Pointing Offset',
        settings=[{'offset_u': val} for val in linspace(0,2.0,10)]),
    TestSettings(
        name='point_eccentricity', desc='Pointing Offset & Eccentricity',
        settings=[{'offset_u': offset_u, 'ell_u': ell_u} 
            for offset_u, ell_u in zip(
                    linspace(0,2.0,10), 
                    [1.00, 1.01, 1.02, 1.03, 1.04, 1.05])])
]

def run_test(requested_tests=[]):
    for test in available_tests:
        if not requested_tests or test.name in requested_tests:
            test.run()

