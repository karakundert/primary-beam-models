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
    def __init__(self, name, desc, values=[], settings=None):
        self.name = name
        self.desc = desc
        self.values = values
        self.runnum = 0

    def run(self):
        for i, val in enumerate(self.values):
            print self.desc
            print "Run number =", self.runnum
            print "%s: %s" % (self.name, val[i])
            #makeMS(runnum="%03d" % self.runnum, **val[i])
            print "makeMS(runnum=%03d, %s) % (self.runnum, val[i])
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
        settings=[{'rotate': val} for val in linspace(0,90,19)]),
    TestSettings(
        name='ell_u', desc='Eccentricity',
        settings=[{'ell_u': val} for val in linspace(1.0,1.05,6)]),
    TestSettings(
        name='offset_u', desc='Pointing Offset',
        settings=[{'ell_u': val} for val in linspace(0,2.0,10)]),
    )
    TestSettings(
        name='rot_eccentricity', desc='Rotation & Eccentricity',
        settings=[{'ell_u': ell_u, 'theta': theta} for ell_u, theta in zip(linspace(0,2.0,10)]), [1.00, 1.01, 1.02, 1.03, 1.04, 1.05])
    )
]

def run_test(requested_tests=[]):
    for test in available_tests:
        if not requested_tests or test.name in requested_tests:
            test.run()

'''
def run_test(runnum=0):

        for i in linspace(0,45,6):
            eccen = [1.00, 1.01, 1.02, 1.03, 1.04, 1.05]
            print "Rotated + Eccentricity"
            print "Run Number = ", runnum
            theta = i
            ell_u = eccen[runnum]
            print "rotation angle = "+ str(i)
            print "eccentricity = "+ str(ell_u)
            makeMS(runnum,noise=False,ell_u=ell_u,theta=theta)
            runnum += 1

        os.system('rm -rf rot_eccentricity')
        os.system('mkdir rot_eccentricity')
        os.system('mv Data* rot_eccentricity')
'''
