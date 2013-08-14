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
        for val in self.values:
            print self.desc
            print "Run number =", self.runnum
            self.runnum += 1
            print self.name, "=", val
            arg = {self.name: val}
            if settings:
                settings.add(arg)
            else
                settings = arg
            makeMS(runnum="%03d" % self.runnum, **args)

        os.system('rm -rf %s' % self.name)
        os.system('mkdir %s' % self.name)
        os.system('mv Data* %s' % self.name)

available_tests = [
    TestSettings(name='noise', desc='Noise', values=linspace(0.0,0.5,6)),
    TestSettings(name='supports', desc='Support Beams On/Off', values=[False]),
    TestSettings(name='rotate', desc='Rotation Angle', values=linspace(0,90,19)),
    TestSettings(name='ell_u', desc='Eccentricity', values=linspace(1.0,1.05,6)),
    TestSettings(name='offset_u', desc='Pointing Offset',
        values=linspace(0,2.0,10)),
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
