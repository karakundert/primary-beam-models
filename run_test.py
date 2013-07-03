execfile('simsky.py')

from numpy import *
from scipy import *
from scipy import fftpack
from scipy import signal
from numpy import fft
from scipy import ndimage

import matplotlib.pylab as pl
import random


def run_test(runnum=0):

        #print "Noise On/Off"
        #print "Run Number = ", runnum

        #makeMS(runnum,noise=False)
        #runnum += 1

        #print "Support Beams On/Off"
        #print "Run Number = ", runnum

        #makeMS(runnum,supports=False)
        #runnum += 1

        for i in linspace(1.0,1.05,5):
            print "Eccentricity"
            print "Run Number = ", runnum
            ell_u = i
            print i
            makeMS(runnum,noise=False,ell_u=ell_u)
            runnum += 1
