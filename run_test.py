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


def run_test(runnum=0):

        for i in linspace(0.0,0.5,6):
            print "Noise"
            print "Run Number = ", runnum
            noise = i
            print "gaussian noise = "+ str(i)
            makeMS(runnum,noise=noise)
            runnum += 1

        os.system('mv Data* noise')

        print "Support Beams On/Off"
        print "Run Number = ", runnum

        makeMS(runnum,supports=False)
        runnum += 1

        os.system('mv Data* noise')

        runnum = 20

        for i in linspace(0,45,8):
            print "Rotation Angle"
            print "Run Number = ", runnum
            theta = i
            print "rotation angle = "+ str(i)
            makeMS(runnum,noise=False,theta=theta)
            runnum += 1

        os.system('mv Data* noise')

        runnum = 30

        for i in linspace(1.0,1.05,6):
            print "Eccentricity"
            print "Run Number = ", runnum
            ell_u = i
            print "eccentricity = "+ str(i)
            makeMS(runnum,noise=False,ell_u=ell_u)
            runnum += 1

        os.system('mv Data* noise')

        runnum = 40

        for i in linspace(0,2.0,10):
            print "Phase Offset"
            print "Run Number = ", runnum
            offset = str(i)+"arcmin"
            print "offset = "+ offset
            makeMS(runnum,noise=False,offset_u=offset)
            runnum += 1

        os.system('mv Data* noise')
