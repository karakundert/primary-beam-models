from numpy import *
from matplotlib import *
from scipy import *

from scipy import fftpack

import matplotlib.pyplot as pl

# Unit conversion
##  1.0 / ( cell size in radians ) = 
#	( u_max in units of wavelength )
##  1.0 / ( field of view in radians ) = 
#	( cell size in uv-domain in units of wavelength )

cell = 0.1 # arcmin
fov = 20.0 # arcmin
wvlen = 1 * 10**-9 # observed wavelength
d = 25 # dish width

# Image pixels
xvals = arange(-1*fov/2,fov/2,cell)
N = xvals.shape[0]

# UV cell size and extent - in units of wavelengths
uvcell = 1/ (fov/60.0 * pi/180.0)
uvmax = N/2 * uvcell

uvals = arange(-1*uvmax, uvmax, uvcell)
N2 = uvals.shape[0]

res = ( 1.0 / uvmax ) * ( 180.0 / pi ) * ( 60.0 / 1.0 ) 

print 'cell (arcmin): ', cell
print 'fov (arcmin): ', fov
print 'res (arcmin): ', res
print 'uvcell (lambda): ',uvcell
print 'uvmax (lambda): ', uvmax

print N, N2  # To make sure they're equal !

########################

s = 0.5
x0 = 3.0

# data
#data = exp ( -1 * (xvals-3.0)**2 / (s*s) )
data = zeros(len(xvals))
data[27] = 1.0
data[70] = 1.0
data[100] = 1.0
data[112] = 1.0
data[130] = 1.0

N = data.shape[0]

# aperture function
aper = zeros(len(uvals))
d_uv = len(uvals) / 100.0
low = ( len(uvals) / 2 ) - d_uv
high = ( len(uvals) / 2 ) + d_uv

for i in arange(low,high,1):
	aper[i] = 1.0

aper = aper / aper.sum()

# aperture autocorrelation
auto_corr = numpy.convolve(aper,aper,'same')
auto_corr = fftpack.fftshift(auto_corr)

pl.clf()
pl.plot(auto_corr)

# power pattern function
power = N * fftpack.ifft(auto_corr)
power = fftpack.ifftshift(power)

# multiply true brightness with power pattern
# save to data
power_data = data * power

# Sampling function
spoints = (len(xvals) / 2) + 1 
samples = zeros ( data.shape )
locs = random.uniform(0,1,spoints)*100
for i in range(0,spoints):
    samples[ int(locs[i]) ] = 1.0

samples[N/2:N] = samples[N/2:0:-1]

# FT of power pattern modified data
fdata = 1.0/N * fftpack.fft( power_data )
fdata = fftpack.fftshift( fdata )

# Sample FT of data
sfdata = samples * fdata

# Inverse FT of sampled data
ifdata = N * fftpack.ifft( fftpack.ifftshift( sfdata ) )

# PSF
psf = 1.0/N * fftpack.fft(samples)
psf = fftpack.fftshift(psf)

# Plot everything...
pl.clf()
pl.subplot(321)
pl.title("True Signal and Power-Modified Data")
pl.plot(xvals, data)

pl.subplot(322)
pl.title("Power Pattern Modified Data")
pl.plot(xvals, power_data)

pl.subplot(323)
pl.title("Antenna Power Pattern")
pl.plot(xvals, power)

pl.subplot(324)
pl.title("Fourier Signal")
pl.plot(uvals, abs(fdata))
pl.plot(uvals, abs(sfdata))

pl.subplot(325)
pl.title("Point Spread Function")
pl.plot(xvals, psf)

pl.subplot(326)
pl.title("Observed Signals")
pl.plot(xvals, ifdata)

pl.show()
