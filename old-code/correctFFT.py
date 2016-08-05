from scipy import *
from scipy import fftpack

import matplotlib.pyplot as pl

########## Unit conversion
##  1.0 / ( cell size in radians ) = ( u_max in units of wavelength )
##  1.0 / ( field of view in radians ) = ( cell size in uv-domain in units of wavelength )

cell = 0.1 # arcmin
fov = 20.0 # arcmin

# Image pixels
xvals = arange(-1*fov/2,fov/2,cell)
N = xvals.shape[0]

# UV cell size and extent - in units of wavelengths
uvcell = 1/ (fov/60.0 * 3.14158/180.0)
uvmax = N/2 * uvcell

uvals = arange(-1*uvmax, uvmax, uvcell)
N2 = uvals.shape[0]

print 'cell (arcmin): ', cell
print 'fov (arcmin): ', fov
print 'uvcell (lambda): ',uvcell
print 'uvmax (lambda): ', uvmax

print N, N2  # To make sure they're equal !

# Source size
s = 0.05

### Source location
### Move this around 0.0, 0.1, 0.2, etc and see what happens to the phase...
x0 = 0.5

# data
data = exp ( -1 * (xvals-x0)**2 / (s*s) )

#data.fill(0.0)
#data[100+x0/cell]=1.0

# FT of data
fdata = fftpack.fftshift( fftpack.fft( fftpack.ifftshift( data ) ) )

# Inverse FT
ifdata = fftpack.fftshift( fftpack.ifft( fftpack.ifftshift( fdata ) ) )

# Plot everything...
pl.clf()
pl.subplot(311)
pl.title("Signal domain")
pl.plot(xvals,data,'b')
pl.plot(xvals,real(ifdata),'g')

pl.subplot(312)
pl.title("Fourier domain - amp")
pl.plot(uvals, abs(fdata) )

pl.subplot(313)
pl.title("Fourier domain - phase (radians)")
pl.plot(uvals, arctan2( imag(fdata), real(fdata) ) )
pl.axis( [ min(uvals), max(uvals), -3.2, 3.2] )

pl.show()
