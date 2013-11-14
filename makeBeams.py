###
### Functions for making apertures and primary beams for use in simulation.
###

from numpy import *
from scipy import *
from scipy import fftpack
from scipy import signal
from numpy import fft
from numpy import max
from scipy import ndimage

import matplotlib.pylab as pl
import random
import time

###############################################

def imageFromArray(arr,outfile,coord={},linear=F):
    newia = casac.image()
    newia.fromarray(outfile=outfile,pixels=arr,csys=coord,linear=F,overwrite=T)
    newia.close()

###############################################

def perturbAperture(image="model",imsize=256,cellsize='8.0arcsec',
                    reffreq='1.5GHz', d = 25.0,
                    aper_rreal = 'reaperture_pol5.im',
                    aper_rimag = 'imaperture_pol5.im',
                    aper_lreal = 'reaperture_pol8.im',
                    aper_limag = 'imaperture_pol8.im',
                    noise = 0.0,
                    pointing = True):

        # noise == True: there is Gaussian noise in the aperture function
        # supports == True: there are shadows from the support beams

        # ell_u and ell_v describe how elliptical the aperture is. If both are
        # equal to each other, the aperture will be a circle. Any relative
        # change in ell_u and ell_v will result in an elliptical aperture.

        cell = qa.convert(qa.quantity(cellsize),'arcmin')['value']
        fov = cell*imsize
        freq = qa.convert(qa.quantity(reffreq),'Hz')['value']
        c = 3e8
        j = sqrt(-1)
        wvlen = c / freq # observed wavelength
        d = 25 # width of a dish in m
        spat_lam = d / wvlen

        # Image pixels
        xvals = arange(-1*fov/2.0,fov/2.0,cell)
        yvals = arange(-1*fov/2.0,fov/2.0,cell)
        Nxy = xvals.shape[0]

        # UV cell size and extent - in units of wavelengths
        uvcell = 1 / (fov/60.0 * pi/180.0)
        uvmax = Nxy/2.0 * uvcell
        d_uv = spat_lam / uvcell

        uvals = arange(-1*uvmax, uvmax, uvcell)
        vvals = arange(-1*uvmax, uvmax, uvcell)
        Nuv = uvals.shape[0]

        print 'cell (arcmin): ', cell
        print 'fov (arcmin): ', fov
        print 'uvcell (lambda): ',uvcell
        print 'uvmax (lambda): ', uvmax

        print Nxy, Nuv  # To make sure they're equal !

        ###################

        # Aperture Function
        uu, vv = mgrid[:Nuv, :Nuv]
       
        aper_r = zeros((Nuv,Nuv), 'complex')
        aper_l = zeros((Nuv,Nuv), 'complex')

        ia.open(aper_rreal)
        aper_r_real = ia.getchunk()
        ia.close()

        ia.open(aper_rimag)
        aper_r_imag = ia.getchunk()
        ia.close()

        ia.open(aper_lreal)
        aper_l_real = ia.getchunk()
        ia.close()

        ia.open(aper_limag)
        aper_l_imag = ia.getchunk()
        ia.close()

        real_aper_r = zeros((2048,2048))
        imag_aper_r = zeros((2048,2048))
        real_aper_l = zeros((2048,2048))
        imag_aper_l = zeros((2048,2048))

        real_aper_r = aper_r_real[:,:,0,0]
        imag_aper_r = aper_r_imag[:,:,0,0]
        real_aper_l = aper_l_real[:,:,0,0]
        imag_aper_l = aper_l_imag[:,:,0,0]

        #low = Nuv / 2 - Nuv / 4
        #high = Nuv / 2 + Nuv / 4
        #print low, high
        aper_r[:,:] = real_aper_r + j * imag_aper_r
        aper_l[:,:] = real_aper_l + j * imag_aper_l
        
        imageFromArray(imag(aper_r),"squint_H")
        imageFromArray(imag(aper_l),"squint_V")

        # Add phase ramp to aperture function
        # offset measured in arcseconds. please remember before choosing offset
        # values for VLA.
        max_offset_VLA = 60*0.10*wvlen/d
        max_offset_ALMA_uncor = 2.0
        max_offset_ALMA_cor = 0.5
        max_offset = max_offset_ALMA_cor
        if pointing == True:
            offset_u = random.uniform(-max_offset,max_offset)
            offset_v = random.uniform(-max_offset,max_offset)
        else:
            offset_u = 0.0
            offset_v = 0.0
        shift_u = offset_u # in arcsecs
        shift_v = offset_v # in arcsecs
        phs = zeros((Nuv,Nuv))
        phs_u = -1 * (uvals[uu] + uvcell/2.0) * ( shift_u / 60.0 / 60.0 * ( pi / 180.0 ) ) * 2 * pi
        phs_v = -1 * (vvals[vv] + uvcell/2.0) * ( shift_v / 60.0 / 60.0 * ( pi / 180.0 ) ) * 2 * pi

        real_noise = noise * randn(Nuv,Nuv) + 1
        imag_noise = noise * randn(Nuv,Nuv) + 1
        time1 = time.time()
        aper_r = aper_r * abs(real_noise)
        aper_l = aper_l * abs(real_noise)
        time2 = time.time()
        phs_u = phs_u * imag_noise
        phs_v = phs_v * imag_noise
        time2 = time.time()
        print "noise time = " + str(time2-time1)
        phs = exp(j*(phs_u + phs_v))

        phs_aper_r = aper_r * phs
        phs_aper_l = aper_l * phs

        imageFromArray(real(phs_aper_r),image+"H")
        imageFromArray(real(phs_aper_l),image+"V")

        volt_r = fftpack.ifftshift(fftpack.fft2(fftpack.fftshift(phs_aper_r)))
        volt_l = fftpack.ifftshift(fftpack.fft2(fftpack.fftshift(phs_aper_l)))
        
        imageFromArray(real(volt_r),image+"Hreal")
        imageFromArray(imag(volt_r),image+"Himag")
        imageFromArray(real(volt_l),image+"Vreal")
        imageFromArray(imag(volt_l),image+"Vimag")

        del aper_r,aper_l,phs,phs_u,phs_v,volt_r,volt_l,\
            phs_aper_r,phs_aper_l,uu, vv, uvals, vvals

        return Nxy

###############################################

def makePrimaryBeam_VLA(imsize=256,cellsize='8.0arcsec',coord="coord.torecord()",
                    reffreq='1.5GHz', pbname = "model",
                    aper1_Rreal = "aper00Hreal", aper1_Rimag = "aper00Himag",
                    aper2_Rreal = "aper01Hreal", aper2_Rimag = "aper01Himag",
                    aper1_Lreal = "aper00Vreal", aper1_Limag = "aper00Vimag",
                    aper2_Lreal = "aper01Vreal", aper2_Limag = "aper01Vimag",
                    Nxy = 200, area = -1):

        power_r = zeros((imsize,imsize))
        power_l = zeros((imsize,imsize))

        i = 0
        while i < 2:
            if i == 0:
                aper1_real = aper1_Rreal
                aper1_imag = aper1_Rimag
                aper2_real = aper2_Rreal
                aper2_imag = aper2_Rimag
            else:
                aper1_real = aper1_Lreal
                aper1_imag = aper1_Limag
                aper2_real = aper2_Lreal
                aper2_imag = aper2_Limag
                
            # Combine aperture files to make complex aperture array
            ia.open(aper1_real)
            a1_real = ia.getchunk()
            ia.close()

            ia.open(aper1_imag)
            a1_imag = ia.getchunk()
            ia.close()

            ia.open(aper2_real)
            a2_real = ia.getchunk()
            ia.close()

            ia.open(aper2_imag)
            a2_imag = ia.getchunk()
            ia.close()

            time1 = time.time()
            aper1 = a1_real + (sqrt(-1) * a1_imag)
            aper2 = a2_real + (sqrt(-1) * a2_imag)
            power = aper1 * aper2
            time2 = time.time()
            print "convolution time = " + str(time2-time1)

            # Power pattern function
            power = power / max(power)

            if i == 0:
                power_r = power
            else:
                power_l = power

            i = i + 1

        stokes_power = zeros((imsize,imsize,4,1))
        stokes_power[:,:,0,0] = (power_r[:,:] + power_l[:,:]) / 2
        stokes_power[:,:,3,0] = (power_r[:,:] - power_l[:,:]) / 2

        imageFromArray(real(stokes_power),pbname+"real",coord)
        imageFromArray(imag(stokes_power),pbname+"imag",coord)

        aper_area = -1

        del aper1_real, aper1_imag, aper2_real, aper2_imag, power,\
            power_r, power_l, stokes_power
        
###############################################
