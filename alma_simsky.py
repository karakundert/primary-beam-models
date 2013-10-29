###
### Simulate Sky and Data...
###

from numpy import *
from scipy import *
from scipy import fftpack
from scipy import signal
from numpy import fft
from numpy import max
from scipy import ndimage
from simutil import *

import matplotlib.pylab as pl
import random
import time

execfile('alma_mksrc.py')
execfile('parang.py')
execfile('writeAps.py')
#execfile('makeBeams.py')

util = simutil("")

###############################################

def makeMS(runnum=0, makeBeams = True,
           noise=0.0, supports=True,
           ell_u = 1.0, ell_v = 1.0,
           pointing = False, theta = 0.0,
           rot = False):


  for i in xrange(2):
      dirname = "Data"+str(runnum)
      if i == 0:
          dirname = dirname+"-centered"
      else:
          dirname = dirname+"-multi"
      basename = dirname+"/points"
      msname = basename + '.ms';
      imname = basename+'.true.im';
      predict_imname = imname + '.predict'
      resname = dirname+"/pb-residuals-"+str(runnum)
      pbname = dirname+"/primary-beams"
      apname = dirname+"/apertures"
      clname = "mysources"
      ra0="19:59:28.500";
      dec0="-23.44.01.50";
      nchan=1;
      imsize=2048;
      predict_imsize = 1024
      cellsize='1.38arcsec';
      reffreq='100.0GHz';
      stokesvals=[1.0,1.0,0.0,0.0]
      ftm='ft'

      num_ant = 10

      rot_init = rot

      #addNoise(msname);
      if makeBeams == True:
          os.system('rm -rf '+dirname)
          os.system('rm -rf '+resname)
          os.system('rm -rf theresult.*')

          clname = clname+str(i)+'.cl'
          d = makeMSFrame(dirname=dirname,msname=msname,
                      ra0=ra0,dec0=dec0,nchan=nchan,numants=num_ant);

          os.system('mkdir '+apname)
          os.system('mkdir '+pbname)

          makeImage(msname=msname,imname=imname,clname=clname,
                    imsize=imsize,cellsize=cellsize,ra0=ra0,dec0=dec0,
                    nchan=nchan,reffreq=reffreq)

          ia.open(imname)
          coord = ia.coordsys()
          ia.close()

          #os.system('rm -rf ap.ms')
          #os.system('cp -r '+msname+' ap.ms')
          #writeApertures(msname='ap.ms',imname=imname)
          #print "apertures written"
          
          Nxy_list = []
          for i in xrange(num_ant):
              if i % 2:
                  rot = False
              print "aper" + str(i)
              print rot
              image = apname+"/aper%02d" % i
              Nxy = makeAperture(image=image,imsize=imsize,
                              cellsize=cellsize,reffreq=reffreq,
                              d=d[i],noise=noise,supports=supports,
                              ell_u=ell_u,ell_v=ell_v,
                              pointing=pointing,rot=rot)
              Nxy_list.append(Nxy)
              rot = rot_init
          for i in xrange(num_ant):
              for j in xrange(num_ant):
                  if i < j:
                      pbimage = pbname+"/pb%02d&&%02d" % (i,j)
                      print "primary beam"+str(i)+"&&"+str(j)
                      makePrimaryBeam(imsize=imsize,cellsize=cellsize,
                                    coord = coord.torecord(),
                                    reffreq=reffreq,pbname=pbimage,
                                    aper1_Rreal=apname+"/aper%02dXreal" % i,
                                    aper1_Rimag=apname+"/aper%02dXimag" % i,
                                    aper2_Rreal=apname+"/aper%02dXreal" % j,
                                    aper2_Rimag=apname+"/aper%02dXimag" % j,
                                    aper1_Lreal=apname+"/aper%02dYreal" % i,
                                    aper1_Limag=apname+"/aper%02dYimag" % i,
                                    aper2_Lreal=apname+"/aper%02dYreal" % j,
                                    aper2_Limag=apname+"/aper%02dYimag" % j,
                                    Nxy=Nxy_list[i]);
                  else:
                      continue
      else:
          clname = clname+str(i)+'.cl'
          makeMSFrame(dirname=dirname,msname=msname,
                      ra0=ra0,dec0=dec0,nchan=nchan,numants=num_ant);
          makeImage(msname=msname,imname=imname,clname=clname,
                    imsize=imsize,cellsize=cellsize,ra0=ra0,dec0=dec0,
                    nchan=nchan,reffreq=reffreq)

      makeImage(msname=msname,imname=predict_imname,clname=clname,
                imsize=predict_imsize,cellsize=cellsize,ra0=ra0,dec0=dec0,
                nchan=nchan,reffreq=reffreq)

      for i in xrange(num_ant):
          for j in xrange(num_ant):
              if i < j:
                  pair = "%02d&&%02d" % (i,j)
                  newimname = imname+"%02d&&%02d" % (i,j)
                  pbimage = pbname+"/pb%02d&&%02d" % (i,j)
                  pbreal = pbimage+"real"
                  pbimag = pbimage+"imag"
                  print pair
                  makeTrueImage(stokesvals=stokesvals,imname=predict_imname,
                                newimname=newimname,
                                pb_real_file = pbreal,pb_imag_file = pbimag,
                                msname=msname,ftm=ftm,imsize=imsize,
                                predict_imsize=predict_imsize,
                                cellsize=cellsize,ra0=ra0, dec0=dec0,
                                nchan=nchan, reffreq=reffreq, num_ant=num_ant,
                                pair=pair, model=False);
              else:
                  continue

      makeDataTable(msname=msname)
      makeResidualImage(msname, resname, imsize, cellsize, ra0, dec0, nchan,
              reffreq)
      clname = 'mysources'

      # The sources are located in one quadrant of the sky, and therefore one
      # quadrant of the primary beam. The statistics are chosen such that off
      # source is on the far side of the opposite quadrant from the sources, and
      # the near source stats are on the axis in the adjacent quadrant to the
      # sources. If you change the image size (imsize), then you must change the
      # statistics regions to match.

      ia.open(resname)
      stats0 = ia.statistics(logfile=dirname+'/all_stats.txt')
      qq = rg.box(blc=[1500,100,0,0],trc=[1600,200,0,0])
      stats1 = ia.statistics(region=qq,logfile=dirname+'/off_source_stats.txt')
      qq = rg.box(blc=[1030,1040,0,0],trc=[1050,1060,0,0])
      stats2 = ia.statistics(region=qq,logfile=dirname+'/near_source_stats.txt')
      ia.close()


  return [stats0,stats1,stats2]

###############################################

def makeBeamDifference(pb1="primary-beam-model",pb2="primary-beam-perturbed"):

    ia.open(pb1)
    pix1 = ia.getchunk()
    ia.close()

    ia.open(pb2)
    pix2 = ia.getchunk()
    ia.close()

    os.system('cp -r '+pb1+' pbdiff')

    ia.open("pbdiff")
    ia.putchunk(pix1-pix2)
    ia.close()

###############################################


def makeMSFrame(dirname,msname,ra0,dec0,nchan,numants):
  msstokes='XX YY';
  feedtype='perfect X Y';


  ## Directory for the MS
  if(not os.path.exists(dirname)):
    cmd = 'mkdir ' + dirname;
    os.system(cmd);

  #vx = [41.1100006,  -34.110001,  -268.309998,  439.410004,  -444.210022]
  #vy = [3.51999998, 129.8300018,  +102.480003, -182.149994, -277.589996]
  #vz = [0.25,       -0.439999998, -1.46000004, -3.77999997, -5.9000001]
  #d = [12.0,       12.0,         12.0,         12.0,       12.0]
  #an = ['ALMA1','ALMA2','ALMA3','ALMA4','ALMA5'];
  #nn = len(vx)*2.0;
  #x = (vx - (sum(pl.array(vx))/(nn)));
  #y = (vy - (sum(pl.array(vy))/(nn)));
  #z = (vz - (sum(pl.array(vz))/(nn)));

  ####  This call will get locations for all 27 vla antennas.
  #d, an, x, y, z = getAntLocations()

  ####  This call will get locations for all ALMA antennas.
  #x = zeros(numants, 'float')
  #y = zeros(numants, 'float')
  #z = zeros(numants, 'float')
  xx,yy,zz,d,pnames,nant,telescopename = util.readantenna('alma_config_files/alma.cycle2.5.cfg')
  an = pnames[0:numants]
  d = d[0:numants]
  nn = len(xx);
  #i = 0
  #while i < numants:
  #    x[i] = (xx[i] - (sum(pl.array(xx[i]))/(nn)));
  #    y[i] = (yy[i] - (sum(pl.array(yy[i]))/(nn)));
  #    z[i] = (zz[i] - (sum(pl.array(zz[i]))/(nn)));
  #    i = i + 1
  #x = (xx - (sum(pl.array(xx))/(nn)));
  #y = (yy - (sum(pl.array(yy))/(nn)));
  #z = (zz - (sum(pl.array(zz))/(nn)));

  pl.plot(xx,yy,"o")
  pl.savefig('ants')

  x = xx[0:numants]
  y = yy[0:numants]
  z = zz[0:numants]


  obspos = me.observatory('ALMA');
  #obspos = me.position('ITRF', '-0.0m', '0.0m', '3553971.510m');

  ## Make MS Frame
  sm.open(ms=msname);
  sm.setconfig(telescopename='ALMA',x=x,y=y,z=z,dishdiameter=d,
               mount=['alt-az'], antname=an,
               coordsystem='global',referencelocation=obspos);
  #sm.setconfig(telescopename='ALMA',x=x.tolist(),y=y.tolist(),z=z.tolist(),dishdiameter=d,
  #             mount=['alt-az'], antname=an,
  #             coordsystem='global',referencelocation=obspos);
  sm.setspwindow(spwname="3Band",freq="100.0GHz",deltafreq='500MHz',
                 freqresolution='2MHz',nchannels=nchan,stokes=msstokes);
  sm.setfeed(mode=feedtype,pol=['']);
  sm.setfield( sourcename="fake",sourcedirection=me.direction(rf='J2000',v0=ra0,v1=dec0) );
  sm.setlimits(shadowlimit=0.01, elevationlimit='10deg');
  sm.setauto(autocorrwt=0.0);
  sm.settimes(integrationtime='1800s', usehourangle=True,
                       referencetime=me.epoch('UTC','2013/05/10/00:00:00'));
  # Every 30 minutes, from -1h to +1h
  ostep = 0.5
  for loop in pl.arange(-1.0,+1.0,ostep):
    starttime = loop
    stoptime = starttime + ostep
    print starttime, stoptime
    for ch in range(0,nchan):
        sm.observe(sourcename="fake",spwname='3Band',
                   starttime=str(starttime)+'h', stoptime=str(stoptime)+'h');
  sm.close();

  listobs(vis=msname)
  
  return d

###############################################

def imageFromArray(arr,outfile,coord={},linear=F):
    newia = casac.image()
    newia.fromarray(outfile=outfile,pixels=arr,csys=coord,linear=F,overwrite=T)
    newia.close()

###############################################

def makeAperture(image="model",imsize=256,cellsize='8.0arcsec',
                    reffreq='1.5GHz', d = 12.0, noise = 0.0, supports = True,
                    ell_u = 1.0, ell_v = 1.0,
                    pointing = True, rot = False):

        # noise == True: there is Gaussian noise in the aperture function
        # supports == True: there are shadows from the support beams

        # ell_u and ell_v describe how elliptical the aperture is. If both are
        # equal to each other, the aperture will be a circle. Any relative
        # change in ell_u and ell_v will result in an elliptical aperture.

        # theta = rotation of the beam, in degrees.

        cell = qa.convert(qa.quantity(cellsize),'arcmin')['value']
        fov = cell*imsize
        freq = qa.convert(qa.quantity(reffreq),'Hz')['value']
        c = 3e8
        j = sqrt(-1)
        wvlen = c / freq # observed wavelength
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
        aper = zeros((Nuv,Nuv))
        uu, vv = mgrid[:Nuv, :Nuv]
        circle = (ell_u * (uu - ((Nuv/2.0) - 0.5)) ** 2) + (ell_v * (vv -
                 ((Nuv/2.0) - 0.5)) ** 2)
        disk = circle < (d_uv/2.0)**2
        aper[disk] = 1.0
        for u in xrange(Nuv):
            for v in xrange(Nuv):
                if aper[u][v]==1.0:
                    if disk[u][v] == True:
                        if supports == True:
                            # to get secondary reflector support beam shadows
                            if rot == True:
                                if u == v:
                                    aper[u][v] = 0
                                elif v == Nuv - u - 1:
                                    aper[u][v] = 0
                                else:
                                    aper[u][v] = 1.0 + random.gauss(0,noise)
                            else:
                                if u == (Nuv/2) or v == (Nuv/2):
                                    aper[u][v] = 0
                                elif u == (Nuv/2 - 1) or v == (Nuv/2 - 1):
                                    aper[u][v] = 0
                                else:
                                    aper[u][v] = 1.0 + random.gauss(0,noise);
                        else:
                            aper[u][v] = 1.0 + random.gauss(0,noise);
       
        aper_r = zeros((Nuv,Nuv), 'complex')
        aper_l = zeros((Nuv,Nuv), 'complex')
        phs_r = zeros((Nuv,Nuv))
        phs_l = zeros((Nuv,Nuv))

        offset_r = (0.06*wvlen/d)
        phs_ru = -1 * (uvals[uu] + uvcell/2.0) * (offset_r) * 2 * pi
        phs_rv = -1 * (vvals[vv] + uvcell/2.0) * (offset_r) * 2 * pi
        phs_r = exp(j * (phs_ru + phs_rv))
        phs_l = exp(-1 * j * (phs_ru + phs_rv))
        
        aper_r = aper * phs_r
        aper_l = aper * phs_l

        #ia.open('reaperture_pol9.im')
        #aper_r_real = ia.getchunk()
        #ia.close()

        #ia.open('imaperture_pol9.im')
        #aper_r_imag = ia.getchunk()
        #ia.close()

        #ia.open('reaperture_pol12.im')
        #aper_l_real = ia.getchunk()
        #ia.close()

        #ia.open('imaperture_pol12.im')
        #aper_l_imag = ia.getchunk()
        #ia.close()

        #real_aper_r = zeros((1024,1024))
        #imag_aper_r = zeros((1024,1024))
        #real_aper_l = zeros((1024,1024))
        #imag_aper_l = zeros((1024,1024))

        #real_aper_r = aper_r_real[:,:,0,0]
        #imag_aper_r = aper_r_imag[:,:,0,0]
        #real_aper_l = aper_l_real[:,:,0,0]
        #imag_aper_l = aper_l_imag[:,:,0,0]

        #aper_r[1536:2560,1536:2560] = real_aper_r + j * imag_aper_r
        #aper_l[1536:2560,1536:2560] = real_aper_l + j * imag_aper_l
        
        imageFromArray(imag(aper_r),"squint_x")
        imageFromArray(imag(aper_l),"squint_y")

        # Add phase ramp to aperture function
        if pointing == True:
            offset_u = str(random.uniform(-0.5,0.5))+'arcmin'
            offset_v = str(random.uniform(-0.5,0.5))+'arcmin'
        else:
            offset_u = '0.0arcmin'
            offset_v = '0.0arcmin'
        phs_off_u = qa.quantity(offset_u)['value']
        phs_off_v = qa.quantity(offset_v)['value']
        shift_u = phs_off_u # in arcmins
        shift_v = phs_off_v # in arcmins
        phs = zeros((Nuv,Nuv))
        phs_u = -1 * (uvals[uu] + uvcell/2.0) * ( shift_u / 60.0 * ( pi / 180.0 ) ) * 2 * pi
        phs_v = -1 * (vvals[vv] + uvcell/2.0) * ( shift_v / 60.0 * ( pi / 180.0 ) ) * 2 * pi

        #real_noise = zeros((Nuv,Nuv))
        #imag_noise = zeros((Nuv,Nuv))
        #for i in xrange(Nuv):
        #    for j in xrange(Nuv):
        #        real_noise[i][j] = random.gauss(1,noise)
        #        imag_noise[i][j] = random.gauss(1,noise)
        #print aper_r.shape, aper_l.shape
        #aper_r = aper_r * abs(real_noise)
        #aper_l = aper_l * abs(real_noise)
        #phs_u = phs_u * imag_noise
        #phs_v = phs_v * imag_noise
        phs = exp(j*(phs_u + phs_v))

        phs_aper_r = aper_r * phs
        phs_aper_l = aper_l * phs

        imageFromArray(real(phs_aper_r),image+"X")
        imageFromArray(real(phs_aper_l),image+"Y")

        volt_r = fftpack.ifftshift(fftpack.fft2(fftpack.fftshift(phs_aper_r)))
        volt_l = fftpack.ifftshift(fftpack.fft2(fftpack.fftshift(phs_aper_l)))
        
        imageFromArray(real(volt_r),image+"Xreal")
        imageFromArray(imag(volt_r),image+"Ximag")
        imageFromArray(real(volt_l),image+"Yreal")
        imageFromArray(imag(volt_l),image+"Yimag")

        del aper_r,aper_l,phs,phs_u,phs_v,volt_r,volt_l,\
            phs_aper_r,phs_aper_l,uu, vv, uvals, vvals

        return Nxy

###############################################

def makePrimaryBeam(imsize=256,cellsize='8.0arcsec',coord="coord.torecord()",
                    reffreq='1.5GHz', pbname = "model",
                    aper1_Rreal = "aper00Rreal", aper1_Rimag = "aper00Rimag",
                    aper2_Rreal = "aper01Rreal", aper2_Rimag = "aper01Rimag",
                    aper1_Lreal = "aper00Lreal", aper1_Limag = "aper00Limag",
                    aper2_Lreal = "aper01Lreal", aper2_Limag = "aper01Limag",
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
            #auto_corr = signal.fftconvolve(aper1, aper2, 'same')
            time2 = time.time()
            print "convolution time = " + str(time2-time1)

            # Power pattern function
            #power = (Nxy**2) * fftpack.ifft2(fftpack.fftshift(auto_corr))
            #time3 = time.time()
            #print "power time = " + str(time3-time2)
            #power = fftpack.ifftshift(power)
            power = power / max(power)

            if i == 0:
                power_r = power
            else:
                power_l = power

            i = i + 1

        stokes_power = zeros((imsize,imsize,4,1))
        stokes_power[:,:,0,0] = (power_r[:,:] + power_l[:,:]) / 2
        stokes_power[:,:,1,0] = (power_r[:,:] - power_l[:,:]) / 2

        imageFromArray(real(stokes_power),pbname+"real",coord)
        imageFromArray(imag(stokes_power),pbname+"imag",coord)

        aper_area = -1

        del aper1_real, aper1_imag, aper2_real, aper2_imag, power,\
            power_r, power_l, stokes_power
        
###############################################

def makeResidualImage(msname='',resname='',imsize=256,cellsize='8.0arcsec',
                      ra0='', dec0='', nchan=1, reffreq='1.5GHz'):
  ## Make model image
    im.open(msname);
    im.selectvis()
    im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,
                      stokes='IQUV',spw=[0],
                      phasecenter=me.direction(rf='J2000',v0=ra0,v1=dec0),
                      mode='channel',nchan=nchan,start=0,step=1,
                      restfreq=reffreq);
    im.makeimage(type="residual",image=resname);
    im.done();

###############################################

def makeImage(msname='',imname='',clname='mysource.cl',
              imsize=256,cellsize='8.0arcsec',ra0='',dec0='',
              nchan=1,reffreq='1.5GHz'):

  ## Make model image
  os.system('rm -rf '+imname+'*');
  im.open(msname);
  im.selectvis()
  im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,
                 stokes='IQUV',spw=[0],
                 phasecenter=me.direction(rf='J2000',v0=ra0,v1=dec0),
                 mode='channel',nchan=nchan,start=0,step=1,
                 restfreq=reffreq);
  im.make(image=imname);
  im.done();

  cl.open(clname)
  ia.open(imname);
  ia.modify(model=cl.torecord(),subtract=False);
  cl.close();

###############################################

def makeTrueImage(stokesvals=[1.0,0.0,0.0,1.0],imname='',newimname='',
                   pb_real_file = "pb00&&01real",
                   pb_imag_file = "pb00&&01imag",
                   msname='', ftm='ft',imsize=256,predict_imsize=256,
                   cellsize='8.0arcsec',ra0='', dec0='',
                   nchan=1, reffreq='1.5GHz',
                   num_ant = 28, pair="00&&01",model = True):

  for scanid in range(4):
      print "scan #" + str(scanid)
      time1 = time.time()
      cmd = 'cp -r '+imname+' '+newimname
      cmd = cmd.replace('&','\&')
      os.system(cmd)

      ia.open(pb_real_file)
      pb_real = ia.getchunk()
      ia.close()

      ia.open(pb_imag_file)
      pb_imag = ia.getchunk()
      ia.close()

      time2 = time.time()
      print "open files: " + str(time2-time1)
        
      tb.open(msname)
      tb1 = tb.query('SCAN_NUMBER=='+str(scanid+1))
      timelist = tb1.getcol('TIME')
      print len(timelist)
      tb1.close()
      tb.close()
      theta = parang(msname,0,timelist[0])
      if scanid == 0:
          start_theta = theta
      theta = theta - start_theta
      print "time = " + str(timelist[0])
      print "theta = " + str(theta)

      time3 = time.time()
      low = imsize / 2 - predict_imsize / 2 + 1
      high = imsize / 2 + predict_imsize / 2
      pb_real=ndimage.rotate(input=pb_real[low:high,low:high,:,:],angle=theta,reshape=False)
      pb_imag=ndimage.rotate(input=pb_imag[low:high,low:high,:,:],angle=theta,reshape=False)
      time4 = time.time()
      print "rotation: " + str(time4-time3)

      pb = pb_real + sqrt(-1) * pb_imag

      ia.open(imname)
      vals = ia.getchunk()
      ia.close()

      # Fill in from the componentlist 
      ia.open(newimname)
      vals = vals[1:predict_imsize,1:predict_imsize,:,:] * pb_real
      ia.putchunk(vals)
      ia.close();
      time5 = time.time()
      print "new image: " + str(time5-time4)

      ### Predicting
      im.open(msname,usescratch=True);
      im.selectvis(baseline=pair,nchan=nchan,start=0,step=1,scan=str(scanid+1));
      im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,
                     stokes='IQUV',spw=[0],
                     phasecenter=me.direction(rf='J2000',v0=ra0,v1=dec0),
                     mode='channel',nchan=nchan,start=0,step=1,
                     restfreq=reffreq);
      if(ftm=="awp"):
         im.setoptions(cache=imsize*imsize*6*nchan,ftmachine=ftm,
                       applypointingoffsets=False,
                       dopbgriddingcorrections=False,
                       cfcachedirname=imname+".cfcache",
                       pastep=360.0,pblimit=0.001);
      else:
         im.setoptions(ftmachine="ft");
      #im.setvp(dovp=True,usedefaultvp=True,telescope='EVLA');
      im.ft(model=newimname,incremental=False);
      im.done();
      time6 = time.time()
      print "prediction: " + str(time6-time5)

#####################################

def makeDataTable(msname=''):
  
  ### Copy to the data and corrected-data columns
  tb.open(msname,nomodify=False);
  moddata = tb.getcol(columnname='MODEL_DATA');
  #if model == True:
  tb.putcol(columnname='DATA',value=moddata);
  tb.putcol(columnname='CORRECTED_DATA',value=moddata);
  moddata.fill(0.0);
  tb.putcol(columnname='MODEL_DATA',value=moddata);
  tb.close();

#####################################

def addNoise(msname=''):
  sm.openfromms(msname);
  sm.setnoise(mode='simplenoise',simplenoise='0.3Jy');
  sm.corrupt();
  sm.close();

###############################################

def getAntLocations(configuration='D'):
   vx = [41.1100006,  134.110001,   268.309998,  439.410004,  644.210022,
         880.309998,  1147.10999,  1442.41003,  1765.41003,  -36.7900009,
         -121.690002,  -244.789993, -401.190002, -588.48999,  -804.690002,
         -1048.48999, -1318.48999, -1613.98999,  -4.38999987,-11.29,  
         -22.7900009, -37.6899986, -55.3899994, -75.8899994, -99.0899963, 
         -124.690002, -152.690002];
   vy = [3.51999998, -39.8300018,  -102.480003, -182.149994, -277.589996,
        -387.839996, -512.119995, -649.76001,  -800.450012, -2.58999991,
        -59.9099998,  -142.889999, -248.410004, -374.690002, -520.599976,
        -685,        -867.099976, -1066.42004,   77.1500015,  156.910004, 
        287.980011,  457.429993,  660.409973,  894.700012,  1158.82996, 
        1451.43005,  1771.48999];
   vz = [0.25,       -0.439999998, -1.46000004, -3.77999997, -5.9000001,
        -7.28999996, -8.48999977, -10.5,       -9.56000042,      0.25, 
        -0.699999988, -1.79999995, -3.28999996, -4.78999996, -6.48999977,
        -9.17000008, -12.5299997, -15.3699999,  1.25999999,   2.42000008, 
        4.23000002,  6.65999985,  9.5,         12.7700005,  16.6800003, 
        21.2299995,  26.3299999];
   d = [25.0,       25.0,         25.0,         25.0,       25.0,
        25.0,       25.0,         25.0,         25.0,       25.0,
        25.0,       25.0,         25.0,         25.0,       25.0,
        25.0,       25.0,         25.0,         25.0,       25.0,
        25.0,       25.0,         25.0,         25.0,       25.0,
        25.0,       25.0];
   if(configuration=='D'):
      scale = 3.0;
   else:
     if(configuration=='A'):
       scale = 1.0/9.0;
     else:
       scale = 1.0;
       print 'Using VLA C-array coords';
   nn = len(vx);
   x = (vx - (sum(pl.array(vx))/(nn)))/scale;
   y = (vy - (sum(pl.array(vy))/(nn)))/scale;
   z = (vz - (sum(pl.array(vz))/(nn)))/scale;
   an=[];
   for i in range(0,nn):
           an.append("VLA"+str(i));
   return d,an,x,y,z




###############
####  Imaging
###############
