###
### Simulate Sky and Data...
###

from numpy import *
from scipy import *
from scipy import fftpack
from scipy import signal
from numpy import fft
from scipy import ndimage

import matplotlib.pylab as pl
import random

###############################################

def makeMS(runnum=0, noise=0.0, supports=True,
           offset_u = '0.0arcmin', offset_v = '0.0arcmin',
           ell_u = 1.0, ell_v = 1.0, theta = 0.0):


  for i in xrange(3):
      dirname = "Data"+str(runnum)
      basename = dirname+"/points"
      if i == 0:
          basename = basename+"-centered"
      elif i == 1:
          basename = basename+"-half-power"
      else:
          basename = basename+"-low-power"
      msname = basename + '.ms';
      imname = basename+'.true.im';
      resname = dirname+"/pb-residuals-"+str(runnum)
      pbname = dirname+"/primary-beam"
      clname = "mysources"
      ra0="19:59:28.500";
      dec0="+40.44.01.50";
      nchan=1;
      imsize=2048;
      cellsize='0.07arcmin';
      reffreq='6.0GHz';
      stokesvals=[1.0,0.0,0.0,0.0]
      ftm='ft'

      os.system('rm -rf '+dirname)
      os.system('rm -rf '+resname)
      os.system('rm -rf theresult.*')

      clname = clname+str(i)+'.cl'
      makeMSFrame(dirname=dirname,msname=msname,ra0=ra0,dec0=dec0,nchan=nchan);
      addNoise(msname);
      area = makeTrueImage(stokesvals=stokesvals,msname=msname,imname=imname,
                    pbname=pbname+"-model",clname=clname,imsize=imsize,cellsize=cellsize,
                    ra0=ra0, dec0=dec0, nchan=nchan, reffreq=reffreq,
                    noise=0.0, supports=True,offset_u='0.0arcmin',
                    offset_v='0.0arcmin',ell_u=1.0,ell_v=1.0,theta=0.0);
      predictTrueImage(msname=msname,ftm=ftm,imname=imname,
                       imsize=imsize,cellsize=cellsize,ra0=ra0, dec0=dec0,
                       nchan=nchan, reffreq=reffreq, model=True);
      makeTrueImage(stokesvals=stokesvals,msname=msname,imname=imname,
                    pbname=pbname+"-perturbed",clname=clname,imsize=imsize,cellsize=cellsize,
                    ra0=ra0, dec0=dec0, nchan=nchan, reffreq=reffreq,
                    noise=noise, supports=supports,offset_u=offset_u,
                    offset_v=offset_v,ell_u=ell_u,ell_v=ell_v,theta=theta,
                    area=area);
      predictTrueImage(msname=msname,ftm=ftm,imname=imname,
                       imsize=imsize,cellsize=cellsize,ra0=ra0, dec0=dec0,
                       nchan=nchan, reffreq=reffreq, model=False);
      makeResidualImage(msname,resname,imsize,cellsize,ra0, dec0, nchan, reffreq);
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
  qq = rg.box(blc=[1030,1030,0,0],trc=[1040,1040,0,0])
  stats2 = ia.statistics(region=qq,logfile=dirname+'/near_source_stats.txt')
  ia.close()

  return [stats0,stats1,stats2]

###############################################


def makeMSFrame(dirname,msname,ra0,dec0,nchan):
  msstokes='RR LL';
  feedtype='perfect R L';


  ## Directory for the MS
  if(not os.path.exists(dirname)):
    cmd = 'mkdir ' + dirname;
    os.system(cmd);

  #vx = [41.1100006,  -34.110001,  -268.309998,  439.410004,  -444.210022]
  #vy = [3.51999998, 129.8300018,  +102.480003, -182.149994, -277.589996]
  #vz = [0.25,       -0.439999998, -1.46000004, -3.77999997, -5.9000001]
  #d = [25.0,       25.0,         25.0,         25.0,       25.0]
  #an = ['VLA1','VLA2','VLA3','VLA4','VLA5'];
  #nn = len(vx)*2.0;
  #x = 5.0*(vx - (sum(pl.array(vx))/(nn)));
  #y = 5.0*(vy - (sum(pl.array(vy))/(nn)));
  #z = 5.0*(vz - (sum(pl.array(vz))/(nn)));

  ####  This call will get locations for all 27 vla antennas.
  d, an, x, y, z = getAntLocations()


  obspos = me.observatory('VLA');
  #obspos = me.position('ITRF', '-0.0m', '0.0m', '3553971.510m');

  ## Make MS Frame
  sm.open(ms=msname);
  sm.setconfig(telescopename='VLA',x=x.tolist(),y=y.tolist(),z=z.tolist(),dishdiameter=d,
               mount=['alt-az'], antname=an,
               coordsystem='local',referencelocation=obspos);
  sm.setspwindow(spwname="LBand",freq="1.0GHz",deltafreq='500MHz',
                 freqresolution='2MHz',nchannels=nchan,stokes=msstokes);
  sm.setfeed(mode=feedtype,pol=['']);
  sm.setfield( sourcename="fake",sourcedirection=me.direction(rf='J2000',v0=ra0,v1=dec0) );
  sm.setlimits(shadowlimit=0.01, elevationlimit='10deg');
  sm.setauto(autocorrwt=0.0);
  sm.settimes(integrationtime='1000s', usehourangle=True,
                       referencetime=me.epoch('UTC','2013/05/10/00:00:00'));
  sm.observe(sourcename="fake",spwname='LBand', starttime='-3.0h', stoptime='+2.0h');
  sm.close();

###############################################

def makePrimaryBeam(imsize=256,cellsize='8.0arcsec',reffreq='1.5GHz', noise =
                    0.0, supports=True, offset_u='0.0arcmin',
                    offset_v='0.0arcmin', ell_u = 1.0, ell_v = 1.0,
                    theta = 0.0, area = -1):

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
        aper1 = zeros((Nuv,Nuv))
        d_uv = spat_lam / uvcell
        uu, vv = mgrid[:Nuv, :Nuv]
        circle = (ell_u * (uu - ((Nuv/2.0) - 0.5)) ** 2) + (ell_v * (vv -
                 ((Nuv/2.0) - 0.5)) ** 2)
        disk = circle < 4*(d_uv/2.0)**2
        aper[disk] = 1.0
        for u in xrange(Nuv):
            for v in xrange(Nuv):
                if aper[u][v]==1.0:
                    if disk[u][v] == True:
                        if supports == True:
                            # to get secondary reflector support beam shadows
                            if u == (Nuv/2) or v == (Nuv/2):
                                aper1[u][v] = 0
                            elif u == (Nuv/2 - 1) or v == (Nuv/2 - 1):
                                aper1[u][v] = 0
                            else:
                                aper1[u][v] = 1.0 + random.gauss(0,noise);
                        else:
                            aper1[u][v] = 1.0 + random.gauss(0,noise);

        aper=ndimage.rotate(input=aper1,angle=theta,reshape=False)

        pl.figure(1)
        pl.clf()
        pl.imshow(aper)
        
        # Add phase ramp to aperture function
        phs_off_u = qa.quantity(offset_u)['value']
        phs_off_v = qa.quantity(offset_v)['value']
        shift_u = phs_off_u # in arcmins
        shift_v = phs_off_v # in arcmins
        phs = zeros((Nuv,Nuv))
        uu,vv = mgrid[:Nuv, :Nuv]
        j = sqrt(-1)
        phs_u = -1 * (uvals[uu] + uvcell/2.0) * ( shift_u / 60.0 * ( pi / 180.0 ) ) * 2 * pi
        phs_v = -1 * (vvals[vv] + uvcell/2.0) * ( shift_v / 60.0 * ( pi / 180.0 ) ) * 2 * pi
        phs = exp(j*(phs_u + phs_v))
        if noise == True:
            imag_noise = random.gauss(1,0.3)
            phs_u = phs_u * imag_noise
            phs_v = phs_v * imag_noise

        phs_aper = aper * phs

        #pl.figure(4)
        #pl.clf()
        #pl.subplot(121)
        #pl.imshow(imag(phs_aper))
        #pl.colorbar()
        #pl.subplot(122)
        #pl.imshow(real(phs_aper))
        #pl.colorbar()
        
        # Aperture autocorrelation
        auto_corr = signal.fftconvolve(phs_aper, phs_aper, 'same')
        aper_area = auto_corr.sum()
        print aper_area
        print area
        if area != -1:
            auto_corr = auto_corr / area
        else:
            auto_corr = auto_corr / aper_area

        # Power pattern function
        power = (Nxy**2) * fftpack.ifft2(fftpack.fftshift(auto_corr))
        power = fftpack.ifftshift(power)

        pl.figure(2)
        pl.clf()
        pl.imshow(real(power))
        pl.colorbar()
        #pl.figure(3)
        #pl.clf()
        #pl.imshow(imag(power))
        #pl.colorbar()
        
        return power,aper_area

###############################################

def makeResidualImage(msname='',resname='',imsize=256,cellsize='8.0arcsec',
                      ra0='', dec0='', nchan=1, reffreq='1.5GHz'):
  ## Make model image
    im.open(msname);
    im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,
                      stokes='IQUV',spw=[0],
                      phasecenter=me.direction(rf='J2000',v0=ra0,v1=dec0),
                      mode='channel',nchan=nchan,start=0,step=1,
                      restfreq=reffreq);
    im.makeimage(type="residual",image=resname);
    im.close();

###############################################

def makeTrueImage(stokesvals=[1.0,0.0,0.0,0.0],msname='',
                   imname='',pbname='',clname='mysources.cl',
                   imsize=256,cellsize='8.0arcsec',
                   ra0='', dec0='', nchan=1, reffreq='1.5GHz',
                   noise=0.0, supports=True,
                   offset_u = '0.0arcmin', offset_v = '0.0arcmin',
                   ell_u = 1.0, ell_v = 1.0, theta = 0.0, area=-1):

    # noise == True: there is Gaussian noise in the aperture function
    # supports == True: there are shadows from the support beams

  ## Make model image
  os.system('rm -rf '+imname);
  im.open(msname);
  im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,
                 stokes='IQUV',spw=[0],
                 phasecenter=me.direction(rf='J2000',v0=ra0,v1=dec0),
                 mode='channel',nchan=nchan,start=0,step=1,
                 restfreq=reffreq);
  im.make(image=imname);
  im.make(image=pbname)
  im.close();

  # Fill in from the componentlist 
  cl.open(clname);
  ia.open(imname);
  ia.modify(model=cl.torecord(),subtract=False);
  cl.close();
  vals = ia.getchunk();
  pb,aper_area = makePrimaryBeam(imsize,cellsize,reffreq, noise, supports, offset_u,
                       offset_v, ell_u, ell_v, theta, area);
  vals[:,:,0,0] = vals[:,:,0,0] * real(pb[:,:]);
  ia.putchunk(vals);
  ia.close();

  ia.open(pbname)
  vals = ia.getchunk()
  vals[:,:,0,0] = real(pb[:,:])
  ia.putchunk(vals)
  ia.close()

  # Fill in model image with 'stokesvals'
  #ia.open(imname);  
  #vals = ia.getchunk();
  #vals.fill(0.0);
  #shp = vals.shape;
  #for stokes in range(0,shp[2]):
  #  for chan in range(0,shp[3]):
  #    xpix = shp[0]/2.0+5;
  #    ypix = shp[0]/2.0;
  #    vals[xpix,ypix,stokes,chan] = stokesvals[stokes];
  #    xpix = shp[0]/2.0;
  #    vals[xpix,ypix,stokes,chan] = stokesvals[stokes];
  #pb = makePrimaryBeam(imsize,cellsize,reffreq);
  #vals = vals[:,:,0,0] * abs(pb[:,:]);
  #ia.putchunk(vals);
  #ia.close();

  return aper_area


###############################################

def predictTrueImage(msname='', ftm='ft',imname='',imsize=256,
                                 cellsize='8.0arcsec',ra0='', dec0='',
                                 nchan=1, reffreq='1.5GHz', model = True):
  ### Predicting
  im.open(msname,usescratch=True);
  im.selectvis(nchan=nchan,start=0,step=1);
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
  im.ft(model=imname,incremental=False);
  im.close();

  ### Copy to the data and corrected-data columns
  tb.open(msname,nomodify=False);
  moddata = tb.getcol(columnname='MODEL_DATA');
  if model == True:
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
