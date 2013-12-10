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

import matplotlib.pylab as pl
import random
import time

sourcedir = '/lustre/kkundert/Code/'

execfile(sourcedir+'mksrc.py')
execfile(sourcedir+'parang.py')
execfile(sourcedir+'writeAps.py')
execfile(sourcedir+'makeBeams.py')

###############################################

def makeMS(runnum=0, makeBeams = True,
           noise=0.0, supports=True,
           ell_u = 1.0, ell_v = 1.0,
           pointing = False, theta = 0.0):


  for i in xrange(1):
      dirname = "Data"+str(runnum)
      if i == 0:
          dirname = dirname+"-centered"
      elif i == 1:
          dirname = dirname+"-half-power"
      else:
          dirname = dirname+"-low-power"
      basename = dirname+"/points"
      msname = basename + '.ms';
      imname = basename+'.true.im';
      predict_imname = imname + '.predict'
      resname = dirname+"/pb-residuals-"+str(runnum)
      pbname = dirname+"/primary-beams"
      apname = dirname+"/apertures"
      clname = "mysources"
      ra0="19:59:28.500";
      dec0="+40.44.01.50";
      nchan=1;
      imsize=4096;
      predict_imsize = 1024
      cellsize='3.8arcsec';
      reffreq='6.0GHz';
      stokesvals=[1.0,0.0,0.0,1.0]
      ftm='ft'

      num_ant = 5

      #addNoise(msname);
      if makeBeams == True:
          os.system('rm -rf '+dirname)
          os.system('rm -rf '+resname)
          os.system('rm -rf theresult.*')

          clname = clname+str(i)+'.cl'
          d = makeMSFrame(dirname=dirname,msname=msname,
                      ra0=ra0,dec0=dec0,nchan=nchan);

          os.system('mkdir '+apname)
          os.system('mkdir '+pbname)

          makeImage(msname=msname,imname=imname,clname=clname,
                    imsize=imsize,cellsize=cellsize,ra0=ra0,dec0=dec0,
                    nchan=nchan,reffreq=reffreq)

          ia.open(imname)
          coord = ia.coordsys()
          ia.close()

          os.system('rm -rf ap.ms')
          os.system('cp -r '+msname+' ap.ms')
          writeApertures(msname='ap.ms',imname=imname)
          print "apertures written"
          
          Nxy_list = []
          for i in xrange(num_ant):
              print "aper" + str(i)
              image = apname+"/aper%02d" % i
              Nxy = makeAperture(image=image,imsize=imsize,
                              cellsize=cellsize,reffreq=reffreq,
                              d=d[i],noise=noise,supports=supports,
                              ell_u=ell_u,ell_v=ell_v,
                              pointing=pointing)
              Nxy_list.append(Nxy)
          for i in xrange(num_ant):
              for j in xrange(num_ant):
                  if i < j:
                      pbimage = pbname+"/pb%02d&&%02d" % (i,j)
                      print "primary beam"+str(i)+"&&"+str(j)
                      makePrimaryBeam_VLA(imsize=imsize,cellsize=cellsize,
                                    coord = coord.torecord(),
                                    reffreq=reffreq,pbname=pbimage,
                                    aper1_Rreal=apname+"/aper%02dRreal" % i,
                                    aper1_Rimag=apname+"/aper%02dRimag" % i,
                                    aper2_Rreal=apname+"/aper%02dRreal" % j,
                                    aper2_Rimag=apname+"/aper%02dRimag" % j,
                                    aper1_Lreal=apname+"/aper%02dLreal" % i,
                                    aper1_Limag=apname+"/aper%02dLimag" % i,
                                    aper2_Lreal=apname+"/aper%02dLreal" % j,
                                    aper2_Limag=apname+"/aper%02dLimag" % j,
                                    Nxy=Nxy_list[i]);
                  else:
                      continue
      else:
          clname = clname+str(i)+'.cl'
          makeMSFrame(dirname=dirname,msname=msname,
                      ra0=ra0,dec0=dec0,nchan=nchan);
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
      clname = 'mysources'

      # The sources are located in one quadrant of the sky, and therefore one
      # quadrant of the primary beam. The statistics are chosen such that off
      # source is on the far side of the opposite quadrant from the sources, and
      # the near source stats are on the axis in the adjacent quadrant to the
      # sources. If you change the image size (imsize), then you must change the
      # statistics regions to match.

      #ia.open(resname)
      #stats0 = ia.statistics(logfile=dirname+'/all_stats.txt')
      #qq = rg.box(blc=[1500,100,0,0],trc=[1600,200,0,0])
      #stats1 = ia.statistics(region=qq,logfile=dirname+'/off_source_stats.txt')
      #qq = rg.box(blc=[1030,1040,0,0],trc=[1050,1060,0,0])
      #stats2 = ia.statistics(region=qq,logfile=dirname+'/near_source_stats.txt')
      #ia.close()


  ##return [stats0,stats1,stats2]

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


def makeMSFrame(dirname,msname,ra0,dec0,nchan):
  msstokes='RR LL';
  feedtype='perfect R L';


  ## Directory for the MS
  if(not os.path.exists(dirname)):
    cmd = 'mkdir ' + dirname;
    os.system(cmd);

  vx = [41.1100006,  -34.110001,  -268.309998,  439.410004,  -444.210022]
  vy = [3.51999998, 129.8300018,  +102.480003, -182.149994, -277.589996]
  vz = [0.25,       -0.439999998, -1.46000004, -3.77999997, -5.9000001]
  d = [25.0,       25.0,         25.0,         25.0,       25.0]
  an = ['VLA1','VLA2','VLA3','VLA4','VLA5'];
  nn = len(vx)*2.0;
  x = 0.5*(vx - (sum(pl.array(vx))/(nn)));
  y = 0.5*(vy - (sum(pl.array(vy))/(nn)));
  z = 0.5*(vz - (sum(pl.array(vz))/(nn)));

  ####  This call will get locations for all 27 vla antennas.
  #d, an, x, y, z = getAntLocations()


  obspos = me.observatory('EVLA');
  #obspos = me.position('ITRF', '-0.0m', '0.0m', '3553971.510m');

  ## Make MS Frame
  sm.open(ms=msname);
  sm.setconfig(telescopename='EVLA',x=x.tolist(),y=y.tolist(),z=z.tolist(),dishdiameter=d,
               mount=['alt-az'], antname=an,
               coordsystem='local',referencelocation=obspos);
  sm.setspwindow(spwname="CBand",freq="6.0GHz",deltafreq='500MHz',
                 freqresolution='2MHz',nchannels=nchan,stokes=msstokes);
  sm.setfeed(mode=feedtype,pol=['']);
  sm.setfield( sourcename="fake",sourcedirection=me.direction(rf='J2000',v0=ra0,v1=dec0) );
  sm.setlimits(shadowlimit=0.01, elevationlimit='10deg');
  sm.setauto(autocorrwt=0.0);
  sm.settimes(integrationtime='1800s', usehourangle=True,
                       referencetime=me.epoch('UTC','2013/05/10/00:00:00'));
  # Every 30 minutes, from -3h to +3h
  ostep = 0.5
  for loop in pl.arange(-3.0,+3.0,ostep):
    starttime = loop
    stoptime = starttime + ostep
    print starttime, stoptime
    for ch in range(0,nchan):
        sm.observe(sourcename="fake",spwname='CBand',
                   starttime=str(starttime)+'h', stoptime=str(stoptime)+'h');
  sm.close();

  listobs(vis=msname)

  return d

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

  for scanid in range(12):
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
