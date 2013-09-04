from numpy import *
from scipy import *
import random

def makeComponentList(clname="sim.cl",reffreq="1.42GHz",
		      ra0="19:59:28.500",dec0="+40.44.01.50",
		      ncomps=1,shapes=["Gaussian"],
		      fluxvals=[1.0],minaxes=["1.0arcmin"],
		      majaxes=["1.0arcmin"],ras=["5.0arcmin"],
		      decs=["5.0arcmin"],spxs=[-0.5],posangles=["0.0deg"]):
  """ Make a component-list
      params :
  """
  if(ncomps != len(fluxvals) or ncomps != len(minaxes) or  
     ncomps != len(majaxes) or  ncomps != len(ras) or  
     ncomps != len(decs) or  ncomps != len(spxs) or  
     ncomps != len(posangles) or  ncomps != len(shapes) ):
	  print "Please enter all params for each component"
	  return;
  os.system('rm -rf '+clname)
  refRA = qa.unit(ra0);
  refDEC = qa.unit(dec0);
  for comp in range(0,ncomps):
    if(shapes[comp]=="point"):
      cl.addcomponent(flux=fluxvals[comp],fluxunit="Jy",
                  dir=me.direction(rf='J2000',
                                   v0=qa.add(refRA,ras[comp]),
                                   v1=qa.add(refDEC,decs[comp])),
                  shape="point",freq=reffreq,
                  spectrumtype="spectral index",index=spxs[comp]); #[spxs[comp],0,0,0]);
    else:
      cl.addcomponent(flux=fluxvals[comp],fluxunit="Jy",
                  dir=me.direction(rf='J2000',
                                   v0=qa.add(refRA,ras[comp]),
                                   v1=qa.add(refDEC,decs[comp])),
                  shape="Gaussian",freq=reffreq,
                  majoraxis=majaxes[comp],minoraxis=minaxes[comp],
                  positionangle=posangles[comp],
                  spectrumtype="spectral index",index=spxs[comp]); #[spxs[comp],0,0,0]);
  cl.rename(filename=clname);
  cl.close();


###################################################

## Define the locations of your sources ( check that these locations are within the image
### field of you that you are simulating ! )

clname = 'mysources';
ra0="19:59:28.500";
dec0="+40.44.01.50";
reffreq = '12.0GHz'

#ncomps = 1000;
#shapes = ["point" for x in range(1000)]
#fluxvals = ones(1000,float)
#minaxes = empty(1000)
#minaxes = ["" for x in range(1000)]
#majaxes = ["" for x in range(1000)]
#posangles = ["" for x in range(1000)]
#ras = ["" for x in range(1000)]
#decs = ["" for x in range(1000)]
#spxs = zeros(1000,float)
#random.seed(1)
#for i in arange(1000):
#    ra = 100*random.random()
#    dec = 100*random.random()
#    ra = str(ra)+"arcmin"
#    dec = str(dec)+"arcmin"
#    ras[i] = ra
#    decs[i] = dec

## Make a list of components.
# components placed such that there is one source per simulation, placed at
# full-power, half-power, and around 10-20% power points of the beam.
makeComponentList(clname=clname+'0.cl',ncomps=1,shapes=['point'],fluxvals=[1.0],
        minaxes=[""], majaxes=[""], posangles=[""],
        ras=["0.0arcmin"],decs=["0.0arcmin"],
        spxs=[0.0],reffreq=reffreq);
makeComponentList(clname=clname+'1.cl',ncomps=3,
        shapes=['point','point','point'],fluxvals=[1.0,1.0,1.0],
        minaxes=["","",""], majaxes=["","",""],
        posangles=["","",""],
        ras=["0.0arcmin","1.22arcmin","2.02arcmin"],
        decs=["0.0arcmin","1.22arcmin","2.02arcmin"],
        spxs=[0.0,0.0,0.0],reffreq=reffreq);
#makeComponentList(clname=clname+'2.cl',ncomps=1,shapes=['point'],fluxvals=[1.0],
#        minaxes=[""], majaxes=[""],
#        posangles=[""],ras=["2.02arcmin"],decs=["2.02arcmin"],
#        spxs=[0.0],reffreq=reffreq);


## Un-comment and put these few lines into the  makeTrueImage() function - in the place where you've been modifying the source positions. This will evaluate the components onto the image grid. You should be able to see the Gaussians in the 'true' image after this. 

## Fill in from the componentlist 
#cl.open(clname);
#ia.open(imname);
#ia.modify(model=cl.torecord(),subtract=False);
#ia.close();
#cl.close();

###################################################

