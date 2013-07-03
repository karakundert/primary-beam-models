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
clname = 'mysources.cl';
ra0="19:59:28.500";
dec0="+40.44.01.50";
ncomps = 300;
#shapes = ["point","point","point","gaussian","gaussian"];
#fluxvals = [1.0,1.0,1.0,0.0,0.0];
#minaxes = ["","","","0.5arcmin","0.3arcmin"];
#majaxes = ["","","","0.5arcmin","0.5arcmin"];
#posangles = ["","","","0.0deg","45deg"];
#ras =  ["11.3arcmin","0.0arcmin","-2.7arcmin","0.0arcmin","-11.0arcmin"];
#decs = ["2.7arcmin","0.0arcmin","-11.3arcmin","0.0arcmin","5.0arcmin"];
#spxs = [0.0,0.0,0.0,0.0,0.0];
reffreq = '1.5GHz'

shapes = ["point" for x in range(300)]
fluxvals = ones(300,float)
minaxes = empty(300)
minaxes = ["" for x in range(300)]
majaxes = ["" for x in range(300)]
posangles = ["" for x in range(300)]
ras = ["" for x in range(300)]
decs = ["" for x in range(300)]
spxs = zeros(300,float)
random.seed(1)
for i in arange(300):
    ra = 40*random.random()
    dec = 40*random.random()
    ra = str(ra)+"arcmin"
    dec = str(dec)+"arcmin"
    ras[i] = ra
    decs[i] = dec

## Make a list of components.
makeComponentList(clname=clname,ncomps=ncomps,shapes=shapes,fluxvals=fluxvals, minaxes=minaxes, majaxes=majaxes, posangles=posangles,ras=ras,decs=decs,spxs=spxs,reffreq=reffreq);


## Un-comment and put these few lines into the  makeTrueImage() function - in the place where you've been modifying the source positions. This will evaluate the components onto the image grid. You should be able to see the Gaussians in the 'true' image after this. 

## Fill in from the componentlist 
#cl.open(clname);
#ia.open(imname);
#ia.modify(model=cl.torecord(),subtract=False);
#ia.close();
#cl.close();

###################################################

