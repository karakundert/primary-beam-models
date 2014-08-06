def makeComponentList(clname="sim.cl",reffreq="1.5GHz",
		      ra0="19:59:28.500",dec0="-23.44.01.50",
		      ncomps=1,shapes=["Gaussian"],
		      fluxvals=[1.0],minaxes=["1.0arcmin"],
		      majaxes=["1.0arcmin"],ras=["5.0arcmin"],
		      decs=["5.0arcmin"],spxs=[-0.5],posangles=["0.0deg"]):
 ### """ Make a component-list
 ###     params :
 ### """
  if(ncomps != len(fluxvals) or ncomps != len(minaxes) or  
     ncomps != len(majaxes) or  ncomps != len(ras) or  
     ncomps != len(decs) or  ncomps != len(spxs) or  
     ncomps != len(posangles) or  ncomps != len(shapes) ):
	  print "Please enter all params for each component"
	  return;
  os.system('rm -rf '+clname)
  refRA = float(qa.formxxx(ra0,'rad'))
  refDEC = float(qa.formxxx(dec0,'rad'))
  for comp in range(0,ncomps):
    thisdec = refDEC +  float(qa.formxxx( decs[comp] ,'rad'))
    thisra = refRA +  float(qa.formxxx( ras[comp] ,'rad'))   / cos(refDEC)
    raFormat  = qa.formxxx(str(thisra)+'rad','hms');
    decFormat = qa.formxxx(str(thisdec)+'rad','dms');
    d= 'J2000 '+raFormat+' '+decFormat;
    #print 'Comp ', comp , d
    if(shapes[comp]=="point"):
      cl.addcomponent(flux=fluxvals[comp],fluxunit="Jy", 
		      dir=str(d),
		      shape="point",freq=reffreq,
		      spectrumtype="spectral index",index=spxs[comp]); #[spxs[comp],0,0,0]);
    else:
      cl.addcomponent(flux=fluxvals[comp],fluxunit="Jy", 
		     dir=str(d),
		      shape="Gaussian",freq=reffreq,
                      majoraxis=qa.quantity(majaxes[comp], 'arcsec'),
                      minoraxis=qa.quantity(minaxes[comp], 'arcsec'),
                      positionangle=qa.quantity(posangles[comp], 'rad'),
                      spectrumtype="spectral index",index=spxs[comp]); #[spxs[comp],0,0,0]);
  cl.rename(filename=clname);
  cl.close();

###################################################

## Define the locations of your sources ( check that these locations ar
### field of you that you are simulating ! )

clname = 'mysources';
ra0="13:29:53.933";
dec0="-47.11.42.205";
reffreq = '100.0GHz'


# Make a list of components.
makeComponentList(clname=clname+'0.cl',ncomps=1,shapes=['point'],
        fluxvals=[[1.0,1.0,0.0,0.0]], ra0=ra0, dec0=dec0,
        minaxes=[""], majaxes=[""], posangles=[""],
        ras=["15.0arcsec"],decs=["15.0arcsec"],
        spxs=[0.0],reffreq=reffreq);

i = 0
max_srcs = 4
fluxvals = []
shapes = []
minaxes = []
majaxes = []
posangles = []
spxs = []

ras = ["0.31arcmin","0.29arcmin","0.33arcmin","0.32arcmin"]
decs = ["0.32arcmin","0.30arcmin","0.33arcmin","0.28arcmin"]

while i < max_srcs:
    fluxvals.append([200.0,200.0,0.0,0.0])
    shapes.append('gaussian')
    spxs.append(0.0)
    i = i + 1

minaxes=["1.4arcsec","1.2arcsec","1.2arcsec","1.0arcsec"]
majaxes=["2.0arcsec","1.5arcsec","2.8arcsec","2.5arcsec"]
posangles=["40.0deg","0.0deg","-60.0deg","-75.0deg"]

print ras
print decs

makeComponentList(clname=clname+'1.cl', ncomps=max_srcs,
        shapes=shapes, fluxvals=fluxvals, ra0=ra0, dec0=dec0,
        minaxes=minaxes, majaxes=majaxes, posangles=posangles,
        ras=ras, decs=decs, spxs=spxs, reffreq=reffreq);
                            
