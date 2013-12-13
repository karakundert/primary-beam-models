def makeComponentList(clname="sim.cl",reffreq="1.5GHz",
                        ra0="19:59:28.500",dec0="+40.44.01.50",
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
ra0="19:59:28.500";
dec0="+40.44.01.50";
reffreq = '6.0GHz'


# Make a list of components.
# components placed such that there is one source per simulation, place
# full-power, half-power, and around 10-20% power points of the beam.
#makeComponentList(clname=clname+'0.cl',ncomps=1,shapes=['point'],
#        fluxvals=[[1.0,0.0,0.0,1.0]],
#        minaxes=[""], majaxes=[""], posangles=[""],
#        ras=["0.0arcmin"],decs=["0.0arcmin"],
#        spxs=[0.0],reffreq=reffreq);
makeComponentList(clname=clname+'0.cl',ncomps=1,shapes=['point'],
        fluxvals=[[1.0,0.0,0.0,1.0]],
        minaxes=[""], majaxes=[""], posangles=[""],
        ras=["2.0arcmin"],decs=["2.0arcmin"],
        spxs=[0.0],reffreq=reffreq);
                            
