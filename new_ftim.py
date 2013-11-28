import numpy as np;
def ftim(imagename,outimagename,tbt=tb,npt=np):
    tbt.open(imagename,nomodify=false);
    data=tbt.getcol('map');
    tbt.close();
    idata = np.fft.ifftshift(data)
    fdata=(np.fft.fftn(idata));
    sfdata=np.fft.fftshift(fdata);
    tbt.open(outimagename,nomodify=false);
    tbt.putcol('map',sfdata);
    tbt.close();


def ftim2(imagename='MeasuredVoltage_DA_H_LSB.im',outimagename='MeasuredAperture_DA_H_LSB.im', csystemplate='kimtruecollapse.im'):

    ## Clear the way...
    if os.path.exists(outimagename):
        os.system('rm -rf '+ outimagename)
    if os.path.exists('tmp.real'):
        os.system('rm -rf tmp.real*')

    # Make a real image with the same coordinate system
    ## Read shape and cell size
    tb.open(imagename)
    csysrec = tb.getkeyword('coords')
    shp = (tb.getcol('map')).shape
    tb.close()
    shp = shp[0:len(shp)-1]
    cdelt = csysrec['linear0']['cdelt']
    cunit = csysrec['linear0']['units']
    crpix = csysrec['linear0']['crpix']
    print shp, cdelt, cunit

    # Read a valid csys and modify it with the above values
    ia.open(csystemplate)
    tsys = ia.coordsys()
    tmpshp=ia.shape()
    ia.close()

    cdelt[0] = (qa.convert( qa.quantity( cdelt[0], cunit[0] ) , 'rad' ))['value']
    cdelt[1] = (qa.convert( qa.quantity( cdelt[1], cunit[1] ) , 'rad' ))['value']

    tinc = tsys.increment()['numeric']
    tinc[0] = cdelt[0]
    tinc[1] = cdelt[1]
    tinc[3] = 1e+10

    tsys.setincrement(tinc)

    refval = tsys.referencevalue()['numeric']
    refval[3] = 1e+11
    tsys.setreferencevalue( refval )

    trpix = tsys.referencepixel()['numeric']
    trpix[0] = crpix[0]
    trpix[1] = crpix[1]
    tsys.setreferencepixel( trpix)

    for axx in range(0,len(shp)):
        tmpshp[axx] = shp[axx]

    # Set Stokes naxis = 1
    tmpshp[2]=1

    tsys.setstokes(['I'])

    # Make an empty real image.
    print tmpshp

    newia = casac.image()
    newia.fromshape(outfile='tmp.real', shape=tmpshp, csys=tsys.torecord(),overwrite=True)
    newia.close()

#    ia.open('tmp.real')
#    ia.summary()
#    ia.close()
#    print tb.showcache()

    # Check the header
#    imhead('tmp.real')

    ## Take its FT and pull out its coordsys.
    ia.open('tmp.real')
    ia.fft(real='tmp.real.ft') #,axes=[0,1])
    ia.close()
    ia.open('tmp.real.ft')
    uvcsys = ia.coordsys()
    uvshp = ia.shape()
    ia.close()

    print 'uv cell size : ', uvcsys.increment()

    ## Take the FT of the actual complex voltage pattern
    tb.open(imagename,nomodify=false);
    data=tb.getcol('map');
    tb.close();
    idata = np.fft.ifftshift(data)
    fdata=(np.fft.fftn(idata));
    sfdata=np.fft.fftshift(fdata);

    print sfdata.shape, uvshp

    ssf = sfdata.reshape( uvshp )

    print ssf.shape

    # Save real and imag parts separately using the UV coordinate system
    newia = casac.image()
    newia.fromarray(outfile=outimagename+'.real',pixels=np.real(ssf),csys=uvcsys.torecord(),linear=F,overwrite=T)
    newia.close()

    newia = casac.image()
    newia.fromarray(outfile=outimagename+'.imag',pixels=np.imag(ssf),csys=uvcsys.torecord(),linear=F,overwrite=T)
    newia.close()


#####################################################

def prepareuvgrid( targetgrid='reaperture_pol9.im', outfile='templateuvgrid.im' ):
    os.system('rm -rf '+outfile)
    imcollapse(imagename=targetgrid, function='mean', outfile=outfile,axes=[2,3])
    ia.open(outfile)
    csys = ia.coordsys()
    csys.setstokes('I')
    incr = csys.increment()['numeric']
    incr[3]=1e+10
    csys.setincrement(incr)
    ia.setcoordsys( csys.torecord() )
    ia.close()


def uvregrid(sourceimage='MeasuredAperture_DA_H_LSB.im.real', outimage='MeasuredAperture_DA_H_LSB.im.real.regridded', templateimage='templateuvgrid.im'):

    
    ia.open(templateimage)
    csys = ia.coordsys()
    shp = ia.shape()
    ia.close()

    ia.open(sourceimage)
    qq = ia.regrid(outfile=outimage, shape=shp,csys=csys.torecord(),axes=[0,1],overwrite=True)
    qq.close()
    ia.close()


################################################

####  Convert Measured ALMA beams into Aperture functions with correct UV coords.

def convert_almabeams_to_aperture(mypath='/lustre/rurvashi/ALMAsim/DiahPBs/DA'):
    from os import listdir
    from os.path import isfile, join
    voltlist = [ f for f in listdir(mypath) if ( not isfile(join(mypath,f)) and f.count('_complex') ) ]
    
    print len(voltlist)

    cnt=0
    for voltage in voltlist:
        print cnt , ' : Converting ', voltage
        ftim2(imagename=mypath+'/'+voltage ,outimagename=mypath+'/'+voltage+'.aperture', csystemplate='kimtruecollapse.im')


def convert_almabeam_to_aperture_and_regrid(voltage='',aperture=''):
    os.system('rm -rf tmpap*')
    ftim2(imagename=voltage ,outimagename='tmpap', csystemplate='kimtruecollapse.im')
    uvregrid(sourceimage='tmpap.real', outimage=aperture+'.real', templateimage='templateuvgrid.im')
    uvregrid(sourceimage='tmpap.imag', outimage=aperture+'.imag', templateimage='templateuvgrid.im')


def drawfigs(mypath='/lustre/rurvashi/ALMAsim/DiahPBs/DA',outname='Figures/DA'):
    from os import listdir
    from os.path import isfile, join
    powerlist = [ f for f in listdir(mypath) if ( not isfile(join(mypath,f)) and f.count('real') ) ]
    
    print len(powerlist)

    print powerlist

    if(1):
        ia.open(mypath+'/'+powerlist[0])
        powshp = ia.shape()
        ia.close()
        npix = powshp[0]
        lmin=int(npix/2 - npix/4)
        lmax=int(npix/2 + npix/4)
        
        print powshp, lmin, lmax
    
        ### Make figs for powers
        outpow = ''

        pl.figure(figsize=(5,5))
        for pfig in powerlist:
            print 'Opening ' + pfig
            ia.open(mypath+'/'+pfig)
            pix = ia.getchunk()
            ia.close()
            pl.clf()
            print pix.shape
            pl.imshow( pix[lmin:lmax,lmin:lmax,0,0] , aspect='equal',  vmin=0.0,vmax=1.0, origin='upper' ,cmap='jet')
            pl.axis('off')
            pl.savefig(outname+'.power.'+pfig+'.png')
            outpow = outpow + ' ' + outname+'.power.'+pfig+'.png' 

        os.system('convert -delay 20 -loop 0 ' + outpow + ' ' + outname + '.power.gif')
            


######################



