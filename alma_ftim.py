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
    ia.fromshape(outfile='tmp.real', shape=tuple(tmpshp), csys=tsys.torecord(),overwrite=True)
    ia.close()

    # Check the header
    imhead('tmp.real')

    ## Take its FT and pull out its coordsys.
    ia.open('tmp.real')
    ia.fft(real='tmp.real.ft',axes=[0,1])
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

