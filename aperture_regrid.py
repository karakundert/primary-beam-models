import numpy as np

def aperRegrid():

    aptype='DV'
    pathname = '/lustre/kkundert/Code/alma_apertures/thesis_apertures/'
    templateimage =\
    '/lustre/kkundert/Code/alma_apertures/DV/uid___A002_X67509f_X660_APC-DV05-H-LSB-beam_complex.aperture.real.regridded'

    danum = [41,42,43,44,45,46,47,48,52,53,54,55,56,57,58]
    dvnum = [1,2,3,5,6,7,8,9,10,11,13,14,15,17,18,19,20,23,24,25]

    num_ant = len(dvnum)

    for i in xrange(num_ant):
        imname = pathname+aptype+"%02d-H.imag" % dvnum[i]
        imregrid(imagename=imname,template=templateimage,output=imname+'.regridded',axes=[0,1],overwrite=True)
        imname = pathname+aptype+"%02d-H.real" % dvnum[i]
        imregrid(imagename=imname,template=templateimage,output=imname+'.regridded',axes=[0,1],overwrite=True)
        imname = pathname+aptype+"%02d-V.imag" % dvnum[i]
        imregrid(imagename=imname,template=templateimage,output=imname+'.regridded',axes=[0,1],overwrite=True)
        imname = pathname+aptype+"%02d-V.real" % dvnum[i]
        imregrid(imagename=imname,template=templateimage,output=imname+'.regridded',axes=[0,1],overwrite=True)
