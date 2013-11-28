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
