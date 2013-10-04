
def parang(msname,ifld,t):


    # Mean earth position
    tb.open(msname+'/ANTENNA')
    pos=tb.getcol('POSITION')
    meanpos=pl.mean(pos,1)
    frame=tb.getcolkeyword('POSITION','MEASINFO')['Ref']
    units=tb.getcolkeyword('POSITION','QuantumUnits')
    mpos=me.position(frame,
                     str(meanpos[0])+units[0],
                     str(meanpos[1])+units[1],
                     str(meanpos[2])+units[2])
    me.doframe(mpos)

    # _geodetic_ latitude
    latr=me.measure(mpos,'WGS84')['m1']['value']
    
    
    # Field directions
    tb.open(msname+'/FIELD')
    dirs=tb.getcol('DELAY_DIR',ifld,1)[:,0,0]  # get selected field
    tb.close()

    rah=dirs[0]*12.0/pi  # hours
    decr=dirs[1]         # rad
    
    tm=me.epoch('UTC',str(t)+'s')
    last=me.measure(tm,'LAST')['m0']['value']
    last-=floor(last)  # days
    last*=24.0  # hours
    ha=last-rah  # hours
    har=ha*2.0*pi/24.0
    
    parang=atan2( (cos(latr)*sin(har)),
                  (sin(latr)*cos(decr)-cos(latr)*sin(decr)*cos(har)) )


    return parang*180/pi



def printparang():
    tb.open('points.ms')
    tmlist = tb.getcol('TIME')
    tb.close()

    for tm in tmlist:
        val = parang(msname='points.ms', ifld=0, t=tm)
        if( val < 0 ):
            val = val + 360.0
        print tm, val
