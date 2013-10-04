
## Make a copy of one of your MSs, to use for this.
## It will fill in values in the MODEL_DATA column..... 
## which you don't need.

def writeApertures(msname='', imname='', pblimit=1e-04):

  im.open(msname,usescratch=True);
  # Select a tiny piece of data. 
  # This is just so it has something to run through.
  im.selectvis(spw='*',baseline='0&&1',scan='1');

  # Do the prediction
  im.defineimage();
  im.setoptions(ftmachine='mawproject',
                applypointingoffsets=False,
                dopbgriddingcorrections=False,
                cfcachedirname=msname+'.cfcache',
                pastep=360.0, rotpastep=360.0,pblimit=pblimit,
                psterm=False, aterm=True, wbawp=False,mterm=False,conjbeams=False);
  im.ft(model=imname,incremental=False);
  im.close();



#############################

#writeApertures(msname='MsForAp.ms',  imname='ModForAp.im')


