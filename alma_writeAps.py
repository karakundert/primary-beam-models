
## Make a copy of one of your MSs, to use for this.
## It will fill in values in the MODEL_DATA column..... which you don't need.

##writeApertures(msname='almapoints.ms', imname='almatrue.im'

def writeApertures(msname='', imname='',anttype='DA'):

  # Alma antenna types are sent in as a fudge....
  # Possible types : DA, DV, PM, CM.

  fp = open('anttype.txt','w')
  fp.writelines( anttype )
  fp.close()

  # This needs to be an ALMA dataset
  im.open(msname,usescratch=True);
  # Select a tiny piece of data. This is to set the timestep/parangle to use.
  im.selectvis(spw='*',baseline='0&&1',scan='1');

  # Do the prediction
  im.defineimage();
  im.setoptions(ftmachine='mawproject',
                applypointingoffsets=False, 
                dopbgriddingcorrections=False, 
                cfcachedirname=msname+'.cfcache', 
                pastep=360.0, rotpastep=360.0,pblimit=1e-04, 
                psterm=False, aterm=True, wbawp=False,mterm=False,conjbeams=False);
  im.ft(model=imname,incremental=False);
  im.close();

  dirname = 'APdir.'+anttype
  if not os.path.exists(dirname):
    os.system('mkdir ' + dirname)

  os.system('mv reaperture_pol9.im '+dirname)
  os.system('mv imaperture_pol9.im '+dirname)
  os.system('mv reaperture_pol12.im '+dirname)
  os.system('mv imaperture_pol12.im '+dirname)

  os.system('mv ap_voltage_pol9.im '+dirname)
  os.system('mv ap_voltage_pol12.im '+dirname)


#############################

#writeApertures(msname='MsForAp.ms',  imname='ModForAp.im')


