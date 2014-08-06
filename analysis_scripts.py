import numpy as np

def autoAnalysis(testrun='',
                 pathname='/lustre/kkundert/Code/results/new_results/thesis_runs/'):

    f = open(pathname + testrun + '/analysis_results.txt', 'w')

    for u in xrange(3):
        if u == 0:
            dirname = pathname + testrun + '/Data000-points'
        elif u == 1:
            dirname = pathname + testrun + '/Data000-extended'
        else:
            dirname = pathname + testrun + '/Data000-M51'

        imname = dirname + '/points.true.im'
        trueimage = imname + '.predict'
        cleanimage = imname + '.clean.image'
        #cleanimage = dirname + '/simple.image'
        smoothimage = imname + '.smooth'
        diffimage = imname + '.diffmap'
        pbdir = dirname + '/primary-beams'
        numants = 34
        avgpb = pbdir + '/avgpb'
        pbcor = dirname + '/points.true.im.pbcor'
        diff = 0.0
        chinorm = 0.0
        rms = 0.0

        averagePB(dirname=pbdir,numants=numants,avgpb=avgpb)
        pbCorImage(cleanimage=cleanimage,beam=avgpb,outfile=pbcor)
        smoothTrueImage(trueimage=trueimage,beam=avgpb,tempimage=dirname+'/temp',
                        smoothimage=smoothimage,major='1.0arcsec',minor='1.0arcsec',
                        pa='0deg')
        makeDifferenceImage(trueimage=smoothimage,cleanimage=cleanimage,outfile=diffimage)
        if u == 0:
            diff = peakDiff(imtrue=smoothimage,imobs=cleanimage)
            print 'peak difference is ' + str(diff)
            f.write(str(diff)+'\n')
        else:
            rms = rootMeanSquare(diffimage=diffimage)
            print 'rms is: ' + str(rms)
            f.write(str(rms))
            #sigma = sigma(imtrue=smoothimage,imobs=cleanimage,beam=avgpb)
            #print 'sigma is: ' + str(sigma)
            #f.write(str(sigma))
        
    f.close()
            

def imageFidelity(imtrue='',imobs=''):

    # calculate image fidelity comparing model image to observed image
    
    ia.open(imtrue)
    true = ia.getchunk()
    ia.close()

    ia.open(imobs)
    obs = ia.getchunk()
    ia.close()

    fidelity = 0.0
    if len(true) == len(obs):
        N = len(true)
        for u in xrange(N):
            for v in xrange(N):
                val = abs(true[u][v]) / abs(obs[u][v] - true[u][v])
                fidelity += val
                val = 0

    else:
        print "array sizes don't match, please choose different images"

    return fidelity[0][0]

def chiSquareNormalized(imtrue='',imobs='',beam=''):

    # calculate chi-square distribution comparing model image to observed image
    # normalized by expected noise levels
    
    ia.open(imtrue)
    true = ia.getchunk()
    ia.close()

    ia.open(imobs)
    obs = ia.getchunk()
    ia.close()

    ia.open(beam)
    pb = ia.getchunk()
    ia.close()

    difference = 0.0
    pbsum = 0.0
    noise = 1.0 * 10**-3.5
    if len(true) == len(obs):
        N = len(true)
        for u in xrange(N):
            for v in xrange(N):
                val = pb[u,v,0,0] * (obs[u,v,0,0] - true[u,v,0,0])
                difference += val**2
                pbsum += pb[u,v,0,0]
                val = 0
        variance = difference / ((N**2)-1)
        chinorm = variance / (pbsum * noise**2)

    else:
        print "array sizes don't match, please choose different images"

    return chinorm

def sigma(imtrue='',imobs='',beam=''):

    # calculate chi-square distribution comparing model image to observed image
    # normalized by expected noise levels
    
    ia.open(imtrue)
    true = ia.getchunk()
    ia.close()

    ia.open(imobs)
    obs = ia.getchunk()
    ia.close()

    ia.open(beam)
    pb = ia.getchunk()
    ia.close()

    difference = 0.0
    pbsum = 0.0
    noise = 1.0 * 10**-2
    if len(true) == len(obs):
        N = len(true)
        for u in xrange(N):
            for v in xrange(N):
                val = pb[u,v,0,0] * (obs[u,v,0,0] - true[u,v,0,0])
                difference += val**2
                pbsum += pb[u,v,0,0]
                val = 0
        variance = difference / ((N**2)-1)
        sigma = variance / (pbsum)

    else:
        print "array sizes don't match, please choose different images"

    return sigma

def rootMeanSquare(diffimage=''):

    ia.open(diffimage)
    arr = ia.getchunk()
    ia.close()

    N = len(arr)
    testarr = arr[N/2+N/4:N,N/2+N/4:N,0,0]

    rms = sqrt(np.mean(testarr**2))
    return rms
    
def pbCorImage(cleanimage='',beam='',outfile=''):

    os.system('cp -r ' + cleanimage + ' ' + outfile)
    
    ia.open(cleanimage)
    image = ia.getchunk()
    ia.close()

    ia.open(beam)
    pb = ia.getchunk()
    ia.close()

    N = len(pb)

    #implement cut-off point so that corrected image doesn't explode at edge
    for u in xrange(N):
        for v in xrange(N):
            if pb[u][v][0][0] < 0.05:
                pb[u][v][0][0] = 1

    pbcorimage = image[:,:,0,0] / pb[:,:,0,0]

    ia.open(outfile)
    vals = ia.getchunk()
    vals[:,:,0,0] = pbcorimage
    ia.putchunk(vals)
    ia.close()

def makeDifferenceImage(trueimage='',cleanimage='',outfile=''):

    os.system('cp -r ' + cleanimage + ' ' + outfile)
    
    ia.open(trueimage)
    true = ia.getchunk()
    ia.close()

    ia.open(cleanimage)
    clean = ia.getchunk()
    ia.close()

    true = true[:,:,0,0]
    clean = clean[:,:,0,0]

    diff = np.zeros(clean.shape)
    if len(true) == len(clean):
       diff = clean - true
       ia.open(outfile)
       vals = ia.getchunk()
       vals[:,:,0,0] = diff
       ia.putchunk(vals)
       ia.close()
    else:
        print "array sizes don't match, please choose different images"
        return

def averagePB(dirname='',numants=34,avgpb=''):

    os.system('cp -r ' + dirname+'/pb00\&\&01real ' + avgpb)

    ia.open(avgpb)
    vals = ia.getchunk()
    vals = np.zeros(vals.shape)
    ia.putchunk(vals)
    ia.close()

    avgdenom = 0
    for u in xrange(numants):
        for v in xrange(numants):
            if u < v:
                imname = dirname+'/pb%02d&&%02dreal' % (u, v)
                ia.open(imname)
                arr = ia.getchunk()
                ia.close()
                vals += arr
                avgdenom += 1
                
    vals = vals / float(avgdenom)

    ia.open(avgpb)
    ia.putchunk(vals)
    ia.close()

def smoothTrueImage(trueimage='',beam='',tempimage='',smoothimage='',
                    major='1.0arcsec',minor='1.0arcsec',pa='0deg'):

    os.system('cp -r ' + trueimage + ' ' + tempimage)

    ia.open(trueimage)
    true = ia.getchunk()
    ia.close()

    ia.open(beam)
    pb = ia.getchunk()
    ia.close()

    ia.open(tempimage)
    smooth = ia.getchunk()
    smooth = true * pb
    ia.putchunk(smooth)
    smoothed = ia.convolve2d(outfile=smoothimage,major=major,minor=minor,pa=pa,type='gaussian',scale=1,overwrite=True)
    ia.close()
    smoothed.close()

    os.system('rm -rf ' + tempimage)

def peakDiff(imtrue='',imobs=''):

    ia.open(imtrue)
    true = ia.getchunk()
    ia.close()

    ia.open(imobs)
    obs = ia.getchunk()
    ia.close()

    max_true = np.amax(true)
    max_obs = np.amax(obs)
    peak_diff = max_true - max_obs

    return peak_diff

