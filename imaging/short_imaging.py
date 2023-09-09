#!/usr/bin/env sh 

import os
import sys
import numpy as np

vis = sys.argv[-3]
imagename = sys.argv[-2]  # recommend in format of SBxxx_beamxx
step = sys.argv[-1]  # images time length - in unit of seconds
print "** NOTICE  ** the args passed to casa was:" + vis + imagename + step

step = float(step)
# corrected, data
datacolumn = "corrected"
imsize = 2800  # 3000? 4096? 2048?
cell = ['2.5arcsec']
uvrange = '>200m'
# mtmfs, multiscale
deconvolver = 'multiscale'
scales = [0]
nterms = 1
# none, modelcolumn
savemodel = "none"
niter = 0
# weigithing
weighting = 'briggs'
robust = 0.5


# open table and read time column
tb.open(vis)

############################
# Something wrong with if statement in CASA environment???? - need to go back to this later 
############################

# the ASKAP resolution is 9.97s, so we need a different workaround for 10s images
# if ( step <= 10 ):
#     from collections import Counter
#     times = Counter(tb.getcol('TIME')).keys()
#     times.sort()
#     # time in fits is the middle time (but it doesn't matter as 10s is the sampling time, so it's the start, middle, and end time as well)
#     times = np.array(times) - step/2
#     print "Input time length is {}s, using a slightly different way to get short images".format(step)

# else:
#     times = tb.getcol('TIME')
#     interval = tb.getcol('INTERVAL')[0]
#     # get the times array, change unit to second
#     # time in fits is therefore the START time 
#     times = np.arange(start=np.min(times)-interval/2,
#                       stop=np.max(times)+interval/2, step=step)
#     print "Input time length is {}s, i.e., {:.1f}m".format(step, step/60)

from collections import Counter
times = Counter(tb.getcol('TIME')).keys()
times.sort()
# time in fits is the middle time (but it doesn't matter as 10s is the sampling time, so it's the start, middle, and end time as well)
times = np.array(times) - step/2
print "Input time length is {}s, using a slightly different way to get short images".format(step)


tb.close()


# convert MJD to standard UTC
for j in range(times.shape[0]):
    stime = qa.time(qa.quantity(times[j], 's'), form="ymd")[0]
    etime = qa.time(qa.quantity(times[j]+step, 's'), form='ymd')[0]
    imagename_j = '%s_%s' % (imagename, j)
    timerange = '%s~%s' % (stime, etime)

    print '%s, %s' % (os.getcwd(), timerange)

    tclean(vis=vis, selectdata=True, field='', spw='', timerange=timerange, uvrange=uvrange, antenna='', scan='', observation='', intent='', datacolumn=datacolumn, imagename=imagename_j, imsize=imsize, cell=cell, phasecenter='', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', gridder='widefield', facets=1, chanchunks=1, wprojplanes=-1, vptable='', aterm=True,
           psterm=False, wbawp=True, conjbeams=False, cfcache='', computepastep=360.0, pblimit=-0.2, normtype='flatnoise', deconvolver=deconvolver, scales=scales, nterms=nterms, smallscalebias=0.0, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting=weighting, robust=robust, npixels=0, uvtaper=[], niter=niter, gain=0.1, threshold=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.02, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, savemodel=savemodel)

    # change to fits format
    exportfits(imagename='%s.image' % imagename_j,
               fitsimage='%s.fits' % imagename_j)

print '%s Finished. ' % imagename
