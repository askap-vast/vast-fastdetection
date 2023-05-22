#!/usr/bin/env sh

import os
import sys


VIS = sys.argv[-4]
DIR = sys.argv[-3]
SBID = sys.argv[-2] # in format of SBxxx
BEAM = sys.argv[-1] # in format of beam00
print("** NOTICE  ** path passed as arg:", VIS, DIR, SBID, BEAM)


imsize = 10000
cell = ['2.5arcsec']
deconvolver = 'mtmfs'
nterm = 2
scales = [0, 5, 10, 15, 25]
niter = 10000
weighting = 'briggs' # natural, uniform, briggs
robust = 0.5 
uvrange = '>1.0m' # >1.0m
pbcor = False # True, False
pblimit = -0.2 # 0.1, -0.2


os.chdir(DIR)

imagename = '%s_%s' % (SBID, BEAM)
print('Reading %s, saving in %s as %s' % (VIS, DIR, IMAGE))


    # make a deep image
tclean(vis=VIS, selectdata=True, field='', spw='', timerange='', uvrange=uvrange, antenna='', scan='', observation='', intent='', datacolumn='data', imagename=imagename, imsize=imsize, cell=cell, phasecenter='', stokes='I', projection='SIN', specmode='mfs', reffreq='', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', gridder='widefield', facets=1, chanchunks=1, wprojplanes=-1, vptable='', aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', computepastep=360.0, pblimit=pblimit, normtype='flatnoise', deconvolver=deconvolver, scales=scales, nterms=nterm, smallscalebias=0.0, restoration=True, restoringbeam=[], pbcor=pbcor, outlierfile='', weighting=weighting, robust=robust, npixels=0, uvtaper=[], niter=niter, gain=0.1, threshold=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.02, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, savemodel='modelcolumn', startmodel='', parallel=False)
print('Made a deep image finished %s' % beam_nbr)

# subtract model
uvsub(vis=filename)
#print 'Model-subtraction %s Finished.' % beam_split[rank][i]
print('Model-subtraction %s Finished.' % beam_nbr)

# convert image to fits file
exportfits(imagename="%s.image.tt0"%imagename, fitsimage="%s.image.tt0.fits"%imagename)

