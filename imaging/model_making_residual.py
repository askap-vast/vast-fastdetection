
import os
import sys


vis = sys.argv[-2]
imagename = sys.argv[-1]  # recommend in format of SBxxx_beamxxx

print("** NOTICE  ** path passed as arg:", vis, imagename)

imsize = 6144
cell = ['2.5arcsec']
deconvolver = 'multiscale'
nterm = 2
scales = [0, 10, 30]
niter = 0
weighting = 'briggs'  # natural, uniform, briggs
robust = 0.5
uvrange = '>200m'  # >1.0m
pbcor = False  # True, False
pblimit = -0.2  # 0.1, -0.2
threshold="0.75mJy"



# make a deep image
tclean(
    vis=vis,
    selectdata=True,
    field='',
    spw='',
    timerange='',
    uvrange=uvrange,
    antenna='',
    scan='',
    observation='',
    intent='',
    datacolumn='corrected',
    imagename=imagename,
    imsize=imsize,
    cell=cell,
    phasecenter='',
    stokes='I',
    projection='SIN',
    specmode='mfs',
    reffreq='',
    nchan=-1,
    start='',
    width='',
    outframe='LSRK',
    veltype='radio',
    restfreq=[],
    interpolation='linear',
    gridder='widefield',
    facets=2,
    wprojplanes=-1,
    vptable='',
    aterm=True,
    psterm=False,
    wbawp=True,
    conjbeams=False,
    cfcache='',
    computepastep=360.0,
    pblimit=pblimit,
    normtype='flatnoise',
    deconvolver=deconvolver,
    scales=scales,
    nterms=nterm,
    smallscalebias=0.0,
    restoration=True,
    restoringbeam=[],
    pbcor=pbcor,
    outlierfile='',
    weighting=weighting,
    robust=robust,
    npixels=0,
    uvtaper=[],
    niter=niter,
    gain=0.2,
    threshold=threshold,
    cycleniter=-1,
    cyclefactor=1.0,
    minpsffraction=0.02,
    maxpsffraction=0.8,
    interactive=False,
    usemask='user',
    mask='',
    pbmask=0.0,
    savemodel='modelcolumn',
    startmodel='',
    parallel=False
)
print('Made a deep image finished %s' % imagename)


print('nterm is', nterm)
exportfits(imagename="%s.image" % imagename, fitsimage="%s.image.fits" % imagename)


