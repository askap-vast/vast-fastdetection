#!/usr/bin/env sh

import sys

ms_file = sys.argv[-1]
print("** NOTICE  ** the arg3 passed to casa was :"+ms_file)

beam_idx_b = ms_file.find('beam') 
beam_idx_e = ms_file.find('_', beam_idx_b + 4)
beam_t = ms_file[beam_idx_b:beam_idx_e]
beam_list = [beam_t]
#epoch_list = ['15191']
epoch_tag = 'scienceData_'
epoch_idx_b = ms_file.find(epoch_tag)
epoch_idx_e = ms_file.find('_', epoch_idx_b + len(epoch_tag))
epoch_t = ms_file[epoch_idx_b+len(epoch_tag):epoch_idx_e]
epoch_list = [epoch_t]

# don't need to change part
# corrected, data
datacolumn = "corrected"
imsize = 3000 # 4000
cell = ['2.5arcsec']
uvrange = '>1m'
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

# python pakeage
import os
import numpy as np
from mpi4py import MPI

# parallelly
comm = MPI.COMM_WORLD
num_process = comm.Get_size()
rank = comm.Get_rank()

beam_split = np.array_split(beam_list, num_process)

# started to 15 minute imaging
for i in range(len(beam_split[rank])):
    beam = beam_split[rank][i]
    print 'total %s, now %s, in rank %s (e.g. %s), length of this rank %s, num_process %s' % (len(beam_list), i, rank, beam, len(beam_split[rank]), num_process)

    path_idx_e = ms_file.find("corrected_data")
    path = ms_file[:path_idx_e]
    #workdir1 = path + beam + "/"
    workdir1 = path + "images/"
    if not os.path.exists(workdir1):
        os.mkdir(workdir1)

    for epoch in epoch_list:
        #filename = '/o9000/ASKAP/VAST/ligo_variables/corrected_data/scienceData_SB%s_S190814bv.%s_averaged_cal.ms/' % (epoch, beam)
        filename = ms_file
        #workdir2 = '/o9000/ASKAP/VAST/%s/SB%s' % (beam, epoch)
        workdir2 = workdir1 + epoch + "/"
#        if not os.path.exists(workdir2):
#            os.mkdir(workdir2)
        os.chdir(workdir2)

        # 15 minute imaging!!!!
        # open table and read time column
        tb.open(filename)
        times = tb.getcol('TIME')
        interval = tb.getcol('INTERVAL')[0]
        # get the times array, change unit to second
        times = np.arange(start=np.min(times)-interval/2, stop=np.max(times)+interval/2, step=15*60)
        tb.close()

        # convert MJD to standard UTC
        for j in range(times.shape[0]):
            stime = qa.time(qa.quantity(times[j], 's'), form="ymd")[0]
            etime = qa.time(qa.quantity(times[j]+15*60, 's'), form='ymd')[0]
            imagename = '%s_%s' % (beam, j)
            timerange = '%s~%s' % (stime, etime)

            print '%s, %s' % (os.getcwd(), timerange)

            tclean(vis=filename, selectdata=True, field='', spw='', timerange=timerange, uvrange=uvrange, antenna='', scan='', observation='', intent='', datacolumn=datacolumn, imagename=imagename, imsize=imsize, cell=cell, phasecenter='', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', gridder='widefield', facets=1, chanchunks=1, wprojplanes=-1, vptable='', aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', computepastep=360.0, pblimit=-0.2, normtype='flatnoise', deconvolver=deconvolver, scales=scales, nterms=nterms, smallscalebias=0.0, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting=weighting, robust=robust, npixels=0, uvtaper=[], niter=niter, gain=0.1, threshold=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.02, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, savemodel=savemodel)
            # change to fits format
            exportfits(imagename='%s.image'%imagename, fitsimage='%s.fits'%imagename)

    print '%s Finished. ' % beam

