#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 20:41:50 2023
@author: ywan3191
Get ready for various scripts. 
For each beam, we want: 
1. Query CASDA to get visibility download url + selavy catalogue download url. 
Generate **bash_GETDATA_beam??.sh** and run the script to download visibilities 
and catalogues. (Double check the completion of visibilities downloading and 
untar visibilities) 
2. Run slurm_FIXDATA_beam??.sh (including rescale and fix dir)
3. Run slurm_MODELING_beam??.sh (deep modeling)
4. Run slurm_IMGFAST_beam??.sh (fast imaging, in either time interval)
5. Run slurm_SELCAND_beam??.sh (select candidates)
"""


from astroquery.utils.tap.core import TapPlus
from vaster.structure import DataBasic
from vaster.vtools import measure_running_time

import time
import os
import sys
import yaml
import getpass
import requests
import xmltodict
import argparse 

import logging
logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='prepare_scripts', 
        description='prepare various scripts to run VASTER', 
        epilog='Example usage: prepare_scripts 60114', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('sbids', nargs='+', help='input sbid for prepare scripts, number only')
    parser.add_argument('--outdir', type=str, default='.', help='output directory')
    parser.add_argument('--config', type=str, default='config.yml', help='configuration file')
    parser.add_argument('--opal-from-env', action='store_true', 
                        help='read opal information from environment variables $OPAL_USER and $OPAL_PWD')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()
    make_verbose(args)
    logger.info(args)
    logger.info('Total of %s SBIDs to prepare', len(args.sbids))

    config = read_config(args.config)
    args.config = os.path.abspath(args.config)
    logger.info('Reading configuration file %s', args.config)
    logger.info('Running in %s mode', config['MACHINE'])
    args.username, args.password = get_opal(args)

    

    for i, sbid in enumerate(args.sbids):
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))

        databasic = DataBasic(sbid, args.outdir)
        args.paths = databasic.paths
        args.steps = databasic.steps
        logger.info('steps: %s', args.steps)
        build_folder_structure(args.paths)
        copy_config(args, )

        vis, cat, img = query_casda(sbid)
        if len(vis) != 36:
            logger.warning('Number of visibilities is not 36. Skip SB%s', sbid)
            logger.warning('You can download the visibility manually. ')
            continue
        
        format_casa_modeling(args, config)
        format_casa_imgfast(args, config)
        prepare_downloads(args, sbid, cat, img, vis)

        if config['MACHINE'] == 'bash':
            format_bash(args, config, sbid, vis, cat)
        elif config['MACHINE'] == 'ozstar':
            format_ozstar(args, config, sbid, vis, cat)

        logger.info('SB%s preparation finish.', sbid)

    end_time = time.time()
    measure_running_time(start_time, end_time, 3)

  

def read_config(fname):
    # read configuration file 
    if os.path.isfile(fname):
        with open(fname, 'r') as yaml_file:
            config = yaml.safe_load(yaml_file)
        return config
    
    else:
        logger.error('Cannot find configuration file')
        logger.error('%s does not exists', fname)
        sys.exit()

def copy_config(args, ):
    config_in = args.config
    config_out = os.path.join(args.paths['path'], 'config.yml')
    logger.info('copy configuration file')
    txt = f'cp {config_in} {config_out}'
    logger.info('Executing "%s"', txt)
    os.system(txt)


def make_verbose(args):
    if args.verbose:
        logger.warning("verbose mode")
        logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')



def build_folder_structure(paths):
    ############################
    # Build file saving system structure 
    ############################
    for key, dir in paths.items():
        if not os.path.exists(dir):
            os.makedirs(dir, exist_ok=True)
            logger.info('%s: Create new directory %s', key, dir)
        else:
            logger.warning('%s: Directory %s exists.', key, dir)


def query_casda(sbid):
    ############################
    # Find visibilities, selavy catalogues and fits images from CASDA
    ############################
    tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
    job = tap.launch_job_async("SELECT * FROM ivoa.obscore WHERE obs_id='{}' AND (dataproduct_type='{}' OR dataproduct_subtype='{}' OR dataproduct_subtype='{}') ".format(
        'ASKAP-'+str(sbid), 'visibility', 'catalogue.continuum.component', 'cont.restored.t0'))

    r = job.get_results()

    vis = r[r['dataproduct_type'] == 'visibility']
    cat = r[r['dataproduct_subtype'] == 'catalogue.continuum.component']
    img = r[r['dataproduct_subtype'] == 'cont.restored.t0']

    vis.sort('filename')

    logger.info('Found {} visibilities'.format(len(vis)))
    logger.info(vis[['filename', 'access_url']])
    logger.info('Found {} selavy catalogues'.format(len(cat)))
    logger.info(cat[['filename', 'access_url']])
    logger.info('Found {} images'.format(len(img)))
    logger.warning(img[['filename', 'access_url']])

    return vis, cat, img
    

def get_opal(args):
    ############################
    # Get download urls
    ############################
    if args.opal_from_env:
        logger.info('Read OPAL from environment variables')
        username = os.getenv('OPAL_USER')
        logger.info('OPAL username: %s', username)
        password = os.getenv('OPAL_PWD')
    else: 
        username = input("Enter OPAL username: ")
        password = getpass.getpass(str("Enter OPAL password: "))
    return username, password


def get_modeling_params(config):
    if config['SURVEY'] == 'emu':
        params = config['MODELCLEAN']['EMU']
        logger.info('EMU clean mode ')
    elif config['SURVEY'] == 'vast':
        params = config['MODELCLEAN']['VAST']
        logger.info('VAST clean mode')
    else:
        params = config['MODELCLEAN']['OTHER']
        logger.info('User self-defined clean mode')
    return params


def format_bash(args, config, sbid, vis, cat):
    ############################
    # Generate scripts for one beam data
    ############################
    for idx in range(36):
        savename = os.path.join(args.paths['path_scripts'], 'bash_PROCESSING_beam{:02d}.sh'.format(idx))
        oname = f'SB{sbid}_beam{idx:02d}'
        with open(savename, 'w') as fw:
            write_basetxt_bash(fw, sbid, savename)
            url = get_url(vis[idx]['access_url'], args.username, args.password)
            filename = vis[idx]['filename']
            write_download_vis_txt(args, fw, idx, filename, url)
            write_untar_vis_txt(args, fw, idx, filename)
            write_fixdata_txt(args, fw, idx, filename)
            write_run_casa_txt(args, fw, idx, filename, oname, config, mode='modeling')
            write_run_casa_txt(args, fw, idx, filename, oname, config, mode='imaging')

            path_log = os.path.join(args.paths['path_logs'], 'bash_SELCAND_{}.log'.format(oname)) # output log files 
            affix = f' > {path_log} 2>&1'
            write_selcand_txt(args, fw, idx, oname, cat, affix=affix)
            write_clndata_txt(args, fw, idx)
            logger.info('Writing {}'.format(savename))


def format_ozstar(args, config, sbid, vis, cat):
    ############################
    # Generate scripts for one beam data
    ############################
    for idx in range(36):
        oname = f'SB{sbid}_beam{idx:02d}'

        for step in args.steps:
            logger.debug('%s: Run VASTER step %s...', oname, step)

            if step == 'GETDATA' or step == 'UNTAR':
                savename = os.path.join(args.paths['path_scripts'], f'bash_{step}_beam{idx:02d}.sh')
                with open(savename, 'w') as fw:
                    write_basetxt_bash(fw, sbid, savename)
                    url = get_url(vis[idx]['access_url'], args.username, args.password)
                    filename = vis[idx]['filename']
                    if step == 'GETDATA':
                        write_download_vis_txt(args, fw, idx, filename, url)
                    elif step == 'UNTAR':
                        write_untar_vis_txt(args, fw, idx, filename)

            else:
                params, savename = prepare_steps_ozstar(args, idx, config, sbid, oname, step)
                with open(savename, 'w') as fw:
                    write_basetxt_ozstar(fw, sbid, savename, params)
                    if step == 'FIXDATA':
                        write_moduleload_ozstar(fw, config)
                        write_fixdata_txt(args, fw, idx, filename, prefix='srun time ')
                        write_module_unload_ozstar(fw)
                    elif step == 'MODELING':
                        write_run_casa_txt(args, fw, idx, filename, oname, config, mode='modeling', prefix='srun time ')
                    elif step == 'IMGFAST':
                        write_run_casa_txt(args, fw, idx, filename, oname, config, mode='imaging', prefix='srun time ')
                    elif step == 'SELCAND':
                        write_moduleload_ozstar(fw, config)
                        write_selcand_txt(args, fw, idx, oname, cat, prefix='srun time ')
                        write_module_unload_ozstar(fw)
                    elif step == 'CLNDATA':
                        write_clndata_txt(args, fw, idx)
                        
                    write_endtxt_ozstar(fw, sbid, savename, params)
                
            logger.info('Writing {}'.format(savename))



def prepare_downloads(args, sbid, cat, img, vis):
    downloads = ['selavy', 'mosaic_images', 'visibility']
    for step in downloads:
        savename = os.path.join(args.paths['path_scripts'], f'download_{step}.sh')
        with open(savename, 'w') as fw:
            write_basetxt_bash(fw, sbid, savename)
            if step == 'selavy':
                write_downloadtxt(args, fw, cat)
            elif step == 'mosaic_images':
                write_downloadtxt(args, fw, img)
            elif step == 'visibility':
                write_downloadtxt(args, fw, vis)

            logger.info('Writing {}'.format(savename))


def format_casa_modeling(args, config):
    # modeling 
    savename = os.path.join(args.paths['path_scripts'], 'casa_model_making.py')
    params = get_modeling_params(config)
    with open(savename, 'w') as fw:  
        write_casa_params(args, params, fw)
        logger.info('Writing {}'.format(savename))

        if config['RESETMODEL']:
            write_casa_reset_vis(fw)

        write_casa_tclean(fw)

        if config['CAL']:
            write_casa_calib(config['SELFCAL'], fw)

        write_casa_subtract_model(fw)
        write_casa_exportfits(params, fw)


def format_casa_imgfast(args, config):
    # short imaging
    savename = os.path.join(args.paths['path_scripts'], 'casa_short_imaging.py')
    params = config['SHORTIMAGE']
    with open(savename, 'w') as fw:  
        write_casa_params(args, params, fw)
        logger.info('Writing {}'.format(savename))

        write_casa_short_imaging_timestamps(config, fw)
        fw.write('    ')
        write_casa_tclean(fw, imagename='imagename_j')
        fw.write('    ')
        write_casa_exportfits(params, fw, imagename='imagename_j')

        

def write_basetxt_bash(fw, sbid, savename):
    logger.debug('write base txt bash for SB%s saving to %s', sbid, savename)
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write(f'# Processing ASKAP data for SB{sbid}' + '\n')
    fw.write(f'# Run this in terminal with "bash {savename}" ' + '\n')
    fw.write('\n')


def write_basetxt_ozstar(fw, sbid, savename, params):
    logger.debug('write base txt ozstar for SB%s saving to %s', sbid, savename)
    logger.debug(params)
    fw.write("#!/bin/bash" + '\n')
    fw.write('#\n')
    fw.write('#SBATCH --time=' + str(params['TIME'])  + '\n')
    fw.write('#SBATCH --job-name=' + str(params['job_name']) + '\n')
    fw.write('#SBATCH --nodes=' + str(params['NODES']) + '\n')
    fw.write('#SBATCH --ntasks-per-node=' + str(params['NTASKS']) + '\n')
    fw.write('#SBATCH --mem-per-cpu=' + str(params['MEM']) + '\n')
    fw.write('#SBATCH --output='+ str(params['output']) + '\n')
    fw.write('#SBATCH --error='+ str(params['error']) + '\n')
    fw.write('#SBATCH --export=all' + '\n')
    fw.write('\n')

def write_endtxt_ozstar(fw, sbid, savename, params):
    logger.debug('write end txt ozstar for SB%s saving to %s', sbid, savename)
    logger.debug(params)
    fw.write('wait' + '\n')
    fw.write('sacct -j $SLURM_JOB_ID --parsable2 --format=' + params['format'] + ' > ' + params['usage'] + '\n')
    fw.write('\n')


def write_moduleload_ozstar(fw, config):
    fw.write('source ' + config['CONDA'] + '\n')
    fw.write('conda activate ' + config['CONDAENV'] + '\n')
    fw.write('\n')


def write_module_unload_ozstar(fw):
    fw.write('conda deactivate' + '\n')
    fw.write('\n')
    

def get_url(access_url, username, password):
    session = requests.Session()
    session.auth = (username, password)

    s = session.get(access_url)
    data_dict = xmltodict.parse(s.content)

    return data_dict['VOTABLE']['RESOURCE'][0]['TABLE']['DATA']['TABLEDATA']['TR'][0]['TD'][1]


def write_downloadtxt(args, fw, cat):
    for i, access_url in enumerate(cat['access_url']):
        url = get_url(access_url, args.username, args.password)
        filename = cat[i]['filename']
        path_file = os.path.join(args.paths['path_data'], filename)
        logger.debug('write download txt %s', path_file)

        text = f'wget -O {path_file} {url} -t 0 -c'

        fw.write("echo " + '\n')
        fw.write(f"echo progress {i+1}/{len(cat)}" + '\n')
        fw.write("echo " + text + '\n')
        fw.write(text + '\n')
        fw.write("sleep 1s" + '\n')
        fw.write('\n')


def prepare_steps_ozstar(args, idx, config, sbid, oname, step='FIXDATA'):
    savename = os.path.join(args.paths['path_scripts'], f'slurm_{step}_beam{idx:02d}.sh')
    params = config['OZSTAR'][step]
    params['job_name'] = step[:3] + f'-{idx:02d}' + f'-{sbid}'
    params['output'] = os.path.join(args.paths['path_logs'], f'slurm_{step}_{oname}.output')
    params['error'] = os.path.join(args.paths['path_logs'], f'slurm_{step}_{oname}.error')
    params['usage'] = os.path.join(args.paths['path_logs'], f'slurm_{step}_{oname}.usage')
    params['format'] = "JobID,JobName,Partition,NodeList,AllocCPUS,State,ExitCode,Elapsed,MaxRSS,MaxVMSize,CPUTime,TotalCPU,Start,End"
    return params, savename


def write_download_vis_txt(args, fw, idx, filename, url):
    path_file = os.path.join(args.paths['path_data'], filename)
    logger.debug('write download txt %s', path_file)

    text = f'wget -O {path_file} {url} -t 0 -c'
    fw.write(f"echo beam{idx:02d}: Download calibrated visibilities data" + '\n')
    fw.write(text + '\n')
    fw.write('\n')


def write_untar_vis_txt(args, fw, idx, filename):
    path_file = os.path.join(args.paths['path_data'], filename)
    logger.debug('write untar txt %s', path_file)
    path_data = args.paths['path_data']
    text = f'tar xvf {path_file} -C {path_data}'
    fw.write(f"echo beam{idx:02d}: Untar the data to measurement sets" + '\n')
    fw.write(text + '\n')
    fw.write('\n')


def write_fixdata_txt(args, fw, idx, filename, prefix=''):
    filename = filename.replace('.tar', '')
    path_file = os.path.join(args.paths['path_data'], filename)
    logger.debug('write fixdata txt %s', path_file)

    text = f'rm -r {path_file}.corrected'
    fw.write('echo Executing: ' + text + '\n')
    fw.write(text + '\n')
    fw.write('\n')

    text = f'askapsoft_rescale {path_file} {path_file}.corrected'
    fw.write(f"echo beam{idx:02d}: Fix the measurement sets flux scaling" + '\n')
    fw.write(prefix + text + '\n')
    fw.write('\n')

    text = f'fix_dir {path_file}.corrected'
    fw.write(f"echo beam{idx:02d}: Fix the measurement sets pointing" + '\n')
    fw.write(prefix + text + '\n')
    fw.write('\n')


def write_run_casa_txt(args, fw, idx, filename, oname, config, mode='modeling', prefix=''):
    filename = filename.replace('.tar', '.corrected')
    path_file = os.path.join(args.paths['path_data'], filename)
    casa = config['CASA']
    logger.debug('write run casa txt %s output name %s', path_file, oname)

    if mode == 'modeling':
        script_name = 'casa_model_making.py'
        log_name = f'casa_MODELING_{oname}.log'
        fw.write(f"echo beam{idx:02d}: Create sky model and subtract..." + '\n')
        fw.write('cd ' + args.paths['path_models'] + ' \n')
    elif mode == 'imaging':
        script_name = 'casa_short_imaging.py'
        log_name = f'casa_IMGFAST_{oname}.log'
        fw.write(f"echo beam{idx:02d}: Create model-subtracted short images..." + '\n')
        fw.write('cd ' + args.paths['path_images'] + ' \n')

    path_log = os.path.join(args.paths['path_logs'], log_name)
    path_script = os.path.join(args.paths['path_scripts'], script_name)
    text = f'{prefix}{casa} --log2term --logfile {path_log} --nogui -c {path_script} {path_file} {oname}'
    fw.write(text + '\n')
    fw.write('\n')


def write_selcand_txt(args, fw, idx, oname, cat, prefix='', affix=''):
    path_script = 'select_candidates'
    path_deepimage = os.path.join(args.paths['path_models'], oname+'.image.tt0.fits') # deep image
    path_catalogue = os.path.join(args.paths['path_data'], cat[0]['filename']) # selavy catalogue
    path_images = args.paths['path_images']
    path_cand = args.paths['path_cand']

    text = f'{prefix}{path_script} --deepimage {path_deepimage} --catalogue {path_catalogue} '\
        f'--folder {path_images} --beam beam{idx:02d} --outdir {path_cand} --name {oname} '\
        f'--ignore-warning --config {args.config}{affix}'

    fw.write("echo beam{:02d}: Select candidates...".format(idx) + '\n')
    fw.write(text + '\n')
    fw.write('\n')


def write_clndata_txt(args, fw, idx):
    fw.write("echo beam{:02d}: Clean data folder...".format(idx) + '\n')

    for key, value in args.paths.items():
        if key == 'path_data':
            fw.write(f'find {value} -type f -name "*beam{idx:02d}*.tar" | xargs -n 1 -t rm' + '\n')
            fw.write(f'find {value} -type d -name "*beam{idx:02d}*.ms" | xargs -n 1 -t rm -r' + '\n')
        elif key == 'path_models' or key == 'path_images':
            fw.write(f'find {value} -type d -name "*beam{idx:02d}*" | xargs -n 1 -t rm -r' + '\n')
            fw.write(f'find {value} -type f -name "*.last" | xargs -n 1 -t rm' + '\n')
    
    fw.write('\n')


def write_casa_params(args, params, fw):
    fw.write('import os' + '\n')
    fw.write('import sys' + '\n')
    fw.write('import numpy as np' + '\n')
    fw.write('\n')
    fw.write('vis = sys.argv[-2]        # visibility path'  + '\n')
    fw.write('imagename = sys.argv[-1]  # recommend in format of SBxxx_beamxxx' + '\n')
    fw.write('print("** NOTICE  ** path passed as arg:", vis, imagename)' + '\n')
    fw.write('\n')
    
    txt = 'print('
    for key, value in params.items():
        if isinstance(value, str):
            fw.write(key.lower() + ' = ' + f'"{value}"' + '\n' )
        else:
            fw.write(key.lower() + ' = ' + f'{value}' + '\n' )
        txt += f'{key.lower()}, '

    fw.write(txt + ')\n')
    fw.write('\n')


def write_casa_reset_vis(fw):
    txt = '''
# reset previous stored model (if there is)
clearcal(vis=vis)
print('Reset all model and corrected data')

'''
    fw.write(txt)


def write_casa_subtract_model(fw):
    txt = '''# subtract model
uvsub(vis=vis)
print('Model-subtraction Finished', imagename)

'''
    fw.write(txt)

def write_casa_tclean(fw, imagename='imagename'):
    txt = f'''tclean(vis=vis, selectdata=True, field='', spw='', timerange=timerange, uvrange=uvrange, 
antenna='', scan='', observation='', intent='', datacolumn=datacolumn, imagename={imagename}, 
imsize=imsize, cell=cell, phasecenter='', stokes=stokes, projection='SIN', specmode='mfs', 
reffreq='', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], 
interpolation='linear', gridder=gridder, facets=facets, wprojplanes=wprojplanes, 
vptable='', aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', 
computepastep=360.0, pblimit=pblimit, normtype='flatnoise', deconvolver=deconvolver,
scales=scales, nterms=nterm, smallscalebias=0.0, restoration=True, restoringbeam=[],
pbcor=pbcor, outlierfile='', weighting=weighting, robust=robust, npixels=0, uvtaper=[],
niter=niter, gain=gain, threshold=threshold, cycleniter=-1, cyclefactor=1.0, 
minpsffraction=0.02, maxpsffraction=0.8, interactive=False, usemask='user', mask='',
pbmask=0.0, savemodel=savemodel, startmodel='', parallel=False)

'''
    fw.write(txt)

def write_casa_exportfits(params, fw, imagename='imagename'):
    nterm = int(params['NTERM'])
    if nterm == 2:
        txt = f'''exportfits(imagename={imagename}+".image.tt0", fitsimage={imagename}+".image.tt0.fits")
'''
    elif nterm == 1:
        txt = f'''exportfits(imagename={imagename}+".image", fitsimage={imagename}+".image.fits")
'''
    else:
        logger.error("uncompatible nterm number --> %s", nterm)
    fw.write(txt)

def write_casa_calib(params, fw):
    fw.write("# gain phase self-calibration" + '\n')
    for key, value in params.items():
        fw.write(key.lower() + ' = ' + f'"{value}"' + '\n' )

    txt = """gaincal(vis=vis, caltable='calib_{}.pcal'.format(imagename), selectdata=True, 
uvrange=uvrange, solint=solint, calmode=calmode)

# apply calibration
applycal(vis=vis, gaintable='calib_{}.pcal'.format(imagename))

"""
    fw.write(txt)


def write_casa_short_imaging_timestamps(config, fw):
    txt = """# open table and read time column
tb.open(vis)
"""

    if config['TIMESTEP'] == 'int':
        txt += '''
from collections import Counter
times = list(Counter(tb.getcol('TIME')).keys())
step = tb.getcol('INTERVAL')[0]
times.sort()
times = np.array(times) - step/2
print("Integration time length", step, "seconds")
    '''
        
    else:
        txt += 'step = ' + str(config['TIMESTEP']) + '\n'
        txt += '''
times = tb.getcol('TIME')
# get the times array, change unit to second
# time in fits is therefore the START time 
times = np.arange(start=np.min(times), stop=np.max(times)+step, step=step)
print("Input time length", step, "seconds, which is", step/60, 'minutes')
'''

    txt += '''
tb.close()

# convert MJD to standard UTC
for j in range(times.shape[0]):
    stime = qa.time(qa.quantity(times[j], 's'), form="ymd")[0]
    etime = qa.time(qa.quantity(times[j]+step, 's'), form='ymd')[0]
    imagename_j = '{}_t{:04d}'.format(imagename, j)
    timerange = '%s~%s' % (stime, etime)

    print('%s, %s' % (os.getcwd(), timerange))

'''
    fw.write(txt)



if __name__ == "__main__":
    _main()

