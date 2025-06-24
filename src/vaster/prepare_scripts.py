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

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>, Raghav Girgaonkar <raghav.girgaonkar@gmail.com>"


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
        
        if config['MACHINE'] == 'ozstar':
            format_casa_modeling(args, config)
        elif config['MACHINE'] == 'mortimer':
            format_casa_modeling(args, config)
            format_casa_selfcal(args, config)
        format_casa_imgfast(args, config)
        prepare_downloads(args, sbid, cat, img, vis)

        if config['MACHINE'] == 'bash':
            format_bash(args, config, sbid, vis, cat)
        elif config['MACHINE'] == 'ozstar':
            format_ozstar(args, config, sbid, vis, cat)
        elif config['MACHINE'] == 'mortimer':
            format_mortimer(args, config, sbid, vis, cat)

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
    args.self_config = config_out 
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

    # further filtering of visibility
    mask = [True if filename[:11] == 'scienceData' else False for filename in r['filename']]
    vis = r[mask]

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
            write_imager_txt(args, fw, idx, filename, oname, config, mode='modeling')
            write_imager_txt(args, fw, idx, filename, oname, config, mode='imaging')

            path_log = os.path.join(args.paths['path_logs'], 'bash_SELCAND_{}.log'.format(oname)) # output log files 
            affix = f' > {path_log} 2>&1'
            write_selcand_txt(args, fw, idx, oname, config, cat, affix=affix)
            write_clndata_txt(args, fw, idx, config)
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
                        write_fixdata_txt(args, fw, idx, filename, config, prefix='srun ')
                    elif step == 'MODELING':
                        write_imager_txt(args, fw, idx, filename, oname, config, mode='modeling', prefix='srun ')
                    elif step == 'IMGFAST':
                        write_imager_txt(args, fw, idx, filename, oname, config, mode='imaging', prefix='srun ')
                    elif step == 'SELCAND':
                        write_selcand_txt(args, fw, idx, oname, config, cat, prefix='srun ')
                    elif step == 'CLNDATA':
                        write_clndata_txt(args, fw, idx, config)
                        
                    write_endtxt_ozstar(fw, sbid, savename, params)
                
            logger.info('Writing {}'.format(savename))

def format_mortimer(args, config, sbid, vis, cat):
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
                params, savename = prepare_steps_mortimer(args, idx, config, sbid, oname, step)
                with open(savename, 'w') as fw:
                    write_basetxt_mortimer(fw, sbid, savename, params)
                    if step == 'FIXDATA':
                        write_moduleload_mortimer(fw, config)
                        write_fixdata_txt_mortimer(args, fw, idx, filename, config, prefix='')
                    elif step == 'MODELING':
                        write_moduleload_mortimer(fw, config)
                        write_imager_txt_mortimer(args, fw, idx, filename, oname, config, mode='modeling', prefix='', selfcal=False)
                        if config['PHASESELFCAL']:
                            write_casa_selfcal_cmd(args, fw, idx, filename, oname)
                            write_imager_txt_mortimer(args, fw, idx, filename, oname, config, mode='modeling', prefix='', selfcal=True)
                    elif step == 'IMGFAST':
                        write_moduleload_mortimer(fw, config)
                        if config['PHASESELFCAL']:
                            write_imager_txt_mortimer(args, fw, idx, filename, oname, config, mode='imaging', prefix='', selfcal=True)
                        else:
                            write_imager_txt_mortimer(args, fw, idx, filename, oname, config, mode='imaging', prefix='', selfcal=False)
                    elif step == 'SELCAND':
                        write_moduleload_mortimer(fw, config)
                        write_selcand_txt_mortimer(args, fw, idx, oname, config, cat, prefix='')
                    elif step == 'CLNDATA':
                        write_clndata_txt(args, fw, idx, config)
                        
                    # write_endtxt_ozstar(fw, sbid, savename, params)
                
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
    params = config['CASA']['MODELCLEAN']
    with open(savename, 'w') as fw:  
        write_casa_params(args, params, fw)
        logger.info('Writing {}'.format(savename))

        if config['CASA']['RESETMODEL']:
            write_casa_reset_vis(fw)

        write_casa_tclean(fw)

        if config['CASA']['CAL']:
            write_casa_calib(config['CASA']['SELFCAL'], fw)

        write_casa_subtract_model(fw)
        write_casa_exportfits(params, fw)

def format_casa_selfcal(args, config):
    # Create a selfcal script between consecutive wsclean modelings 
    savename = os.path.join(args.paths['path_scripts'], 'casa_phase_selfcal.py')
    with open(savename, 'w') as fw:  
        write_casa_params_selfcal(fw)
        logger.info('Writing {}'.format(savename))

        write_casa_calib(config['CASA']['SELFCAL'], fw)


def format_casa_imgfast(args, config):
    # short imaging
    savename = os.path.join(args.paths['path_scripts'], 'casa_short_imaging.py')
    params = config['CASA']['SHORTIMAGE']
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

def write_basetxt_mortimer(fw, sbid, savename, params):
    logger.debug('write base txt ozstar for SB%s saving to %s', sbid, savename)
    logger.debug(params)
    fw.write("#!/bin/bash" + '\n')
    fw.write('#\n')
    #fw.write('#SBATCH --time=' + str(params['TIME'])  + '\n')
    fw.write('#SBATCH --job-name=' + str(params['job_name']) + '\n')
    fw.write('#SBATCH --nodes=' + str(params['NODES']) + '\n')
    # fw.write('#SBATCH --ntasks-per-node=' + str(params['NTASKS']) + '\n')
    fw.write('#SBATCH --ntasks=' + str(params['NTASKS']) + '\n')
    if 'PARTITION' in params:
        fw.write('#SBATCH --partition=' + str(params['PARTITION']) + '\n')
    
    #fw.write('#SBATCH --mem-per-cpu=' + str(params['MEM']) + '\n')
    fw.write('#SBATCH --output='+ str(params['output']) + '\n')
    fw.write('#SBATCH --error='+ str(params['error']) + '\n')
    fw.write('#SBATCH --export=all' + '\n')
    fw.write('#SBATCH --mail-type=all' + '\n')
    fw.write('#SBATCH --mail-user='+ str(params['email']) + '\n')
    fw.write('\n')

def write_endtxt_ozstar(fw, sbid, savename, params):
    logger.debug('write end txt ozstar for SB%s saving to %s', sbid, savename)
    logger.debug(params)
    fw.write('wait' + '\n')
    fw.write('sleep 10' + '\n')
    fw.write('sacct -j $SLURM_JOB_ID --parsable2 --format=' + params['format'] + ' > ' + params['usage'] + '\n')
    fw.write('\n')

def write_moduleload_mortimer(fw, config):
    fw.write('source ' + config['CONDA'] + '\n')
    fw.write('conda activate ' + config['CONDAENV'] + '\n')
    fw.write('source ' + config['CASAPATH'] + '\n')
    fw.write('\n')


def write_virtual_env_enable(fw, config):
    # write line by line
    for line in config['VIRTUAL_ENV_ENABLE']:
        fw.write(line + '\n')
    fw.write('\n')


def write_virtual_env_disable(fw, config):
    fw.write(config['VIRTUAL_ENV_DISABLE'] + '\n')
    fw.write('\n')


def write_singularity_load(fw):
    fw.write('module load singularity' + '\n')
    fw.write('export SINGULARITY_BINDPATH=$PWD' + '\n')
    fw.write('\n')

def write_singularity_load_mortimer(fw,config):
    fw.write('export SINGULARITY_BINDPATH=' + config['SINGULARITY_BINDPATH'] + '\n')
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

def prepare_steps_mortimer(args, idx, config, sbid, oname, step='FIXDATA'):
    '''Mortimer specific prepare steps'''
    savename = os.path.join(args.paths['path_scripts'], f'slurm_{step}_beam{idx:02d}.sh')
    params = config['MORTIMER'][step]
    params['email'] = config['EMAIL']
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


def write_fixdata_txt(args, fw, idx, filename, config, prefix=''):
    filename = filename.replace('.tar', '')
    path_file = os.path.join(args.paths['path_data'], filename)
    logger.debug('write fixdata txt %s', path_file)

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_enable(fw, config)
        
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

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_disable(fw, config)

def write_fixdata_txt_mortimer(args, fw, idx, filename, config, prefix=''):
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

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_disable(fw, config)


def write_imager_txt(args, fw, idx, filename, oname, config, mode='modeling', prefix=''):
    filename = filename.replace('.tar', '.corrected')
    path_file = os.path.join(args.paths['path_data'], filename)
    imager = config['IMAGER']
    logger.debug('write imaging txt %s output name %s', path_file, oname)

    if imager == 'casa':
        text = write_casa_cmd(args, fw, idx, path_file, oname, mode, prefix, imager)
    elif imager == 'wsclean':
        text = write_wsclean_cmd(args, fw, idx, path_file, oname, config, mode, prefix)
    elif imager == 'flint':
        logger.warning('FLINT HAS NOT BEEN IMPLEMENTED YET - WILL SWITCH TO WSCLEAN')
        text = write_wsclean_cmd(args, fw, idx, path_file, oname, config, mode, prefix)
    else:
        logger.warning('UNIDENTIFIED IMAGER "%s" - WILL SWITCH TO WSCLEAN', config['IMAGER'])
        text = write_wsclean_cmd(args, fw, idx, path_file, oname, config, mode, prefix)

    fw.write('echo ' + text + '\n')
    fw.write(text + '\n')
    fw.write('\n')

def write_imager_txt_mortimer(args, fw, idx, filename, oname, config, mode='modeling', prefix='', selfcal=False):
    filename = filename.replace('.tar', '.corrected')
    path_file = os.path.join(args.paths['path_data'], filename)
    imager = config['IMAGER']
    logger.debug('write imaging txt %s output name %s', path_file, oname)

    if imager == 'casa':
        text = write_casa_cmd(args, fw, idx, path_file, oname, mode, prefix, imager)
    elif imager == 'wsclean':
        text = write_wsclean_cmd_mortimer(args, fw, idx, path_file, oname, config, mode, prefix, selfcal)
    elif imager == 'flint':
        logger.warning('FLINT HAS NOT BEEN IMPLEMENTED YET - WILL SWITCH TO WSCLEAN')
        text = write_wsclean_cmd_mortimer(args, fw, idx, path_file, oname, config, mode, prefix, selfcal)
    else:
        logger.warning('UNIDENTIFIED IMAGER "%s" - WILL SWITCH TO WSCLEAN', config['IMAGER'])
        text = write_wsclean_cmd_mortimer(args, fw, idx, path_file, oname, config, mode, prefix, selfcal)

    fw.write('echo ' + text + '\n')
    fw.write(text + '\n')
    fw.write('\n')


def write_selcand_txt(args, fw, idx, oname, config, cat, prefix='', affix=''):
    path_script = 'select_candidates'

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_enable(fw, config)

    if config['IMAGER'] == 'casa':
        path_deepimage = os.path.join(args.paths['path_models'], oname+'.image.tt0.fits') # deep image
    elif config['IMAGER'] == 'wsclean':
        path_deepimage = os.path.join(args.paths['path_models'], oname+'-MFS-image.fits') # deep image
    else:
        path_deepimage = os.path.join(args.paths['path_models'], oname+'*image*fits') # deep image

    if config['SOURCE_FINDER'] == 'selavy':
        path_catalogue = os.path.join(args.paths['path_data'], cat[0]['filename']) # selavy catalogue
    elif config['SOURCE_FINDER'] == 'aegean':
        # running aegean to search for candidates
        fw.write("echo beam{:02d}: Running aegean to produce deep image catalogues...".format(idx) + '\n')
        fw.write(f"cd {args.paths['path_models']}" + '\n')
        fw.write(f"{prefix}aegean {path_deepimage} --cores 1 --save" + '\n')
        tablename = path_deepimage.replace('fits', 'cat.fits')
        fw.write(f"{prefix}aegean {path_deepimage} --cores 1 --table {tablename}" + '\n')
        fw.write('\n')
        path_catalogue = path_deepimage.replace('fits', 'cat_comp.fits') # selavy catalogue
    else:
        path_catalogue = os.path.join(args.paths['path_data'], cat[0]['filename']) # selavy catalogue
    
    path_images = args.paths['path_images']
    path_cand = args.paths['path_cand']

    text = f'{prefix}{path_script} --deepimage {path_deepimage} --catalogue {path_catalogue} '\
        f'--folder {path_images} --beam beam{idx:02d} --outdir {path_cand} --name {oname} '\
        f'--ignore-warning --config {args.self_config}{affix}'

    fw.write("echo beam{:02d}: Select candidates...".format(idx) + '\n')
    fw.write(text + '\n')
    fw.write('\n')

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_disable(fw, config)

def write_selcand_txt_mortimer(args, fw, idx, oname, config, cat, prefix='', affix=''):
    path_script = 'select_candidates'

    # write_moduleload_mortimer(fw, config)

    if config['IMAGER'] == 'casa':
        path_deepimage = os.path.join(args.paths['path_models'], oname+'.image.tt0.fits') # deep image
    elif config['IMAGER'] == 'wsclean':
        path_deepimage = os.path.join(args.paths['path_models'], oname+'-MFS-image.fits') # deep image
    else:
        path_deepimage = os.path.join(args.paths['path_models'], oname+'*image*fits') # deep image

    if config['SOURCE_FINDER'] == 'selavy':
        path_catalogue = os.path.join(args.paths['path_data'], cat[0]['filename']) # selavy catalogue
    elif config['SOURCE_FINDER'] == 'aegean':
        # running aegean to search for candidates
        fw.write("echo beam{:02d}: Running aegean to produce deep image catalogues...".format(idx) + '\n')
        fw.write(f"cd {args.paths['path_models']}" + '\n')
        fw.write(f"{prefix}aegean {path_deepimage} --cores 1 --save" + '\n')
        tablename = path_deepimage.replace('fits', 'cat.fits')
        fw.write(f"{prefix}aegean {path_deepimage} --cores 1 --table {tablename}" + '\n')
        fw.write('\n')
        path_catalogue = path_deepimage.replace('fits', 'cat_comp.fits') # selavy catalogue
    else:
        path_catalogue = os.path.join(args.paths['path_data'], cat[0]['filename']) # selavy catalogue
    
    path_images = args.paths['path_images']
    path_cand = args.paths['path_cand']

    text = f'{prefix}{path_script} --deepimage {path_deepimage} --catalogue {path_catalogue} '\
        f'--folder {path_images} --beam beam{idx:02d} --outdir {path_cand} --name {oname} '\
        f'--ignore-warning --config {args.self_config}{affix}'

    fw.write("echo beam{:02d}: Select candidates...".format(idx) + '\n')
    fw.write(text + '\n')
    fw.write('\n')

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_disable(fw, config)


def write_clndata_txt(args, fw, idx, config):
    fw.write("echo beam{:02d}: Clean data folder...".format(idx) + '\n')

    for key, value in args.paths.items():
        if key == 'path_data':
            fw.write(f'find {value} -type f -name "*beam{idx:02d}*.tar" | xargs -n 1 -t rm' + '\n')
            fw.write(f'find {value} -type d -name "*beam{idx:02d}*.ms" | xargs -n 1 -t rm -r' + '\n')

        elif key == 'path_models':
            if config['IMAGER'] == 'casa':
                fw.write(f'find {value} -type d -name "*beam{idx:02d}*" | xargs -n 1 -t rm -r' + '\n')
                fw.write(f'find {value} -type f -name "*.last" | xargs -n 1 -t rm' + '\n')
            elif config['IMAGER'] == 'wsclean':
                fw.write(f'find {value} -type f -name "*beam{idx:02d}*" ! -name "*MFS*" | xargs -n 1 -t rm' + '\n')
            else:
                continue

        elif key == 'path_images':
            if config['IMAGER'] == 'casa':
                fw.write(f'find {value} -type d -name "*beam{idx:02d}*" | xargs -n 1 -t rm -r' + '\n')
                fw.write(f'find {value} -type f -name "*.last" | xargs -n 1 -t rm' + '\n')
            elif config['IMAGER'] == 'wsclean':
                fw.write(f'find {value} -type f -name "*beam{idx:02d}*psf.fits" | xargs -n 1 -t rm' + '\n')
                fw.write(f'find {value} -type f -name "*beam{idx:02d}*dirty.fits" | xargs -n 1 -t rm' + '\n')
            else:
                continue

        elif key == 'path_cand':
            fw.write(f'find {value} -type f -name "*beam{idx:02d}*peak.fits" | xargs -n 1 -t rm' + '\n')
            fw.write(f'find {value} -type f -name "*beam{idx:02d}*chisquare.fits" | xargs -n 1 -t rm' + '\n')
            fw.write(f'find {value} -type f -name "*beam{idx:02d}*std.fits" | xargs -n 1 -t rm' + '\n')
            
    
    fw.write('\n')

######## ####### ####### #######
####### CASA starts #######
######## ####### ####### #######
def write_casa_cmd(args, fw, idx, path_file, oname, mode, prefix, imager):
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
    text = f'{prefix}{imager} --log2term --logfile {path_log} --nogui -c {path_script} {path_file} {oname}'
    return text 

def write_casa_selfcal_cmd(args, fw, idx, filename, oname):
    script_name = 'casa_phase_selfcal.py'
    log_name = f'casa_PHASESELFCAL_{oname}.log'
    fw.write(f"echo beam{idx:02d}: Perform phase only self calibration..." + '\n')
    fw.write('cd ' + args.paths['path_models'] + ' \n')

    filename = filename.replace('.tar', '.corrected')
    path_file = os.path.join(args.paths['path_data'], filename)
    path_log = os.path.join(args.paths['path_logs'], log_name)
    path_script = os.path.join(args.paths['path_scripts'], script_name)
    text = f'casa --log2term --logfile {path_log} --nogui -c {path_script} {path_file} {oname}'

    fw.write('echo ' + text + '\n')
    fw.write(text + '\n')
    fw.write('\n')
    # return text

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

def write_casa_params_selfcal(fw):
    fw.write('import os' + '\n')
    fw.write('import sys' + '\n')
    fw.write('import numpy as np' + '\n')
    fw.write('\n')
    fw.write('vis = sys.argv[-2]        # visibility path'  + '\n')
    fw.write('imagename = sys.argv[-1]  # recommend in format of SBxxx_beamxxx' + '\n')
    fw.write('print("** NOTICE  ** path passed as arg:", vis, imagename)' + '\n')
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

######## ####### ####### #######
####### CASA ends #######
######## ####### ####### #######
####### WSCLEAN starts #######
######## ####### ####### #######

def write_wsclean_cmd(args, fw, idx, path_file, oname, config, mode, prefix):
    if config['SINGULARITY'] is True:
        write_singularity_load(fw)
        run_wsclean = 'singularity exec ' + config['WSCLEAN_PATH'] + ' wsclean'
    else:
        run_wsclean = 'wsclean'

    if mode == 'modeling':
        fw.write(f"echo beam{idx:02d}: Create sky model and subtract..." + '\n')
        fw.write('cd ' + args.paths['path_models'] + ' \n')
        params = config['WSCLEAN']['MODELCLEAN']
        text = write_wsclean_params(params)

    elif mode == 'imaging':
        write_intervals_out(args, fw, config, path_file, oname)
        fw.write(f"echo beam{idx:02d}: Create model-subtracted short images..." + '\n')
        fw.write('cd ' + args.paths['path_images'] + ' \n')
        params = config['WSCLEAN']['SHORTIMAGE']
        text = write_wsclean_params(params) + '-intervals-out $intervals_out '
    
    text = f'{prefix}{run_wsclean} {text}-name {oname} {path_file}'
    return text 

def write_wsclean_cmd_mortimer(args, fw, idx, path_file, oname, config, mode, prefix, selfcal):
    if config['SINGULARITY'] is True:
        if mode == 'modeling':
            if selfcal == False:
                write_singularity_load_mortimer(fw,config)
        elif mode == 'imaging':
            write_singularity_load_mortimer(fw,config)
        run_wsclean = 'singularity exec ' + config['WSCLEAN_PATH'] + ' wsclean'
    else:
        run_wsclean = 'wsclean'

    if mode == 'modeling':
        fw.write(f"echo beam{idx:02d}: Create sky model and subtract..." + '\n')
        fw.write('cd ' + args.paths['path_models'] + ' \n')
        params = config['WSCLEAN']['MODELCLEAN']
        if selfcal == True:
            # print('selfcal true')
            params_temp = params.copy()
            params_temp['DATA_COLUMN'] = 'CORRECTED_DATA'
            text = write_wsclean_params(params_temp)
        else:
            # print('selfcal false')
            text = write_wsclean_params(params)

    elif mode == 'imaging':
        write_intervals_out_mortimer(args, fw, config, path_file, oname)
        fw.write(f"echo beam{idx:02d}: Create model-subtracted short images..." + '\n')
        fw.write('cd ' + args.paths['path_images'] + ' \n')
        params = config['WSCLEAN']['SHORTIMAGE']
        if selfcal == True:
            params_temp = params.copy()
            params_temp['DATA_COLUMN'] = 'CORRECTED_DATA'
            text = write_wsclean_params(params_temp) + '-intervals-out $intervals_out '
        else:
            text = write_wsclean_params(params) + '-intervals-out $intervals_out '
    
    text = f'{prefix}{run_wsclean} {text}-name {oname} {path_file}'
    return text 


def write_wsclean_params(params):
    txt = ''
    for key, value in params.items():
        if value is True:
            txt += '-' + key.lower().replace('_', '-') + ' '
        elif value is False:
            txt += ''
        else:
            txt += '-' + key.lower().replace('_', '-') + ' ' + f'{value}' + ' '

    logger.debug('******* wsclean imaging parameters *******')
    logger.debug(txt)
    return txt


def write_intervals_out(args, fw, config, path_file, oname):
    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_enable(fw, config)

    savename = os.path.join(args.paths['path_data'], oname + '_measurements.txt')
    fw.write(f'intervals=($(check_measurements {path_file} --config {args.self_config} --savename {savename}))' + '\n')
    fw.write(r"intervals_out=${intervals[-1]}" + '\n')
    fw.write('\n')

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_disable(fw, config)

def write_intervals_out_mortimer(args, fw, config, path_file, oname):
    
    # write_moduleload_mortimer(fw, config)

    savename = os.path.join(args.paths['path_data'], oname + '_measurements.txt')
    fw.write(f'intervals=($(check_measurements {path_file} --config {args.self_config} --savename {savename}))' + '\n')
    fw.write(r"intervals_out=${intervals[-1]}" + '\n')
    fw.write('\n')

    if config['VIRTUAL_ENV'] is True:
        write_virtual_env_disable(fw, config)
    


if __name__ == "__main__":
    _main()

