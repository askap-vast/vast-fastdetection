#!/usr/bin/env python3
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

import os
import sys
import getpass
import requests
import xmltodict


sbid = sys.argv[-1]  # number only
# path = sys.argv[-1] # output parent location 
path = os.path.join(os.getcwd(), 'SB'+sbid) # create a SB directory in the current location 

loc = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # code location
print("Processing observation", 'SB'+sbid)
print("Saving outptus to", path)
print("Using code in", loc)

# nodes = ['purley-x86-cpu{:02d}'.format(i) for i in range(1, 8) if i not in [1, 4]] + ['hw-x86-cpu{:02d}'.format(j) for j in range(11, 11) if j not in [4, 5, 9]] 
# exclude_nodes = 'purley-x86-cpu[02,08],hw-x86-cpu[01-15]' # hw-x86 is extremely slow!!!
# nodelist = (nodes * 100)[:36]

# print('nodelist', nodelist)

############################
# Build file saving system structure 
############################
path_data = os.path.join(path, 'data') # saving visibilities, selavy catalogues 
path_models = os.path.join(path, 'models') 
path_images = os.path.join(path, 'images')
path_cand = os.path.join(path, 'candidates')
path_fits = os.path.join(path, 'fitsfiles')

path_scripts = os.path.join(path, 'scripts')
path_logs = os.path.join(path, 'logfiles')


def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
        print('Create new directory', dir)
    else:
        print('Directory {} exists.'.format(dir))


create_dir(path)
create_dir(path_data)
create_dir(path_models)
create_dir(path_images)
create_dir(path_cand)
create_dir(path_scripts)
create_dir(path_logs)
create_dir(path_fits)


############################
# Find visibilities, selavy catalogues and fits images from CASDA
############################
tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
job = tap.launch_job_async("SELECT * FROM ivoa.obscore WHERE obs_id='{}' AND (dataproduct_type='{}' OR dataproduct_subtype='{}' OR dataproduct_subtype='{}') ".format(
    sbid, 'visibility', 'catalogue.continuum.component', 'cont.restored.t0'))

r = job.get_results()

vis = r[r['dataproduct_type'] == 'visibility']
cat = r[r['dataproduct_subtype'] == 'catalogue.continuum.component']
img = r[r['dataproduct_subtype'] == 'cont.restored.t0']

vis.sort('filename')

print('Found {} visibilities'.format(len(vis)))
print(vis[['filename', 'access_url']])
print('')
print('Found {} selavy catalogues'.format(len(cat)))
print(cat[['filename', 'access_url']])
print('')
print('Found {} images'.format(len(img)))
print(img[['filename', 'access_url']])
print('')

if len(vis) != 36:
    print("")
    print('Number of visibilities is not 36. Exit with error. ')
    print('You can download the visibility manually. ')
    sys.exit()


############################
# Get download urls
############################
# username = 'wym20131028@gmail.com'
# print('OPAL username:', username)
username = input("Enter OPAL username: ")
password = getpass.getpass(str("Enter OPAL password: "))


def get_url(access_url):
    session = requests.Session()
    session.auth = (username, password)

    s = session.get(access_url)
    data_dict = xmltodict.parse(s.content)

    return data_dict['VOTABLE']['RESOURCE'][0]['TABLE']['DATA']['TABLEDATA']['TR'][0]['TD'][1]


############################
# Generate bash_GETDATA_beam??.sh
############################

for idx in range(36):
    savename = os.path.join(path_scripts, 'bash_GETDATA_beam{:02d}.sh'.format(idx))

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')
        fw.write('# Generate automatically from a python script' + '\n')
        fw.write('# Download and untar visibility for SB{} beam{:02d}'.format(sbid, idx) + '\n')
        fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
            savename) + '\n')
        fw.write('\n\n\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        url = get_url(vis[idx]['access_url'])
        filename = vis[idx]['filename']
        path_file = os.path.join(path_data, filename)

        text = 'wget -O {} {} -t 0'.format(path_file, url)
        fw.write("echo " + '\n')
#        fw.write(text + '\n')
        fw.write(text + ' -c' + '\n')
        fw.write('\n')

#        text = 'tar xvf {} -C {}'.format(path_file, path_data)
#        fw.write("echo " + '\n')
#        fw.write(text + '\n')
#        fw.write('\n')

    print('Writing {}'.format(savename))



############################
# Generate bash_CHECKDATA.sh
############################

# savename = os.path.join(path_scripts, 'bash_CHECKDATA.sh')

# with open(savename, 'w') as fw:
#     fw.write("#!/bin/bash" + '\n')
#     fw.write('\n')
#     fw.write('# Generate automatically from a python script' + '\n')
#     fw.write('# Download visibility for SB{} '.format(sbid) + '\n')
#     fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
#         savename) + '\n')
#     fw.write('\n\n\n')

#     for idx in range(36):

#         if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
#             print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

#         url = get_url(vis[idx]['access_url'])
#         filename = vis[idx]['filename']
#         path_file = os.path.join(path_data, filename)

#         text = 'wget -O {} {} -t 0 -c'.format(path_file, url)
#         fw.write("echo " + '\n')
#         fw.write("echo progress {}/{}".format(idx+1, len(vis)) + '\n')
#         fw.write(text + '\n')
#         fw.write("sleep 1s" + '\n')
#         fw.write('\n')

# print('Writing {}'.format(savename))



############################
# Generate download_selavy.sh
############################
savename = os.path.join(path_scripts, 'download_selavy.sh')

with open(savename, 'w') as fw:
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write('# Download selavy components catalogue for SB{}'.format(sbid) + '\n')
    fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
        savename) + '\n')
    fw.write('\n\n\n')

    for i, access_url in enumerate(cat['access_url']):
        filename = cat[i]['filename']
        try:
            url = get_url(access_url)
            path_file = os.path.join(path_data, filename)
            cat_ind = i

            text = 'wget -O {} {} -t 0'.format(path_file, url)

            fw.write("echo " + '\n')
            fw.write("echo progress {}/{}".format(i+1, len(cat)) + '\n')
            fw.write("echo " + text + '\n')
            fw.write(text + '\n')
            fw.write("sleep 1s" + '\n')
            fw.write('\n')
        except:
            print(i, filename, access_url, 'doesnt work')

print('Writing {}'.format(savename))


############################
# Generate download_mosaic_images.sh
############################
savename = os.path.join(path_scripts, 'download_mosaic_images.sh')

with open(savename, 'w') as fw:
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write('# Download mosaiced fits images for SB{}'.format(sbid) + '\n')
    fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
        savename) + '\n')
    fw.write('\n\n\n')

    for i, access_url in enumerate(img['access_url']):
        filename = img[i]['filename']
        try:
            url = get_url(access_url)
            path_file = os.path.join(path_data, filename)

            text = 'wget -O {} {} -t 0'.format(path_file, url)

            fw.write("echo " + '\n')
            fw.write("echo progress {}/{}".format(i+1, len(img)) + '\n')
            fw.write("echo " + text + '\n')
            fw.write(text + '\n')
            fw.write("sleep 1s" + '\n')
            fw.write('\n')
        except:
            print(i, filename, access_url, 'doesnt work')

print('Writing {}'.format(savename))



############################
# Generate slurm_FIXDATA_??.sh
############################


# def FIXDATA(affix, filename):
#     string = """
#     #!/bin/bash 

#     #SBATCH --partition=all-x86-cpu
#     #SBATCH --time=1:00:00
#     #SBATCH --job-name=FIX-{:02d}
#     #SBATCH --nodes=1
#     #SBATCH --ntasks-per-node=1
#     #SBATCH --mem=10gb
#     #SBATCH --output={}.output
#     #SBATCH --error={}.error
#     #SBATCH --export=all

#     module use /home/app/modulefiles
#     module load casacore/cpu-py3.6.5-3.1.0

#     time -p python {} {} {}

#     time -p python {} {}
#     """.format(
#         idx, 
#         os.path.join(path_logs, 'slurm_FIXDATA_'+affix), 
#         os.path.join(path_logs, 'slurm_FIXDATA_'+affix), 
#         os.path.join(loc, 'askapsoft_rescale.py'), 
#         filename, 
#         filename+'.corrected', 
#         os.path.join(loc, 'fix_dir.py'), 

#     )

#     return string 


for idx in range(36):
    savename = os.path.join(path_scripts, 'slurm_FIXDATA_beam{:02d}.sh'.format(idx))
    affix = 'SB{}_beam{:02d}'.format(sbid, idx)

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')

        # fw.write('#SBATCH --partition=purley-cpu' + '\n')
        fw.write('#SBATCH --time=1:00:00' + '\n')
        fw.write('#SBATCH --job-name=FIX-{:02d}'.format(idx) + '\n')
        fw.write('#SBATCH --nodes=1' + '\n')
        fw.write('#SBATCH --ntasks-per-node=1' + '\n')
        
#        if 'exclude_nodes' in locals():
#            fw.write('#SBATCH --exclude={}'.format(exclude_nodes) + '\n')
#        elif 'nodelist' in locals():
#            fw.write('#SBATCH --nodelist={}'.format(nodelist[idx]) + '\n')
        
        fw.write('#SBATCH --mem=50gb' + '\n')
        fw.write('#SBATCH --output='+os.path.join(path_logs, 'slurm_FIXDATA_{}.output'.format(affix)) + '\n')
        fw.write('#SBATCH --error='+os.path.join(path_logs, 'slurm_FIXDATA_{}.error'.format(affix)) + '\n')
        fw.write('#SBATCH --export=all' + '\n')
        fw.write('\n')

        # fw.write('module use /home/app/modulefiles' + '\n')
        # fw.write('module load casacore/cpu-py3.6.5-3.1.0' + '\n')
        fw.write('source ~/.activate_conda' + '\n')
        fw.write('conda activate vaster' + '\n')

        fw.write('\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        filename = vis[idx]['filename'][:-4]
        path_file = os.path.join(path_data, filename)

        text = 'srun time -p python {} {} {}'.format(
            os.path.join(loc, 'tools', 'askapsoft_rescale.py'), 
            path_file, path_file+'.corrected'
            )        
        fw.write(text + '\n')
        fw.write('\n')

        text = 'srun time -p python {} {}'.format(
            os.path.join(loc, 'tools', 'fix_dir.py'), 
            path_file+'.corrected'
            )    
        fw.write(text + '\n')
        fw.write('\n')

        fw.write('conda deactivate' + '\n')

    print('Writing {}'.format(savename))



############################
# Generate slurm_MODELING_??.sh
############################

for idx in range(36):
    savename = os.path.join(path_scripts, 'slurm_MODELING_beam{:02d}.sh'.format(idx))
    affix = 'SB{}_beam{:02d}'.format(sbid, idx)

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')

        fw.write('#SBATCH --time=50:00:00' + '\n')
        fw.write('#SBATCH --job-name=MOD-{:02d}'.format(idx) + '\n')
        fw.write('#SBATCH --nodes=1' + '\n')
        fw.write('#SBATCH --ntasks-per-node=1' + '\n')
        fw.write('#SBATCH --mem=200gb' + '\n')
        fw.write('#SBATCH --output='+os.path.join(path_logs, 'slurm_MODELING_{}.output'.format(affix)) + '\n')
        fw.write('#SBATCH --error='+os.path.join(path_logs, 'slurm_MODELING_{}.error'.format(affix)) + '\n')
        fw.write('#SBATCH --export=all' + '\n')
        fw.write('\n')

        fw.write('\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        filename = vis[idx]['filename'][:-4] + '.corrected'
        path_file = os.path.join(path_data, filename)

        text = 'srun time casa --logfile {} --nogui -c {} {} {}'.format(
            os.path.join(path_logs, 'casa_MODELING_{}.log'.format(affix)), 
            os.path.join(loc, 'imaging', 'model_making.py'), 
            path_file, 
            affix
        )
        fw.write(text + '\n')
        fw.write('\n')

    print('Writing {}'.format(savename))


############################
# Generate slurm_IMGFAST_??.sh
############################

for idx in range(36):
    savename = os.path.join(path_scripts, 'slurm_IMGFAST_beam{:02d}.sh'.format(idx))
    affix = 'SB{}_beam{:02d}'.format(sbid, idx)

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')

        # fw.write('#SBATCH --partition=all-x86-cpu' + '\n')
        fw.write('#SBATCH --time=30:00:00' + '\n')
        fw.write('#SBATCH --job-name=IMG-{:02d}'.format(idx) + '\n')
        fw.write('#SBATCH --nodes=1' + '\n')
        fw.write('#SBATCH --ntasks-per-node=1' + '\n')

        # if 'exclude_nodes' in locals():
        #     fw.write('#SBATCH --exclude={}'.format(exclude_nodes) + '\n')
        # elif 'nodelist' in locals():
        #     fw.write('#SBATCH --nodelist={}'.format(nodelist[idx]) + '\n') 
        
        fw.write('#SBATCH --mem=50gb' + '\n')
        fw.write('#SBATCH --output='+os.path.join(path_logs, 'slurm_IMGFAST_{}.output'.format(affix)) + '\n')
        fw.write('#SBATCH --error='+os.path.join(path_logs, 'slurm_IMGFAST_{}.error'.format(affix)) + '\n')
        fw.write('#SBATCH --export=all' + '\n')
        fw.write('\n')

        # fw.write('module use /home/app/modulefiles' + '\n')
        # fw.write('module load casa/5.0.0-218.el6' + '\n')
        # fw.write('module load python/cpu-3.6.5' + '\n')
        fw.write('\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        filename = vis[idx]['filename'][:-4] + '.corrected'
        path_file = os.path.join(path_data, filename)

        text = 'srun time casa --logfile {} --nogui -c {} {} {} {}'.format(
            os.path.join(path_logs, 'casa_IMGFAST_{}.log'.format(affix)), 
            os.path.join(loc, 'imaging', 'short_imaging.py'), 
            path_file, 
            affix, 
            900
        )
        fw.write(text + '\n')
        fw.write('\n')

    print('Writing {}'.format(savename))



############################
# Generate slurm_SELCAND_??.sh
############################

for idx in range(36):
    savename = os.path.join(path_scripts, 'slurm_SELCAND_beam{:02d}.sh'.format(idx))
    affix = 'SB{}_beam{:02d}'.format(sbid, idx)

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')

        # fw.write('#SBATCH --partition=all-x86-cpu' + '\n')
        fw.write('#SBATCH --time=10:00:00' + '\n')
        fw.write('#SBATCH --job-name=SEL-{:02d}'.format(idx) + '\n')
        fw.write('#SBATCH --nodes=1' + '\n')
        fw.write('#SBATCH --ntasks-per-node=1' + '\n')
        
        # if 'exclude_nodes' in locals():
        #     fw.write('#SBATCH --exclude={}'.format(exclude_nodes) + '\n')
        # elif 'nodelist' in locals():
        #     fw.write('#SBATCH --nodelist={}'.format(nodelist[idx]) + '\n') 

        fw.write('#SBATCH --mem=50gb' + '\n')
        fw.write('#SBATCH --output='+os.path.join(path_logs, 'slurm_SELCAND_{}.output'.format(affix)) + '\n')
        fw.write('#SBATCH --error='+os.path.join(path_logs, 'slurm_SELCAND_{}.error'.format(affix)) + '\n')
        fw.write('#SBATCH --export=all' + '\n')
        fw.write('\n')

        # fw.write('module use /home/app/modulefiles' + '\n')
        # fw.write('module load python/cpu-3.7.4' + '\n')
        fw.write('source ~/.activate_conda' + '\n')
        fw.write('conda activate vaster' + '\n')
        fw.write('\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        text = 'srun time python {} {} {} {} {} {} {}'.format(
            os.path.join(loc, 'select_candidates.py'), # scripts
            # os.path.join(path_models, affix+'.image.tt0.fits'), # deep image
            os.path.join(path_data, img[0]['filename']), # deep image
            os.path.join(path_data, cat[cat_ind]['filename']), # selavy catalogue
            path_images, # short images location
            'beam{:02d}'.format(idx), # beam number
            path_cand, # output directory 
            affix# affix
        )
        fw.write(text + '\n')
        fw.write('\n')

        fw.write('conda deactivate')

    print('Writing {}'.format(savename))


############################
# Generate slurm_CLNDATA_??.sh
############################

for idx in range(36):
    savename = os.path.join(path_scripts, 'slurm_CLNDATA_beam{:02d}.sh'.format(idx))
    affix = 'SB{}_beam{:02d}'.format(sbid, idx)

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')

        # fw.write('#SBATCH --partition=purley-cpu' + '\n')
        fw.write('#SBATCH --time=1:00:00' + '\n')
        fw.write('#SBATCH --job-name=CLN-{:02d}'.format(idx) + '\n')
        fw.write('#SBATCH --nodes=1' + '\n')
        fw.write('#SBATCH --ntasks-per-node=1' + '\n')
        fw.write('#SBATCH --mem=50gb' + '\n')
        fw.write('#SBATCH --output='+os.path.join(path_logs, 'slurm_CLNDATA_{}.output'.format(affix)) + '\n')
        fw.write('#SBATCH --error='+os.path.join(path_logs, 'slurm_CLNDATA_{}.error'.format(affix)) + '\n')
        fw.write('#SBATCH --export=all' + '\n')
        fw.write('\n')
        fw.write('\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        filename = vis[idx]['filename'][:-4]
        path_file = os.path.join(path_data, filename)

        fw.write('rm {}.tar'.format(path_file) + '\n')
        fw.write('rm -r {}'.format(path_file) + '\n')
        fw.write('\n')

        fw.write('mv {} {}'.format(os.path.join(path_models, '*beam{:02d}*.fits'.format(idx)), path_fits) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_images, '*beam{:02d}*.fits'.format(idx)), path_fits) + '\n')
        fw.write('rm -r {}/*beam{:02d}*'.format(path_models, idx) + '\n')
        fw.write('rm -r {}/*beam{:02d}*'.format(path_images, idx) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_fits, '*beam{:02d}*image*fits'.format(idx)), path_models) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_fits, '*beam{:02d}*fits'.format(idx)), path_images) + '\n')
        fw.write('\n')

    print('Writing {}'.format(savename))