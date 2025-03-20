# VASTER - ASKAP Intra-Observation Transient Search Pipeline

VASTER is a novel pipeline designed for the efficient detection and characterization of intra-observation transients in data from the Australian Square Kilometre Array Pathfinder (ASKAP) telescope. 
This pipeline streamlines the entire process, from data retrieval to candidate selection, providing valuable insights into transient astrophysical phenomena.

See details in the paper [Wang Y. et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523.5661W/abstract)

Some flowcharts can also be found in [this Google slides](https://docs.google.com/presentation/d/1ODIjt0YC_LiqUu84r6AsVh4wcZmS4R0KW523N--9PD0/edit?usp=sharing)

Example outputs: [SB11676](https://unisydneyedu-my.sharepoint.com/:f:/g/personal/ywan3191_uni_sydney_edu_au/Em-RqhrRKdBGmnvJzd1CU4UBZyFnSh4E6jvBzu9icCux2g?e=e6FAow)

**(The making of a PDF containing detailed installation and running instructions with updated documentaton is in progress)**

## Installation

### Recommended General Installation Instructions (Applicable to Clusters or Personal Machines)

1. To install VASTER, it is recommanded to install it in a virtual environment (such as with conda or miniconda) with python==3.9. Download and install miniconda from here 
[miniconda](https://docs.anaconda.com/miniconda/install/) 

2. For ease of adding conda and casa to your current path variable, create two scripts: `.activate_conda` and `.activate_casa`. In `.activate_conda` add the line:
```
eval "$(/path/to/miniconda3/bin/conda shell.bash hook)"
```
and in `.activate_casa` add the line:
```
export PATH="/path/to/casa/bin:$PATH"
```

3. These two scripts can then be used to activate conda and casa for running slurm jobs by adding the lines,
```
source /path/to/.activate_conda
source /path/to/.activate_casa
```
to the slurm batch file or to activate conda or casa in general in your terminal session.

4. Once conda and casa are activated, create a vaster conda environment with python 3.9 as follows (the name of the environment can be different),
```
conda create --name vaster python=3.9
```
and activate it,
```
conda activate vaster
```

5. Download or git clone this repo, 
```
git clone https://github.com/askap-vast/vast-fastdetection/tree/main
```

6. We now have to install the required libraries and the VASTER pipeline, this can be installed via 
```
pip install .
```
An editable version (for development) can be installed as follows,
```
pip install -e .
``` 
7. In case of dependency issues for `aplpy`, first comment out aplpy (by adding a `#` as prefix) in the `requirements.txt` file and run 
```
pip install -r requirements
```
once the rest of the libraries have been installed, install aplpy using conda-forge as follows,
```
conda install -c conda-forge aplpy
```
This takes care of all dependency issues that may arise.

The VASTER pipeline can then be installed by
```
pip install .
```
or an editable version
```
pip install -e .
```

8. Test the installation by running 
```
prepare_scripts -h
```
This should show a help message for the `prepare_scripts` utility.

## Package Requirements

* python                             >=3.9
* numpy                              >=1.23.1
* scipy                              
* astropy                            >=5.3.3  
* aplpy                              
* matplotlib                         <3.6
* scikit-image                          
* astroquery
* xmltodict                          https://github.com/conda-forge/xmltodict-feedstock
* python-casacore 
* pandas
* pillow                             >=10.0.0

## Features

### Workflow 

1. **Automated Data Retrieval:** VASTER can automatically download ASKAP data from [CASDA](https://data.csiro.au/domain/casdaObservation?redirected=true) based on given SBID, simplifying the initial data acquisition process.

2. **Model Creation:** The pipeline generates a deep sky model from calibrated visibilities to serve as a reference for data analysis.

3. **Model Subtraction:** VASTER performs model-subtraction to isolate deviations from the deep sky model, potentially revealing transient sources.

4. **Short Image Creation:** Users can specify timescales (e.g., 10 seconds or 15 minutes), and VASTER creates short residual images for further analysis.

5. **Statistical Map Analysis:** The pipeline generates statistical maps (chi-square map, peak map, standard deviation map, and Gaussian map) based on the short residual images, aiding in the identification of transient objects.

6. **Transient Selection**: VAST can select transient candidates based on the generated statistical maps.

###  Outputs
For each transient candidate, VASTER generates the following outputs:  
* candidate catalogues (csv)
* model image cutouts (png)
* animations created from model-subtracted snapshot images (gif)


## Step-by-Step instructions (Cluster)
Instructions for running VASTER in a slurm supercomputer system - tested on OzStar and UWM's Mortimer clusters.

### Prerequisites (OzStar)

1. To create an account on OzStar refer to [New Account Creation](https://supercomputing.swin.edu.au/account-management/new_account_request). 
One can then login by:
```
ssh -XY <username>@ozstar.swin.edu.au
```

2. Once account is created, you would need to join (with approval) a project in order to get access to and analyze data. Each project has a project code (such as oz330).
A request to join project can be sent here [Account Management System](https://supercomputing.swin.edu.au/account-management/).

3. You would also need an OPAL account, which can be created here [OPAL account creation](https://opal.atnf.csiro.au/register). **Note: This is a required step to run VASTER on any machine.**

### Preparing Scripts 

1. (Optional) In your personal directory, create an empty folder named after the SBID(s) you want to analyze. For eg: 

```
mkdir -p SB50210
```

2. Copy a sample `config.yml` file from [config file](https://github.com/askap-vast/vast-fastdetection/blob/main/prepare/config.yml) into your directory. (Not in SBXXXX/). Update any parameters in this file (such as the machine, inaging parameters etc).

3. We can now prepare the scripts as follows. If the config file is already in your directory, 
```
prepare_scripts <sbid>
```
else
```
prepare_scripts <sbid> --config /path/to/config.yml
```
For a complete list of parameters available for this steps, one can run
```
prepare_scripts -h
```
In runtime, this script may ask you for your OPAL username and password. To bypass this for every time you run this script for different fields, you can set the environment variables `$OPAL_USER` and `$OPAL_PWD` as your username and password respectively. This option is also available in the `config.yml` file.

4. This will create a new folder named "SBxxxx" under current directory. 
Within the new folder, it will generate the following sub directories
* data/ (to store visibilities)
* models/ (to store time-averaged model images)
* images/ (to store model-subtracted short images)
* candidates/ (to store final candidates including csv and png)
* scripts/ (to store scripts that will be used)
* logfiles/ (to store logfiles)
The created scripts for analysis can be found in the SBXXXX/scripts directory.

5. Note: Certain slurm-based clusters may not have the `sacct` utility running or available to monitor jobs. In such cases, the corresponding line in the slurm scripts should be commented out.

### Downloading Visibility Data
To download visibility data, we will use the `prepare_data` utility.
1. To get help on various parameters of this utility run, 
```
prepare_data -h
```
2. By default all 36 beams of the SBID will be downloaded, but one can specify specific beams by:
```
prepare_data <sbid> -b 0 1 2 --untar
```
this will download beam00, beam 01 and beam02 visibility data in the SBXXXX/data/ folder and untar it.

### Running Scripts
Once your data is downloaded the `submit_slurm_jobs` utility can be used to submit multiple slurm jobs corresponding to the different steps of VASTER at once. 
1. Help on this utility can be accessed through:
```
submit_slurm_jobs -h
``` 

2. To slurm scripts for all steps run
```
submit_slurm_jobs <sbid> -b 0 1 2 
```
This will submit jobs for beams 0, 1 and 2.

3. To submit certain steps of the job,
```
submit_slurm_jobs <sbid> -b 0 1 2 --steps FIXDATA MODELING IMGFAST
```
This will only submit jobs for the first 3 steps of the pipeline. 

4. VASTER has 5 main steps: FIXDATA MODELING IMGFAST SELCAND and CLNDATA, **these steps have to be run in this particular order.** 
each step can be run individually for any beam as follows:

**Rescale and fix the data**
```
sbatch /your/output/folder/scripts/slurm_FIXDATA_beamXX.sh
```
**Deep Modeling**
```
sbatch /your/output/folder/scripts/slurm_MODELING_beamXX.sh
```
A .fits image in the form of `SBID_beamXX.image.tt0.fits` should appear in the models folder.

**Short Timescale Imaging**
```
sbatch /your/output/folder/scripts/slurm_IMGFAST_beamXX.sh
```
Short images will be saved in the images folder.
You can check the progress in the `.output` file from the logfiles folder if running on a cluster.

**Candidate Selection**
```
sbatch /your/output/folder/scripts/slurm_SELCAND_beamXX.sh
```
Chi-square and peak fits images of each beam will be saved in the candidates folder. Lightcurve, deep image and snapshot animation of the candidates (if any) are also saved.

**Final candidate list**
```python
python /your/folder/to/pipeline/get_overall_table.py <SBID>
```
Change the `base_folder` to the pathname of the candidates folder.
`base_url` is used when the lightcurve, deep image and snapshot are saved to an online server.

A `SBID.csv` file will be produced at the end.

## Step-by-step instructions: (Personal Machines)
When running VASTER on personal machines and not slurm-based supercomputing systems the following steps can be followed:

1. Before creating scripts (see Running Scripts section), change ```machine=bash``` in config file. 
2. Once required data is retrieved (see Downloading Visibility Data section), the different steps of VASTER can be run as follows:

**Rescale and fix the data**
```
askapsoft_rescale /path/to/SBXXXXX/data/<filename>.beamXX_averaged_cal.leakage.ms /path/to/SBXXXXX/data/<filename>.beam00_averaged_cal.leakage.ms.corrected
fix_dir /path/to/SBXXXXX/data/<filename>.beam00_averaged_cal.leakage.ms.corrected
```
**Deep Modeling**
```
cd /path/to/SBXXXXX/models
casa --log2term --logfile /path/to/SBXXXXX/logfiles/casa_MODELING_SBXXXXX_beamXX.log --nogui -c /path/to/SBXXXXX/scripts/casa_model_making.py /path/to/SBXXXXX/data/<filename>.beamXX_averaged_cal.leakage.ms.corrected SBXXXXX_beamXX
```
**Short Timescale Imaging**
```
cd /path/to/SBXXXXX/images
casa --log2term --logfile /path/to/SBXXXXX/logfiles/casa_IMGFAST_SBXXXXX_beamXX.log --nogui -c /path/to/SBXXXXX/scripts/casa_short_imaging.py /path/to/SBXXXXX/data/<filename>.beamXX_averaged_cal.leakage.ms.corrected SBXXXXX_beamXX
```
**Candidate Selection**
```
select_candidates --deepimage /path/to/SBXXXXX/models/SBXXXX_beamXX.image.tt0.fits --catalogue /path/to/SBXXXXX/data/selavy-image.i.<filename>.SBXXXXX.cont.taylor.0.restored.conv.components.xml --folder /path/to/SBXXXXX/images --beam beamXX --outdir /path/to/SBXXXXX/candidates --name SBXXXX_beamXX --ignore-warning --config /path/to/config.yml
```
3. The above steps can be run manually or one can use the `bash_PROCESSING_beamXX.sh` script that is created (after `prepare_scripts`) to run them sequentially as follows:
```
bash bash_PROCESSING_beamXX.sh
```

## Parameter customisation (OUTDATED)
**Short timescale imaging setting**

10-second images are generated using the task `tclean` from CASA. Paramters can be adjusted in `/vast-fastdetection/imaging/short_imaging.py` to accommodate for the scientific goal. Relevant parameters may include:

`imsize` controls the size of the image in unit of pixels and `cell` sets the angular size of one pixel. Increase `cell` when imaging a larger image to reduce runtime.

`uvrange` sets the uv-range of data to be imaged. Greater value represents more compact object.

`gridder` and `wprojplanes` determine the type of gridder used and the amount of w-values employed for W-projection. `gridder = widefield` and `wprojplanes = -1` accounts for the w-term and generates spatially accurate image but requires longer runtime. `gridder = standard` and `wprojplanes = 1` ignores w-projection and produces inaccurate image, especially when the source is away from the beam center, but the imaging is roughly 5 times faster.

Refer to CASA documentation for further details on `tclean`.

**Candidate selection threshold**

The sigma-level of the chi-square map and peak map can be adjusted in `/vast-fastdetection/run_all.py`.
The limit of candidates plotted is also set in that python code.



## Detailed instruction for candidates selection 

The `select_candidates.py` scripts can automatically build a cube, generate final significance maps, select candidates and generate final products using optimised parameters. 

If you are interested in modifying some parameters, please see below. 

```python
from vastcube.cube import Cube, Filter
from vastcube.plot import Candidates
```

### Generate a (significance) cube

**Load a bunch of short images**

Save the location of a series of short images into `imagelist`. Note the images should be in a correct order (e.g., with time ascending). 

Create a `Cube` class for the following processing. 

```python
imagelist = glob.glob('/folder/to/your/images/*.fits')

cube = Cube(imagelist)
```

**Remove bad images**

Remove images with rms larger than two times median RMS level. 

```python
cube.remove_bad_images()
```

**Generate a Cube**

The default way in `run_cube.py` is to form a cube directly 

```python
cube.icube()
```

If can convolve a spatial kernel with each image, to smooth the noise level and improve the detection. 

```python
cube.icube(ktype='gaussian', nx=19, ny=19)
```

The kernel type can be modified through `ktype='gaussian'` (a 2D gaussian kernel), `ktype='psf'` (dirty beam)

The kernel size can be modified through `nx` and `ny` (pixel values)

The kernel HWFM will be automatically calculated through fits header information (i.e., the synthesised beam size). 

Note this smooth process will take much more time (~40min). 

The generate (smoothed) cube is saved in `cube.sigcube`

### Select transients candidates through a matched filter

**Build a Filter using the generated cube**

```python
f = Filter(cube.sigcube)
```

Do a Gaussian (or other kernel) smooth on time axis 

```python
f.fmap("gaussian", width=4)
```

* Chisquare map - "chisquare"
* Gaussian map - "gaussian"
* Peak map - "peak"
* Standard Deviation map - "std"

**Save the output map to fits file**

```python
f.to_fits(fitsname, imagename=imagelist[0])
```

**Find local maximum (candidates)**

```python
c = Candidates(chisq_map, peak_map, std_map)
```

Or include Gaussian map using argument `gaussian_map=<gaussian_map_name>`

```python
c.local_max(min_distance=30, sigma=5)
```

The detection threshold can be changed using `sigma` (in log space). 

`min_distance` is the minimum distance (pixel) of neighouring local maximum

**Remove potential artefacts**

We need deep source catalogue information to removew potential artifacts - currently it supports `aegean` format table. But it can also change to `selavy` format table. 

```python
c.select_candidates(deepcatalogue)
```

**Save final results to a csv/vot table**

```python
c.save_csvtable(tablename, savevot=True)
```


