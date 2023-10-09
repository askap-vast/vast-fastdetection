# VASTER - ASKAP Intra-Observation Transient Search Pipeline

VASTER is a noval pipeline designed for the efficient detection and characterization of intra-observation transients in data from the Australian Square Kilometre Array Pathfinder (ASKAP) telescope. 
This pipeline streamlines the entire process, from data retrieval to candidate selection, providing valuable insights into transient astrophysical phenomena.

See details in the paper [Wang Y. et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523.5661W/abstract)

Some flowcharts can also be found in [this Google slides](https://docs.google.com/presentation/d/1ODIjt0YC_LiqUu84r6AsVh4wcZmS4R0KW523N--9PD0/edit?usp=sharing)

Example outputs: [SB11676](https://unisydneyedu-my.sharepoint.com/:f:/g/personal/ywan3191_uni_sydney_edu_au/Em-RqhrRKdBGmnvJzd1CU4UBZyFnSh4E6jvBzu9icCux2g?e=e6FAow)

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

## Installation

```
git clone git@github.com:askap-vast/vast-fastdetection.git
```

### Package Requirements

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

## Quick Usage

Running in slurm supercomputer system - only tested for Shanghai SKA regional prototype
```
bash run_everything.sh <sbid>
```
It will then automatically download data from CASDA and then submit a bunch of sbatch jobs 

Rnning in bash
```
bash run_everything_for_bash.sh <sbid>
```

It will create a new folder named "SBxxxx" under current directory. 
Within the new folder, it will generate the following sub directories
* data/ (to store visibilities)
* models/ (to store time-averaged model images)
* images/ (to store model-subtracted short images)
* candidates/ (to store final candidates including csv and png)
* fitsfiles/ (intermeidate folder)
* scripts/ (to store scripts that will be used)
* logfiles/ (to store logfiles)


## Step-by-step insturction (bash shell)

1. prepare data structure 
```
python tools/get_everything_ready_for_bash.py <sbid number> <path to save results>
```
It will create data structure under defined path, with sub directories mentioned above. 


2. download visibilities
```
bash <path>/scripts/download_selavy.sh
bash <path>/scripts/bash_CHECKDATA.sh 
```

**All of the following steps are stored in <path>/scripts/bash_PROCESSING_beam??.sh**, i.e.
```
bash <path>/scripts/bash_PROCESSING_beam??.sh
```
Read the following if you want to do it manually. 

3. untar visibilities
```
tar xvf <path>/data/scienceData.....tar
```

4. fix askap soft flux scale and beam positions
```
python tools/askapsoft_rescale.py <input visibilities> <output visibilities name>
python tools/fix_dir.py <visibilities name>
```

5. create sky model and subtract
```
casa --logfile <path>/logfiles/<name of logfile> -c imaging/model_making.py <input visibilities> <output name>
```

6. create model-subtracted short images
```
casa --logfile <path>/logfiles/<name of logfile> -c imaging/short_imaging.py <input visibilities> <output name>
```

7. candidates selection
```python
python select_candidates.py <time-averaged model image> <deep_source_selavy_catalog> <folder_short_fits_images> <beam_number> <output folder> <output_name>
```


## Step-by-step instruction (cluster)
This instruction assumes the pipeline runs on a computer cluster, so commands like `sbatch <script name>` are used. If it is run on a local machine, copy the line inside the script and run it with python.

**Prepare the scripts for the pipeline**

If running the pipeline on a cluster, adjust the runtime of `MODELING`, `IMGFAST`, and `SELCAND` as appropriate.
```python
python /your/folder/to/pipeline/tools/get_everything_ready.py <SBID> <output_folder>
```

**Download the visibility data**
Download the data one by one for 36 beams.
```bash
bash /your/output/folder/scripts/bash_GETDATA_beamXX.sh
```
Download selavy components and mosaiced fits images.
```bash
bash /your/output/folder/scripts/download_selavy.sh
bash /your/output/folder/scripts/download_mosaic_images.sh
```

**Untar the data**
```bash
cd /your/output/folder/data
tar xvf scienceData.XXX.YYY.ZZZ.beamXX_averaged_cal_leakage.ms.tar
```
Repeat this step for all 36 beams.

A folder of the same name without the tar extension can be seen in the data folder after the process completes.

**Rescale and fix the data**
```bash
sbatch /your/output/folder/scripts/slurm_FIXDATA_beamXX.sh
```

**Deep modelling**
```bash
sbatch /your/output/folder/scripts/slurm_MODELING_beamXX.sh
```
A .fits image in the form of `SBID_beamXX.image.tt0.fits` should appear in the models folder.

**Short timescale imaging**
```bash
sbatch /your/output/folder/scripts/slurm_IMGFAST_beamXX.sh
```
Short images will be saved in the images folder.
You can check the progress in the `.output` file from the logfiles folder if running on a cluster.

**Candidate selection**
```bash
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

## Paramter customisation
**Short timescle imaging setting**

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

**generate a cube**

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


