# vast-fastdetection

From short model-subtracted FITS images to top transients candidates. 

Some flowcharts can be found in [this Google slides](https://docs.google.com/presentation/d/1ODIjt0YC_LiqUu84r6AsVh4wcZmS4R0KW523N--9PD0/edit?usp=sharing)

## Requirements (recommand)

* numpy                              (1.23.1)
* scipy                              (1.5.2)
* astropy                            (5.1)  
* aplpy                              (2.1.0)   
* matplotlib                         (3.5.2)  
* scikit-image                       (0.19.3)      


## Step-by-step instruction
This instruction assumes the pipeline runs on a computer cluster, so commands like `sbatch <script name>` are used. If it is run on a local machine, copy the line inside the script and run it with python.

**Prepare the scripts for the pipeline**

Change the `loc` to the folder where the pipeline is saved.

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

## Instruction 
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


