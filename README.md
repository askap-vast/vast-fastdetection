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


## Simple usage

```python
python run_cube.py <deep_source_catalog> <folder_short_fits_images> <beam_number> <output_name>
```

The `run_cube.py` scripts can automatically build a cube, generate final significance maps, select candidates and generate final products using optimised parameters. 

If you are interested in modifying some parameters, please see below. 

## Instruction 

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


