# vast-fastdetection

From short model-subtracted FITS images to top transients candidates. 

## Example usage

**Load a bunch of short images**

Save the location of a series of short images into `imagelist`. Note the images should be in a correct order (e.g., with time ascending). 

Generate a `Cube` class for the following processing. 

```python
imagelist = glob.glob('/folder/to/your/images/*.fits')

cube = Cube(imagelist)
```

### Generate a significance cube

Convolve a spatial kernel with each image, to smooth the noise level and improve the detection. 

```python
cube.icube(ktype='gaussian', nx=19, ny=19)
```

The kernel type can be modified through `ktype='gaussian'` (a 2D gaussian kernel), `ktype='psf'` (dirty beam)

The kernel size can be modified through `nx` and `ny` (pixel values)

The kernel HWFM will be automatically calculated through fits header information (i.e., the synthesised beam size). 

The generate smoothed cube is saved in `cube.sigcube`

### Select transients candidates through a matched filter

Build a Filter using the generated cube

```python
f = Filter(cube.sigcube)
```

Do a Gaussian (or other kernel) smooth on time axis 

```python
f.fmap("gaussian", width=1)
```

Find the local maximum and get the sky position - need a fits file to provide wcs (can be any of FITS images from the `imagelist`)

```python
f.local_max(imagename=imagelist[0], sigma=3, min_distance=120)
```

The detection threshold can be changed using `sigma`. 

`min_distance` is the minimum distance (pixel) of neighouring local maximum



