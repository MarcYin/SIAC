# A sensor invariant Atmospheric Correction (SIAC)
### Feng Yin
### Department of Geography, UCL
### ucfafyi@ucl.ac.uk


[![PyPI version](https://img.shields.io/pypi/v/siac.svg?longCache=true&style=flat)](https://pypi.org/project/SIAC/)
[![conda](https://anaconda.org/f0xy/siac/badges/version.svg?longCache=true&style=flat)](https://anaconda.org/F0XY/siac)
[![py version](https://img.shields.io/pypi/pyversions/siac.svg?longCache=true&style=flat)](https://pypi.org/project/SIAC/)
[![build](https://img.shields.io/travis/MarcYin/SIAC/master.svg?longCache=true&style=flat)](https://travis-ci.org/MarcYin/SIAC/)
[![Documentation Status](https://readthedocs.org/projects/siac/badge/?version=latest)](https://siac.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/MarcYin/SIAC/branch/master/graph/badge.svg?longCache=true&style=flat)](https://codecov.io/gh/MarcYin/SIAC)
[![Coverage Status](https://coveralls.io/repos/github/MarcYin/SIAC/badge.svg?branch=master)](https://coveralls.io/github/MarcYin/SIAC?branch=master)
[![Lisence](https://img.shields.io/pypi/l/siac.svg?longCache=true&style=flat)](https://pypi.org/project/SIAC/)
[![DOI](https://zenodo.org/badge/117815245.svg)](https://zenodo.org/badge/latestdoi/117815245)

This atmospheric correction method uses MODIS MCD43 BRDF product to get a coarse resolution simulation of earth surface. A model based on MODIS PSF is built to deal with the scale differences between MODIS and Sentinel 2 / Landsat 8. We uses the ECMWF CAMS prediction as a prior for the atmospheric states, coupling with 6S model to solve for the atmospheric parameters. We do not have topography correction and homogeneouse surface is used without considering the BRDF effects.


## Data needed:
* MCD43 : 16 days before and 16 days after the Sentinel 2 / Landsat 8 sensing date
* ECMWF CAMS Near Real Time prediction: a time step of 3 hours with the start time of 00:00:00 over the date, and data from 01/04/2015 are mirrored in UCL server at: http://www2.geog.ucl.ac.uk/~ucfafyi/cams/
* Global DEM: Global DEM VRT file built from ASTGTM2 DEM, and most of the DEM over land are mirrored in UCL server at: http://www2.geog.ucl.ac.uk/~ucfafyi/eles/
* Emulators: emulators for atmospheric path reflectance, total transmitance and single scattering Albedo, and the emulators for Sentinel 2, Landsat 8 and MODIS trained with 6S.V2 can be found at: http://www2.geog.ucl.ac.uk/~ucfafyi/emus/

## Installation:

1. Directly from github 

```bash
pip install https://github.com/multiply-org/atmospheric_correction/archive/master.zip
```


2. Using PyPI

```bash
pip install SIAC
```


3. Using anaconda

```bash
conda install -c f0xy -c conda-forge siac
```


To save your time for installing GDAL:

```bash
conda install -c conda-forge gdal>2.1
```



The typical usage for Sentinel 2 and Landsat 8:
```python
from SIAC import SIAC_S2
SIAC_S2('/directory/where/you/store/S2/data/') # this can be either from AWS or Senitinel offical package
```
```python
from SIAC import SIAC_L8                                                                           
SIAC_L8('/directory/where/you/store/L8/data/') 
``` 

An example of correction for Landsat 5 for a more detailed demostration of the usage is shown [here](https://github.com/MarcYin/Global-analysis-ready-dataset)

## Examples and Map:

A [page](http://www2.geog.ucl.ac.uk/~ucfafyi/Atmo_Cor/index.html) shows some correction samples.

A [map](http://www2.geog.ucl.ac.uk/~ucfafyi/map) for comparison between TOA and BOA.

### LICENSE
GNU GENERAL PUBLIC LICENSE V3
