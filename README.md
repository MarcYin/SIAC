# A sensor invariant Atmospheric Correction (SIAC)
### Feng Yin
### Department of Geography, UCL
### ucfafyi@ucl.ac.uk


[![PyPI version](https://img.shields.io/pypi/v/siac.svg?longCache=true&style=flat)](https://pypi.org/project/SIAC/)
[![conda](https://anaconda.org/f0xy/siac/badges/version.svg?longCache=true&style=flat)](https://anaconda.org/F0XY/siac)
[![py version](https://img.shields.io/pypi/pyversions/siac.svg?longCache=true&style=flat)](https://pypi.org/project/SIAC/)
![Docker Image Size (latest by date)](https://img.shields.io/docker/image-size/marcyin/siac)
[![Documentation Status](https://readthedocs.org/projects/siac/badge/?version=latest)](https://siac.readthedocs.io/en/latest/?badge=latest)
[![Lisence](https://img.shields.io/pypi/l/siac.svg?longCache=true&style=flat)](https://pypi.org/project/SIAC/)
[![DOI](https://zenodo.org/badge/117815245.svg)](https://zenodo.org/badge/latestdoi/117815245)

This atmospheric correction method uses MODIS MCD43 BRDF product to get a coarse resolution simulation of earth surface. A model based on MODIS PSF is built to deal with the scale differences between MODIS and Sentinel 2 / Landsat 8. We uses the ECMWF CAMS prediction as a prior for the atmospheric states, coupling with 6S model to solve for the atmospheric parameters. We do not have topography correction and homogeneouse surface is used without considering the BRDF effects.

## Citation:

Yin, F., Lewis, P. E., & Gómez-Dans, J. L. (2022). Bayesian atmospheric correction over land: Sentinel-2/MSI and Landsat 8/OLI. _EGUsphere_, _2022_, 1–62. doi:10.5194/egusphere-2022-170


## Auxillary data needed (Automatically downloaded by SIAC):
* MCD43 : 
  - 16 days before and 16 days after the Sentinel 2 / Landsat 8 sensing date. 
  - This has been updated to automatically download data from Google Earth Engine (GEE), which is much faster than the preivous way. This means you will need to register to get access to GEE at [here](https://earthengine.google.com).
  - Or you can still use the previous way to download the data by adding the `Gee = False` in the `SIAC_S2` or `SIAC_L8` class, i.e. `SIAC_S2(**kwargs, gee=False)` or `SIAC_L8(**kwargs, gee=False)`.
* ECMWF CAMS Near Real Time prediction: 
  - Time step of 3 hours or 1 hour with the start time of 00:00:00 over the date, and data from 01/01/2015 are mirrored in UK Jasmin server at: https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/
* Global DEM: 
  - Global DEM VRT file built from ASTGTM2 DEM, and most of the DEM over land are mirrored in UK Jasmin server at: https://gws-access.jasmin.ac.uk/public/nceo_ard/DEM_V3/global_dem.vrt
* Emulators: 
  - Emulators for atmospheric path reflectance, total transmittance and single scattering Albedo, and the emulators for Sentinel 2 and Landsat 8 trained with 6S.V2 are packed in the current repository.


## Installation:

You will need to have Gdal and Lightgbm installed and it is suggested to install them with:

- conda:
  ```bash
  conda install -c conda-forge gdal lightgbm
  ```
- mamba:
  ```bash
  mamba install -c conda-forge gdal lightgbm
  ```

Then you can install SIAC:

- Directly from github 

  ```bash
  pip install https://github.com/MarcYin/SIAC/archive/master.zip
  ```

<!-- 
1. Using PyPI

```bash
pip install SIAC
```


3. Using anaconda

```bash
conda install -c f0xy -c conda-forge siac
``` -->


## GEE authenticate:
If you have not used GEE python API before, you will need to authenticate to GEE first after you installed SIAC:

![](https://gist.githubusercontent.com/MarcYin/880d289816b2e8d698f2b6b8d8c514ac/raw/6ef799a4df077598757fe9f7c0fc2cd83e60d372/once.svg)


- In terminal:
  ```bash
  earthengine authenticate --auth_mode=notebook
  ```

- Or in python:
  ```python
  import ee
  ee.Authenticate()
  ```


## Usage:
The typical usage of SIAC for and Landsat 8&9:

- Sentinel 2 
  ```python
  from SIAC import SIAC_S2
  global_dem = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/DEM_V3/global_dem.vrt'
  cams_dir = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/'
  SIAC_S2('/directory/where/you/store/S2/data/', global_dem = global_dem, cams_dir=cams_dir)
  ```

- Landsat 8

  ```python
  from SIAC import SIAC_L8
  global_dem = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/DEM_V3/global_dem.vrt'
  cams_dir = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/'
  SIAC_L8('/directory/where/you/store/L8/data/', global_dem = global_dem, cams_dir=cams_dir) 
  ``` 
- Landsat 9

  ```python
  from SIAC import SIAC_L8
  global_dem = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/DEM_V3/global_dem.vrt'
  cams_dir = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/'
  SIAC_L8('/directory/where/you/store/L9/data/', global_dem = global_dem, cams_dir=cams_dir)
  ```

## Outputs from SIAC_S2

All the outputs from SIAC are specified in the `siac_output.json` within the original S2 L1C folder:

  An example of the [siac_output.json](https://gws-access.jasmin.ac.uk/public/nceo_ard/S2/30/U/UD/S2A_MSIL1C_20170312T114341_N0204_R123_T30UUD_20170312T114341.SAFE/siac_output.json)

The following table specify a list of the outputs from SIAC and their corresponding meanings:

| Abbreviation      | Description | Scale     |  Comments |
| :---           |    :----           |          :--- | :--- |
| `siacLog`        | Siac log file       | N/A  | |
| `toaOvrs`        |  Toa reflectance RGB overviews        | N/A      | |
| `boaOvrs`        |  Surface reflectance RGB overviews        | N/A      | |
| `toaOvrFull`     |  Toa reflectance RGB overview full resolution        | N/A      | |
| `boaOvrFull`    |  Surface reflectance RGB overviews        | N/A      | |
| `viewAngles`     |  View angles for each band       | 0.01    | 2 bands GeoTiff: B1 view azimuth, B2 view zenith  |
| `sunAngles`      |  Sun angles for each band       | 0.01    | 2 bands GeoTiff: B1 sun azimuth, B2 sun zenith  |
| `SurfaceReflectance`        |  Surface reflectance for each band        | 0.0001  | |
| `SurfaceReflectanceUncertainty`        |  Surface reflectance uncertainty for each band   | 0.0001  | |
| `atmoParas`        |  Atmospheric parameters | N/A  | aerosol optical depth[-], total column of water vapour [ $g/cm^2$ ] and total column of ozone [ $cm-atm$ ] |
| `atmoParasUncs`        | Atmospheric parameter uncertainties | N/A  | |
| `Cloud probability`       | Cloud | 0.01  | |
| `Version` | Version of the SIAC software | N/A  | |
| `CleanPixelPercentage` | Clean pixel percentage | N/A  | |
| `ValidPixelPercentage` | Valid pixel percentage | N/A  | |

## Outputs from SIAC_L8

All the outputs from SIAC are within the original L8/L9 L1C folder.


- An example of correction for Landsat 5 for a more detailed demostration of the usage is shown [here](https://github.com/MarcYin/Global-analysis-ready-dataset)


## Examples and Map:

A [page](http://www2.geog.ucl.ac.uk/~ucfafyi/Atmo_Cor/index.html) shows some correction samples.

A [map](http://www2.geog.ucl.ac.uk/~ucfafyi/map) for comparison between TOA and BOA.

## LICENSE
GNU GENERAL PUBLIC LICENSE V3

