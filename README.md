# Sentinel 2 and Landsat 8 Atmospheric correction 
### Feng Yin
### Department of Geography, UCL
### ucfafyi@ucl.ac.uk

This atmospheric correction method uses MODIS MCD43 BRDF product to get a coarse resolution simulation of earth surface. A model based on MODIS PSF is built to deal with the scale differences between MODIS and Sentinel 2 / Landsat 8. We uses the ECMWF CAMS prediction as a prior for the atmospheric states, coupling with 6S model to solve for the atmospheric parameters. We do not have a proper cloud mask at the moment and no topography correction as well. Homogeneouse surface is used without considering the BRDF effects.

## Data needed:
* MCD43 : 16 days before and 16 days after the Sentinel 2 / Landsat 8 sensing date
* ECMWF CAMS Near Real Time prediction: a time step of 3 hours with the start time of 00:00:00 over the date
* global dem: Global DEM VRT file built from ASTGTM2 DEM, and a bash script under eles/ can be used to generate with the individual files.
* emus (`multiply_atmospheric_corection/emus`): emulators for the 6S and wv restrival, can be found at: http://www2.geog.ucl.ac.uk/~ucfafyi/emus/

### LICENSE
GNU GENERAL PUBLIC LICENSE V3
