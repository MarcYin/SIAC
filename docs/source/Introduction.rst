SIAC
****************************

Introduction
==============

Land surface reflectance is the fundamental variable for the most of earth observation (EO) missions, 
and corrections of the atmospheric disturbs from the cloud, gaseous, aerosol help to get accurate spectral 
description of earth surface. Unlike the previous empirical ways of atmospheric correction, we propose a 
data fusion method for atmospheric correction of satellite images, with an initial attempt to include the 
uncertainty information from different data source. It takes advantage of the high temporal resolution of 
MODIS observations to get BRDF description of the earth surface as the prior information of the earth 
surface property, uses the ECMWF CAMS Near-real-time as the prior information of the atmospheric sates, 
to get optimal estimations of the atmospheric parameters. The code is written in python and we have tested 
it with Sentinel 2, Landsat 8, Landsat 5, Sentinel 3 and MODIS data and it shows SIAC can correct the 
atmospheric effects reasonablely well.

Sentinel 2 and Landsat 8 correction examples
============================================

A `page <http://www2.geog.ucl.ac.uk/~ucfafyi/Atmo_Cor/index.html>`_ shows some correction samples.


A `map <http://www2.geog.ucl.ac.uk/~ucfafyi/map>`_ shows the validation results against AERONET measurements.
If you click points on the scatter plot, it will the comparison between the TOA and BOA refletance and you 
can drag to compare between them and click on the image to have spectral comparison over your clicked pixel as well.
