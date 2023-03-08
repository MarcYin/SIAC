# SIAC nadir BRDF-adjusted reflectance (NBAR) and uncertainty calculation

This is a description of the SIAC nadir BRDF-adjusted reflectance (NBAR) product. It closely follows [Roy et al. 2016](https://doi.org/10.1016/j.rse.2016.01.023) paper. The correction is down by normalising the viewing geometry and the solar geometry to the nadir geometry with a given solar zenith angle (SZA). For generating NBAR for time series of images, the SZA can be constant or variable (changing the SZA over the year). 


# Implementation

The implementation of the NBAR algorithm is done in Python. The code is available on GitHub within the [SIAC](https://github.com/MarcYin/SIAC) repository. The NBAR algorithm is within the `get_bar.py` file. The main function is `create_nbar`:

```python
def create_nbar(s2_file_dir, nbar_sza='atan2', logger=None, 
                mosaic_start_date=None, mosaic_end_date=None, mosaic_hour = None, 
                Gee = True, use_VIIRS = False, vnp43_folder=None, temporal_window = 16):
    """
    Create NBAR from Sentinel-2 data.
    
    Solar zenith angle `nbar_sza`:
    
    It can be either the mean of the SZA from S2 ('use_s2') 
    Or the SZA at the subsolar point ('atan2', default) from https://doi.org/10.1016/j.renene.2021.03.047.
    If mosaic over a period of time, use the mean SZA for the whole period, calculated
    from the SZA with the subsolar point ('temporal_average_sza').
    Or to any user defined sza (float number between 0-60 is suggested).

    Args:
        s2_file_dir (str): path to the Sentinel-2 data.
        nbar_sza (str, optional): nbar sza. Defaults to 'atan2'.
        logger (logging, optional): logger for function. Defaults to None.
        mosaic_start_date (datetime, optional): mosaic starting date. Defaults to None.
        mosaic_end_date (datetime, optional): mosaic ending date. Defaults to None.
        mosaic_hour (float, optional): float hour between 0-24. Defaults to None.
        Gee (bool, optional): Whether to use GEE or not. Defaults to True.
        use_VIIRS (bool, optional): whether to use VIIRS. Defaults to False.
        vnp43_folder (str, optional): path for saving VIIRS data. Defaults to None.
        temporal_window (int, optional): days before and after the obervation date. Defaults to 16.
    return:
        None: the outputs are saved in the same folder as the S2 data.
    """    
```

Table of sources of BRDF information:

| Source | Description|
| --- | --- |
| use_VIIRS | VIIRS BRDF product from NASA |
| Gee | MODIS BRDF products hosted by Google Earth Engine |
|  |

Table for different SZA options:

| nbar_sza | Description | fname_postfix |
| --- | --- | --- |
| use_s2 | Use the mean of the SZA from S2 | _nbar_sza_s2 |
| atan2 | Use the SZA at the subsolar point from https://doi.org/10.1016/j.renene.2021.03.047 | _nbar_sza_avg |
| temporal_average_sza | Use the mean SZA for the whole period, calculated from the SZA with the subsolar point | _nbar_sza_temporal_avg |
| {float} | Use any user defined sza (float number between 0-60 is suggested) | \_nbar\_sza_[%02d] |
|  ||





# Theoretical background


## MODIS/VIRRIS BRDF kernel

The MODIS/VIIRS BRDF information is derived with the RossThick-LiSparseReciprocal (RTLSR) kernels model from  atmospherically corrected multiangular reflectance observations, detail description is available [here](https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf). The BRDF kernel is calculated at 7 MODIS bands [1, 2, 3, 4, 5, 6, 7], at wavelengths of 645, 858, 469, 555, 1240, 1640, 2130 nm, respectively. The officail product proivides three BRDF weighting parameters: isotropic ($f_{iso}$), volume ($f_{vol}$) and geometric ($f_{geo}$) for all the 7 MODIS bands. With the BRDF kernel parameters $K_{iso}$, $K_{vol}$ and $K_{geo}$, angular reflectance $\hat{r}^{BRDF}(\theta_v, \theta_s)$ can be calculated as follows:

$$
\begin{align}
\hat{r}^{BRDF}(\theta_v, \theta_s) = f_{iso} + f_{vol}K_{vol}(\theta_v, \theta_s) + f_{geo}K_{geo}(\theta_v, \theta_s)
\end{align}
$$

where $\theta_v$ is the viewing geometry and $\theta_s$ is the solar geometry. 



## Interpolation from MODIS wavelength to Sentinel-2 wavelength

As there is no BRDF information available from Sentinel-2 (S2), the MODIS/VIIRS BRDF information is used. However, the S2 reflectance is calculated at 10 Sentinel-2 bands [2, 3, 4, 5, 6, 7, 8, 8A, 11, 12], at wavelengths of [492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4], respectively. Therefore, the MODIS BRDF kernel needs to be interpolated to the S2 wavelengths. The interpolation is done by linear interpolating the kernel weighting parameters of two MODIS bands closest to each S2 bands:

The interpolation slope is calculated by the wavelengths of the two MODIS bands closest ($\lambda_{left}^{MODIS}$, $\lambda_{right}^{MODIS}$) to each S2 band ($\lambda_{s2}$):

$$
\begin{align}
m = \frac{\lambda_{s2}- \lambda_{left}^{MODIS}}{\lambda_{right}^{MODIS} - \lambda_{left}^{MODIS}}
\end{align}
$$


The kernel weighting parameters for S2 band ($f_{[iso, vol, geo]}^{S2}$) is then calculates from the two MODIS bands ($f_{left_{[iso, vol, geo]}}^{MODIS}$, $f_{right_{[iso, vol, geo]}}^{MODIS}$) as follows:

$$
\begin{align}
f_{[iso, vol, geo]}^{S2} = f_{left_{[iso, vol, geo]}}^{MODIS} + m(f_{right_{[iso, vol, geo]}}^{MODIS} - f_{left_{[iso, vol, geo]}}^{MODIS})
\end{align}
$$

The uncertainty in the interpolated kernel weighting parameters ($\sigma_{f_{[iso, vol, geo]}^{S2}}$) is calculated as follows:

$$
\begin{align}
\sigma_{f_{[iso, vol, geo]}^{S2}} = \sqrt{m^2 \sigma_{f_{right_{[iso, vol, geo]}}^{MODIS}}^2 + (1-m)^2 \sigma_{f_{left_{[iso, vol, geo]}}^{MODIS}}^2}
\end{align}
$$

where $\sigma_{f_{right_{[iso, vol, geo]}}^{MODIS}}$ and $\sigma_{f_{left_{[iso, vol, geo]}}^{MODIS}}$ are the uncertainties in the kernel weighting parameters of the two MODIS bands closest to each S2 band. These uncertainties in the BRDF kernel weighting parameters are calculated using the [SIAC](https://gmd.copernicus.org/articles/15/7933/2022/) gap filled BRDF kernel weighting parameters and their uncertainty. 

For S2 bands that are outside the MODIS bands, the closest MODIS band is used and the uncertainty is set by normilal increasing the uncertainty of the closest MODIS band by 20% for each band outside the MODIS bands.


## Simulation of anglular reflectance for S2

After the interpolation of the BRDF kernel weighting parameters, the angular reflectance $\hat{r}^{BRDF}(\theta_v, \theta_s)$ can be calculated for each S2 band. 

$$
\begin{align}
\hat{r}^{BRDF}(\theta_v, \theta_s) = f_{iso}^{S2} + f_{vol}^{S2}K_{vol}(\theta_v, \theta_s) + f_{geo}^{S2}K_{geo}(\theta_v, \theta_s)
\end{align}
$$

where $\theta_v$ is the viewing geometry and $\theta_s$ is the solar geometry.

And the uncertainty in the simulate surface reflectance ($\sigma_{\hat{r}^{BRDF}(\theta_v, \theta_s)}$) is calculated from uncertainty of the interpolated S2 kernel weight parameters ($\sigma_{f_{[iso, vol, geo]}^{S2}}$):

$$
\begin{align}
\sigma_{\hat{r}^{BRDF}(\theta_v, \theta_s)} = \sqrt{(\frac{\partial \hat{r}^{BRDF}}{\partial f_{iso}^{S2}} \sigma_{f_{iso}^{S2}})^2 + (\frac{\partial \hat{r}^{BRDF}}{\partial f_{vol}^{S2}}\sigma_{f_{vol}^{S2}})^2 + (\frac{\partial \hat{r}^{BRDF}}{\partial f_{geo}^{S2}} \sigma_{f_{geo}^{S2}})^2}
\end{align}
$$


where $\frac{\partial \hat{r}^{BRDF}}{\partial f_{iso}} = 1$, $\frac{\partial \hat{r}^{BRDF}}{\partial f_{vol}} = K_{vol}(\theta_v, \theta_s)$ and $\frac{\partial \hat{r}^{BRDF}}{\partial f_{geo}} = K_{geo}(\theta_v, \theta_s)$ are the partial derivatives of the BRDF kernel weighting parameters with respect to the BRDF kernel weighting parameters ($f_{iso}$, $f_{vol}$ and $f_{geo}$).


So, the uncertainty in the BRDF reflectance ($\sigma_{\hat{r}^{BRDF}(\theta_v, \theta_s)}$) is calculated as follows:

$$
\begin{align}
\sigma_{\hat{r}^{BRDF}(\theta_v, \theta_s)} = \sqrt{(\sigma_{f_{iso}^{S2}})^2 +  K_{vol}(\theta_v, \theta_s)^2 \sigma_{f_{vol}^{S2}}^2 + K_{geo}(\theta_v, \theta_s)^2 \sigma_{f_{geo}^{S2}}^2}
\end{align}
$$


## NBAR calculation

Assuming that a SZA ($\theta_s^{nbar}$) is given, the NBAR for a given band of wavelength $\lambda$ is calculated as follows the [Roy et al. 2016](https://doi.org/10.1016/j.rse.2016.01.023) paper. Since the azimuths angles will not affect the c factor if the view zenith angle is 0, so we can ignore the azimuths angles in the following calculation:

$$
\begin{align}
c= \frac{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}
\end{align}
$$

$$
\begin{align}
r_{nbar} = c\times r_{s2}
\end{align}
$$

Here $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$ is the BRDF reflectance at the S2 view and solar geometry. $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ is the BRDF reflectance at the nadir viewing geometry and the given SZA ($\theta_s^{nbar}$). The BRDF reflectance is calculated using the MODIS/VIIRS BRDF information. The $c$ is the correction factor to normalise the S2 reflectance to the given SZA. The $r_{s2}$ is the surface reflectance of S2. The $r_{nbar}$ is the NBAR reflectance.



### Considerations on SZA

This would cause issues when the surface property is not homogeneous, as within the MODIS 500meterx500meter pixel the BRDF shape could be vastly different at the scale of S2 10metersx10meters pixel. Therefore, it is recommended to not extrapolate the NBAR SZA too far from the initial S2 SZA. If the application requires images from an individual date and can handle the variation in SZA, you may choose to normalise the SZA to the original S2 reported SZA or the SZA calculated from the image latitude and longitude. However, if the application requires mosaic from time series of images, it is recommended to use the mean SZA calculated over time series of images.



## NBAR uncertainty calculation

### Uncerainty in $c$ factor
The uncertainty calculation for NBAR is not given in the [Roy et al. 2016](https://doi.org/10.1016/j.rse.2016.01.023) paper. To provide per-pixel uncertainty for NBAR, we use the SIAC calculated BRDF kernel weighting parameters and their uncertainty to calculate the uncertainty in the c factor ($\sigma_c$). There is a high correlation between the simulated surface reflectance at the nadir viewing geometry and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. It is important to use full uncertainty propagation to calculate the uncertainty in the c factor. The uncertainty in the c factor is calculated as follows:

The standard equation for uncertainty propagation for $y = \frac{A}{B}$ is given as follows:

$$
\begin{align}
\sigma_y = |y| \sqrt{(\frac{\sigma_A}{A})^2 + (\frac{\sigma_B}{B})^2 - 2 \frac{\sigma_{AB}}{AB}}
\end{align}
$$

where $\sigma_y$ is the uncertainty in $y$, $\sigma_A$ is the uncertainty in $A$, $\sigma_B$ is the uncertainty in $B$ and $\sigma_{AB}$ is the covariance between $A$ and $B$. The covariance between $A$ and $B$ can be calculated from the correlation coefficient ($r$) as follows:

$$
\begin{align}
\sigma_{AB} = r \sigma_A \sigma_B
\end{align}
$$

So, the standard equation for uncertainty propagation for $y = \frac{A}{B}$ is given as follows:

$$
\begin{align}
\sigma_y = |y| \sqrt{(\frac{\sigma_A}{A})^2 + (\frac{\sigma_B}{B})^2 - 2 \frac{r \sigma_A \sigma_B}{AB}}
\end{align}
$$

In the case of c factor, the $y$ is the c factor, the $A$ is the $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ and the $B$ is $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. The uncertainty in c factor ($\sigma_c$) is calculated as follows:

$$
\begin{align}
\sigma_c = |c|\sqrt{(\frac{\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}}{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})})^2 + (\frac{\sigma_{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})})^2 - 2 \frac{p \sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})} \sigma_{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar}) \hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}
\end{align}
$$

where $\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}$ is the uncertainty in the $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$, $\sigma_{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}$ is the uncertainty in $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. The Pearson correlation coefficient ($p$) is the correlation coefficient between each pair of $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. However, it is inhibitive to calculate them for all the pixels; the correlation coefficient is then calculated from $\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}$ and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$ over the whole image:

$$
\begin{align}
p = \frac{cov(\hat{R}^{BRDF}(\theta_v=0, \theta_s^{nbar}), \hat{R}^{BRDF}(\theta_v^{s2}, \theta_s^{s2}))}{\sigma_{\hat{R}^{BRDF}(\theta_v=0, \theta_s^{nbar})} \sigma_{\hat{R}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}
\end{align}
$$

where $\hat{R}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ is a vector of all $\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}$ and  $\hat{R}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$ is a vector of all $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$.

### Uncertainty in NBAR reflectance

Since we have now have the uncertainty in the c factor, we can calculate the uncertainty in the NBAR reflectance ($r_{nbar}$) as follows:

$$
\begin{align}
\sigma_{r_{nbar}} = \sqrt{\sigma_{r_{s2}}^2 c^2 + \sigma_c^2 r_{s2}^2}
\end{align}
$$

where $\sigma_{r_{s2}}$ is the uncertainty in the S2 reflectance ($r_{s2}$).


### MODIS BRDF uncertainty (appropriateness)

Apart from the uncertainty in the c factor, there is also uncertainty of whether the MODIS BRDF kernel is appropriate for S2 pixels within the MODIS 500mx500m pixel. This is represented by the difference between S2 reflectance ($r_{s2}$) and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. In this way, S2 pixel reflectance that are significantly different from  $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$  will have a higher uncertainty. This appropriateness uncertainty is calculated as follows:

$$
\begin{align}
\sigma_{app} = \sqrt{(\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2}) - r_{s2})^2}
\end{align}
$$

where the $\sigma_{app}$ is the appropriateness uncertainty.

### Combined uncertainty

Finally, the combined uncertainty in the NBAR reflectance ($\sigma_{r_{nbar}} $) is calculated as follows:

$$
\begin{align}
\sigma_{r_{nbar}}  = \sqrt{\sigma_{r_{s2}}^2 c^2 + \sigma_c^2 r_{s2}^2 + \sigma_{app}^2}
\end{align}
$$

