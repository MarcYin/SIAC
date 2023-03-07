# SIAC nadir BRDF-adjusted reflectance (NBAR) and uncertainty calculation

This is a description of the SIAC nadir BRDF-adjusted reflectance (NBAR) product. It closely follows [Roy et al. 2016](https://doi.org/10.1016/j.rse.2016.01.023) paper. The correction is down by normalising the viewing geometry and the solar geometry to the nadir geometry with a given solar zenith angle (SZA). For generating NBAR for time series of images, the SZA can be constant or variable (changing the SZA over the year). 

## Considerations on SZA

As there is no BRDF information available from Sentinel-2 (S2), the MODIS/VIIRS BRDF information is used. This would cause issues when the surface property is not homogeneous, as within the MODIS 500meterx500meter pixel the BRDF shape could be vastly different at the scale of S2 10metersx10meters pixel. Therefore, it is recommended to not extrapolate the NBAR SZA too far from the initial S2 SZA. If the application requires images from an individual date and can handle the variation in SZA, you may choose to normalise the SZA to the original S2 reported SZA or the SZA calculated from the image latitude and longitude. However, if the application requires mosaic from time series of images, it is recommended to use the mean SZA calculated over time series of images.



## BRDF kernel
Before introducing the NBAR calculation, we first introduce the BRDF kernel weighting parameters and their uncertainty calculation. The BRDF kernel weighting parameters are used to calculate the BRDF reflectance. To simplify the notation in the uncertainty calculation, we dropped the wavelength subscript in the following calculation. The simulated surface reflectance at the viewing geometry ($\theta_v$) and the solar geometry ($\theta_s$) is calculated as follows with the BRDF kernel weighting parameters ($f_{iso}$, $f_{vol}$ and $f_{geo}$) and the BRDF kernels ($K_{vol}$ and $K_{geo}$):

$$
\begin{align}
\hat{r}^{BRDF}(\theta_v, \theta_s) = f_{iso} + f_{vol}K_{vol}(\theta_v, \theta_s) + f_{geo}K_{geo}(\theta_v, \theta_s)
\end{align}
$$


The uncertainty in the simulate surface reflectance ($\hat{r}^{BRDF}_{unc}$) is calculated as follows:

$$
\begin{align}
\hat{r}^{BRDF}_{unc}(\theta_v, \theta_s) = \sqrt{(\frac{\partial \hat{r}^{BRDF}}{\partial f_{iso}} \times f_{iso_{unc}})^2 + (\frac{\partial \hat{r}^{BRDF}}{\partial f_{vol}} \times f_{vol_{unc}})^2 + (\frac{\partial \hat{r}^{BRDF}}{\partial f_{geo}} \times f_{geo_{unc}})^2}
\end{align}
$$

where $\frac{\partial \hat{r}^{BRDF}}{\partial f_{iso}} = 1$, $\frac{\partial \hat{r}^{BRDF}}{\partial f_{vol}} = K_{vol}(\theta_v, \theta_s)$ and $\frac{\partial \hat{r}^{BRDF}}{\partial f_{geo}} = K_{geo}(\theta_v, \theta_s)$ are the partial derivatives of BRDF reflectance with respect to the weighting parameters of the BRDF kernel ($f_{iso}$, $f_{vol}$ and $f_{geo}$). The uncertainty in the BRDF kernel weighting parameters ($f_{iso_{unc}}$, $f_{vol_{unc}}$ and $f_{geo_{unc}}$) is calculated using the [SIAC](https://gmd.copernicus.org/articles/15/7933/2022/) smoothed BRDF kernel weighting parameters and their uncertainty. 

So, the uncertainty in the BRDF reflectance ($\hat{r}^{BRDF}_{unc}$) is calculated as follows:

$$
\begin{align}
\hat{r}^{BRDF}_{unc}(\theta_v, \theta_s) = \sqrt{(f_{iso_{unc}})^2 + (f_{vol_{unc}} \times K_{vol}(\theta_v, \theta_s))^2 + (f_{geo_{unc}} \times K_{geo}(\theta_v, \theta_s))^2}
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

## NBAR uncertainty calculation

### Uncerainty in c factor
The uncertainty calculation for NBAR is not given in the [Roy et al. 2016](https://doi.org/10.1016/j.rse.2016.01.023) paper. To provide per-pixel uncertainty for NBAR, we use the SIAC calculated BRDF kernel weighting parameters and their uncertainty to calculate the uncertainty in the c factor ($c_{unc}$). There is a high correlation between the simulated surface reflectance at the nadir viewing geometry and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. It is important to use full uncertainty propagation to calculate the uncertainty in the c factor. The uncertainty in the c factor is calculated as follows:

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
y_{unc} = |y| \sqrt{(\frac{\sigma_A}{A})^2 + (\frac{\sigma_B}{B})^2 - 2 \frac{r \sigma_A \sigma_B}{AB}}
\end{align}
$$

In the case of c factor, the $y$ is the c factor, the $A$ is the $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ and the $B$ is $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. The uncertainty in c factor ($c_{unc}$) is calculated as follows:

$$
\begin{align}
c_{unc} = |c|\sqrt{(\frac{\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}}{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})})^2 + (\frac{\sigma_{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})})^2 - 2 \frac{r \sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})} \sigma_{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar}) \hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}
\end{align}
$$

where $\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}$ is the uncertainty in the $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ and the given SZA ($\theta_s^{nbar}$), $\sigma_{\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}$ is the uncertainty in $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. The correlation coefficient ($r$) is the correlation coefficient between each pair of $\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$. However, it is inhibitive to calculate them for all the pixels; the correlation coefficient is calculated from $\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}$ and $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$ over the whole image:

$$
\begin{align}
r = \frac{cov(\hat{R}^{BRDF}(\theta_v=0, \theta_s^{nbar}), \hat{R}^{BRDF}(\theta_v^{s2}, \theta_s^{s2}))}{\sigma_{\hat{R}^{BRDF}(\theta_v=0, \theta_s^{nbar})} \sigma_{\hat{R}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})}}
\end{align}
$$

where $\hat{R}^{BRDF}(\theta_v=0, \theta_s^{nbar})$ is a vector of all $\sigma_{\hat{r}^{BRDF}(\theta_v=0, \theta_s^{nbar})}$ and  $\hat{R}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$ is a vector of all $\hat{r}^{BRDF}(\theta_v^{s2}, \theta_s^{s2})$.

### Uncertainty in NBAR reflectance

Since we have now have the uncertainty in the c factor, we can calculate the uncertainty in the NBAR reflectance ($r_{nbar}$) as follows:

$$
\begin{align}
r_{nbar, unc} = \sqrt{\sigma_{r_{s2}}^2 c^2 + c_{unc}^2 r_{s2}^2}
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

Finally, the combined uncertainty in the NBAR reflectance ($r_{nbar}$) is calculated as follows:

$$
\begin{align}
r_{nbar, unc} = \sqrt{\sigma_{r_{s2}}^2 c^2 + c_{unc}^2 r_{s2}^2 + \sigma_{app}^2}
\end{align}
$$


