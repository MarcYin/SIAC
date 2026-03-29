#![allow(non_local_definitions)]
//! PSF Convolution
//!
//! DCT-based Gaussian PSF convolution for scale matching between
//! high-resolution satellite imagery (e.g., Sentinel-2 at 10m) and
//! coarse-resolution BRDF products (e.g., MODIS at 500m).

use ndarray::{Array2, ArrayView2};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;

/// PSF Convolver using DCT-based Gaussian convolution
///
/// Efficiently convolves high-resolution imagery with a Gaussian PSF
/// to match the spatial response of coarser resolution sensors.
#[pyclass]
pub struct PSFConvolver {
    /// PSF standard deviation in x direction (pixels)
    sigma_x: f64,
    /// PSF standard deviation in y direction (pixels)
    sigma_y: f64,
}

#[pymethods]
impl PSFConvolver {
    /// Create a new PSF convolver
    ///
    /// # Arguments
    /// * `sigma_x` - PSF standard deviation in x direction (pixels)
    /// * `sigma_y` - PSF standard deviation in y direction (pixels)
    #[new]
    #[pyo3(signature = (sigma_x=29.75, sigma_y=39.0))]
    fn new(sigma_x: f64, sigma_y: f64) -> Self {
        Self { sigma_x, sigma_y }
    }

    /// Apply PSF convolution to an image
    ///
    /// # Arguments
    /// * `image` - Input 2D array
    ///
    /// # Returns
    /// Convolved image array
    fn convolve<'py>(
        &self,
        py: Python<'py>,
        image: PyReadonlyArray2<f64>,
    ) -> PyResult<&'py PyArray2<f64>> {
        let image = image.as_array();
        let result = self.dct_convolve(image);
        Ok(result.into_pyarray(py))
    }
}

impl PSFConvolver {
    /// DCT-based Gaussian convolution.
    ///
    /// **NOTE**: This is a simplified spatial-domain implementation that applies
    /// the PSF via direct Gaussian-weighted averaging.  For production use,
    /// prefer the Python-side `scipy.ndimage.gaussian_filter` path in
    /// `kernel_model.KernelModelDeriver._convolve_2d` which handles NaN masking
    /// correctly.
    fn dct_convolve(&self, image: ArrayView2<f64>) -> Array2<f64> {
        let (height, width) = (image.shape()[0], image.shape()[1]);

        // Spatial-domain Gaussian convolution (separable, truncated at 4σ).
        // We apply the filter in two passes: first along rows (x), then cols (y).
        let mut tmp = Array2::zeros((height, width));
        let radius_x = (4.0 * self.sigma_x).ceil() as usize;
        let radius_y = (4.0 * self.sigma_y).ceil() as usize;

        // Pre-compute 1-D Gaussian weights for x and y.
        let weights_x: Vec<f64> = (0..=radius_x)
            .map(|d| (-0.5 * (d as f64 / self.sigma_x).powi(2)).exp())
            .collect();
        let weights_y: Vec<f64> = (0..=radius_y)
            .map(|d| (-0.5 * (d as f64 / self.sigma_y).powi(2)).exp())
            .collect();

        // Pass 1: convolve along x (columns) for each row.
        for i in 0..height {
            for j in 0..width {
                let mut sum = 0.0_f64;
                let mut wsum = 0.0_f64;
                let j_lo = if j >= radius_x { j - radius_x } else { 0 };
                let j_hi = (j + radius_x).min(width - 1);
                for jj in j_lo..=j_hi {
                    let d = if jj >= j { jj - j } else { j - jj };
                    let w = weights_x[d];
                    let v = image[[i, jj]];
                    if v.is_finite() {
                        sum += w * v;
                        wsum += w;
                    }
                }
                tmp[[i, j]] = if wsum > 0.0 { sum / wsum } else { f64::NAN };
            }
        }

        // Pass 2: convolve along y (rows) for each column.
        let mut result = Array2::zeros((height, width));
        for j in 0..width {
            for i in 0..height {
                let mut sum = 0.0_f64;
                let mut wsum = 0.0_f64;
                let i_lo = if i >= radius_y { i - radius_y } else { 0 };
                let i_hi = (i + radius_y).min(height - 1);
                for ii in i_lo..=i_hi {
                    let d = if ii >= i { ii - i } else { i - ii };
                    let w = weights_y[d];
                    let v = tmp[[ii, j]];
                    if v.is_finite() {
                        sum += w * v;
                        wsum += w;
                    }
                }
                result[[i, j]] = if wsum > 0.0 { sum / wsum } else { f64::NAN };
            }
        }

        result
    }
}
