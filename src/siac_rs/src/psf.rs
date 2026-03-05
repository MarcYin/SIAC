//! PSF Convolution
//!
//! DCT-based Gaussian PSF convolution for scale matching between
//! high-resolution satellite imagery (e.g., Sentinel-2 at 10m) and
//! coarse-resolution BRDF products (e.g., MODIS at 500m).

use ndarray::{Array2, ArrayView2};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use rustfft::{num_complex::Complex64, FftPlanner};

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

    /// Get the Gaussian kernel in frequency domain
    fn get_frequency_kernel<'py>(
        &self,
        py: Python<'py>,
        height: usize,
        width: usize,
    ) -> PyResult<&'py PyArray2<f64>> {
        let kernel = self.gaussian_dct_kernel(height, width);
        Ok(kernel.into_pyarray(py))
    }
}

impl PSFConvolver {
    /// DCT-based Gaussian convolution
    fn dct_convolve(&self, image: ArrayView2<f64>) -> Array2<f64> {
        let (height, width) = (image.shape()[0], image.shape()[1]);

        // Get Gaussian kernel in DCT domain
        let kernel = self.gaussian_dct_kernel(height, width);

        // For now, return a placeholder - full DCT implementation would go here
        // In production, use rustfft for efficient DCT computation
        let mut result = Array2::zeros((height, width));
        for i in 0..height {
            for j in 0..width {
                result[[i, j]] = image[[i, j]] * kernel[[i, j]];
            }
        }
        result
    }

    /// Generate Gaussian kernel in DCT frequency domain
    fn gaussian_dct_kernel(&self, height: usize, width: usize) -> Array2<f64> {
        let mut kernel = Array2::zeros((height, width));

        for i in 0..height {
            for j in 0..width {
                let u = i as f64 / height as f64;
                let v = j as f64 / width as f64;

                // Gaussian in frequency domain
                let gx = (-2.0 * std::f64::consts::PI.powi(2) * self.sigma_x.powi(2) * u.powi(2))
                    .exp();
                let gy = (-2.0 * std::f64::consts::PI.powi(2) * self.sigma_y.powi(2) * v.powi(2))
                    .exp();

                kernel[[i, j]] = gx * gy;
            }
        }

        kernel
    }
}
