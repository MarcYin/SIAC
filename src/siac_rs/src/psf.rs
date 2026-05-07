#![allow(non_local_definitions)]
//! PSF Convolution
//!
//! Direct spatial-domain Gaussian-weighted convolution for scale matching
//! between high-resolution satellite imagery (e.g., Sentinel-2 at 10m) and
//! coarse-resolution BRDF products (e.g., MODIS at 500m).
//!
//! Despite the legacy `dct_convolve` method name (kept for backwards
//! compatibility with existing call sites), the actual implementation is a
//! separable spatial-domain Gaussian filter truncated at 4σ with NaN-aware
//! weight normalisation — there is no DCT step. See the function comment on
//! `dct_convolve` for details.

use ndarray::{Array2, ArrayView2};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;

/// PSF Convolver using direct spatial-domain Gaussian convolution.
///
/// Convolves high-resolution imagery with a Gaussian PSF to match the
/// spatial response of coarser resolution sensors. Uses a separable
/// (two-pass) implementation truncated at 4σ.
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
    /// Direct spatial-domain Gaussian-weighted convolution (legacy name
    /// `dct_convolve` kept for backwards compatibility).
    ///
    /// Applies the PSF via a separable two-pass Gaussian filter truncated at
    /// 4σ. Each output pixel is a NaN-aware weighted mean of the input pixels
    /// inside the kernel footprint: non-finite inputs are dropped from both
    /// numerator and denominator so a single NaN does not contaminate the
    /// neighbourhood, and an output is NaN only when *every* contributing
    /// input is non-finite.
    ///
    /// The legacy `dct_convolve` name reflects an earlier DCT-based prototype;
    /// the current implementation is purely spatial-domain. For the production
    /// Python path, prefer `kernel_model.KernelModelDeriver._convolve_2d`
    /// which wraps `scipy.ndimage.gaussian_filter` with the same NaN handling.
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
                let j_lo = j.saturating_sub(radius_x);
                let j_hi = (j + radius_x).min(width - 1);
                for jj in j_lo..=j_hi {
                    let d = jj.abs_diff(j);
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
                let i_lo = i.saturating_sub(radius_y);
                let i_hi = (i + radius_y).min(height - 1);
                for ii in i_lo..=i_hi {
                    let d = ii.abs_diff(i);
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array2;

    #[test]
    fn convolution_preserves_constant_field() {
        // A separable Gaussian filter with weight-normalised edge handling
        // must reproduce a uniform input exactly (sum of weights cancels).
        let conv = PSFConvolver::new(1.5, 2.0);
        let image = Array2::from_elem((12, 10), 0.42_f64);
        let out = conv.dct_convolve(image.view());
        for v in out.iter() {
            assert_relative_eq!(*v, 0.42, epsilon = 1e-12);
        }
    }

    #[test]
    fn convolution_preserves_shape() {
        let conv = PSFConvolver::new(1.0, 1.5);
        let image = Array2::<f64>::zeros((7, 5));
        let out = conv.dct_convolve(image.view());
        assert_eq!(out.shape(), image.shape());
    }

    #[test]
    fn nan_pixels_are_ignored_in_weighted_mean() {
        // A single NaN in an otherwise-uniform field must not contaminate
        // output values: finite-only weighted averaging gives back the
        // underlying constant at every output pixel.
        let conv = PSFConvolver::new(1.0, 1.0);
        let mut image = Array2::from_elem((6, 6), 1.0_f64);
        image[[3, 3]] = f64::NAN;
        let out = conv.dct_convolve(image.view());
        for v in out.iter() {
            assert!(v.is_finite(), "output must be finite despite input NaN");
            assert_relative_eq!(*v, 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn all_nan_input_produces_all_nan_output() {
        let conv = PSFConvolver::new(1.0, 1.0);
        let image = Array2::from_elem((4, 4), f64::NAN);
        let out = conv.dct_convolve(image.view());
        assert!(out.iter().all(|v| v.is_nan()));
    }

    #[test]
    fn impulse_broadens_with_larger_sigma() {
        // A point source should spread over a wider area as sigma increases.
        // Compare the central-row sum of two convolutions: the wider kernel
        // will push more mass away from the impulse column.
        let mut image = Array2::<f64>::zeros((11, 11));
        image[[5, 5]] = 1.0;
        let narrow = PSFConvolver::new(0.8, 0.8).dct_convolve(image.view());
        let wide = PSFConvolver::new(2.5, 2.5).dct_convolve(image.view());
        assert!(
            narrow[[5, 5]] > wide[[5, 5]],
            "narrow PSF should leave more mass at the impulse center"
        );
    }
}
