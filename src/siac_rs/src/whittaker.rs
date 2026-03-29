//! Whittaker smoothing utilities for temporal BRDF priors.

use nalgebra::{DMatrix, DVector};
use ndarray::Array4;
use numpy::{IntoPyArray, PyArray4, PyReadonlyArray4};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

fn second_difference_penalty(n_time: usize, lambda: f64) -> DMatrix<f64> {
    if n_time < 3 {
        return DMatrix::<f64>::zeros(n_time, n_time);
    }

    let mut d2 = DMatrix::<f64>::zeros(n_time - 2, n_time);
    for row in 0..(n_time - 2) {
        d2[(row, row)] = 1.0;
        d2[(row, row + 1)] = -2.0;
        d2[(row, row + 2)] = 1.0;
    }
    (d2.transpose() * d2) * lambda
}

/// Apply Whittaker smoothing along the time axis of a `(time, band, y, x)` cube.
#[pyfunction]
pub fn whittaker_smooth_cube<'py>(
    py: Python<'py>,
    values: PyReadonlyArray4<f32>,
    weights: PyReadonlyArray4<f32>,
    lambda_: f64,
) -> PyResult<&'py PyArray4<f32>> {
    if lambda_ < 0.0 {
        return Err(PyValueError::new_err("lambda_ must be >= 0"));
    }

    let values = values.as_array();
    let weights = weights.as_array();
    if values.shape() != weights.shape() {
        return Err(PyValueError::new_err(
            "values and weights must share the same (time, band, y, x) shape",
        ));
    }

    let shape = values.shape();
    let (n_time, n_band, ny, nx) = (shape[0], shape[1], shape[2], shape[3]);
    let penalty = second_difference_penalty(n_time, lambda_);
    // Store penalty as a flat Vec for cheaper cloning per-pixel.
    let penalty_flat: Vec<f64> = penalty.as_slice().to_vec();
    let series_count = n_band * ny * nx;

    let smoothed_series: Vec<Vec<f32>> = (0..series_count)
        .into_par_iter()
        .map(|flat_idx| {
            let band = flat_idx / (ny * nx);
            let rem = flat_idx % (ny * nx);
            let y = rem / nx;
            let x = rem % nx;

            let mut rhs = DVector::<f64>::zeros(n_time);
            // Reconstruct from flat vec – avoids nalgebra clone overhead.
            let mut system = DMatrix::from_column_slice(n_time, n_time, &penalty_flat);
            let mut has_valid = false;

            for t in 0..n_time {
                let value = values[[t, band, y, x]] as f64;
                let weight = weights[[t, band, y, x]] as f64;
                if value.is_finite() && weight.is_finite() && weight > 0.0 {
                    system[(t, t)] += weight;
                    rhs[t] = weight * value;
                    has_valid = true;
                }
            }

            if !has_valid {
                return vec![f32::NAN; n_time];
            }

            // System is symmetric positive-definite → Cholesky is ~2× faster than LU.
            if let Some(chol) = system.clone().cholesky() {
                let solution = chol.solve(&rhs);
                solution.iter().map(|value| *value as f32).collect()
            } else {
                // Fallback to LU if Cholesky fails (e.g. near-zero lambda).
                let lu = system.lu();
                if let Some(solution) = lu.solve(&rhs) {
                    solution.iter().map(|value| *value as f32).collect()
                } else {
                    vec![f32::NAN; n_time]
                }
            }
        })
        .collect();

    let mut out = Array4::<f32>::from_elem((n_time, n_band, ny, nx), f32::NAN);
    for flat_idx in 0..series_count {
        let band = flat_idx / (ny * nx);
        let rem = flat_idx % (ny * nx);
        let y = rem / nx;
        let x = rem % nx;
        let series = &smoothed_series[flat_idx];
        for t in 0..n_time {
            out[[t, band, y, x]] = series[t];
        }
    }

    Ok(out.into_pyarray(py))
}
