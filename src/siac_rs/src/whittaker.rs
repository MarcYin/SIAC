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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn penalty_zero_for_short_series() {
        // With fewer than 3 time steps the second-difference operator is
        // undefined; penalty matrix must collapse to zero so the solver
        // reduces to weighted least squares.
        let p = second_difference_penalty(2, 10.0);
        assert_eq!(p.shape(), (2, 2));
        assert!(p.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn penalty_structure_is_tridiagonal_like() {
        // D2^T D2 produces a pentadiagonal symmetric matrix with known
        // interior row pattern [1, -4, 6, -4, 1] scaled by lambda.
        let p = second_difference_penalty(5, 1.0);
        assert_relative_eq!(p[(2, 0)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(p[(2, 1)], -4.0, epsilon = 1e-12);
        assert_relative_eq!(p[(2, 2)], 6.0, epsilon = 1e-12);
        assert_relative_eq!(p[(2, 3)], -4.0, epsilon = 1e-12);
        assert_relative_eq!(p[(2, 4)], 1.0, epsilon = 1e-12);
        // Symmetric.
        for i in 0..5 {
            for j in 0..5 {
                assert_relative_eq!(p[(i, j)], p[(j, i)], epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn penalty_scales_linearly_with_lambda() {
        let a = second_difference_penalty(6, 1.0);
        let b = second_difference_penalty(6, 7.5);
        for i in 0..6 {
            for j in 0..6 {
                assert_relative_eq!(b[(i, j)], 7.5 * a[(i, j)], epsilon = 1e-12);
            }
        }
    }
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

    // Heavy compute uses only Rust-owned buffers (`values`/`weights` are
    // `PyReadonlyArray` views over the underlying NumPy buffers, safe to
    // access without the GIL because we hold the readonly guard). Release
    // the GIL while the rayon pool runs so other Python threads can progress.
    let smoothed_series: Vec<Vec<f32>> = py.allow_threads(|| {
        (0..series_count)
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
            .collect()
    });

    let mut out = Array4::<f32>::from_elem((n_time, n_band, ny, nx), f32::NAN);
    for (flat_idx, series) in smoothed_series.iter().enumerate().take(series_count) {
        let band = flat_idx / (ny * nx);
        let rem = flat_idx % (ny * nx);
        let y = rem / nx;
        let x = rem % nx;
        for t in 0..n_time {
            out[[t, band, y, x]] = series[t];
        }
    }

    Ok(out.into_pyarray(py))
}
