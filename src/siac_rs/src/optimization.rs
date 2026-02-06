//! Multi-grid Optimization Utilities
//!
//! Helper functions for the multi-grid L-BFGS-B solver used in aerosol retrieval.

use ndarray::{Array2, ArrayView2};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use rayon::prelude::*;

/// Remap high-resolution data to coarse grid by averaging
///
/// # Arguments
/// * `data` - High-resolution 2D array
/// * `coarse_shape` - Target (rows, cols) for coarse grid
///
/// # Returns
/// Coarse-resolution array with averaged values
#[pyfunction]
fn remap_to_coarse_grid<'py>(
    py: Python<'py>,
    data: PyReadonlyArray2<f64>,
    coarse_rows: usize,
    coarse_cols: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let data = data.as_array();
    let (high_rows, high_cols) = (data.shape()[0], data.shape()[1]);

    let mut coarse = Array2::<f64>::zeros((coarse_rows, coarse_cols));
    let mut counts = Array2::<f64>::zeros((coarse_rows, coarse_cols));

    // Map each high-res pixel to coarse grid cell
    for i in 0..high_rows {
        for j in 0..high_cols {
            let ci = (i * coarse_rows) / high_rows;
            let cj = (j * coarse_cols) / high_cols;

            let val = data[[i, j]];
            if val.is_finite() {
                coarse[[ci, cj]] += val;
                counts[[ci, cj]] += 1.0;
            }
        }
    }

    // Average
    for i in 0..coarse_rows {
        for j in 0..coarse_cols {
            if counts[[i, j]] > 0.0 {
                coarse[[i, j]] /= counts[[i, j]];
            }
        }
    }

    Ok(coarse.into_pyarray_bound(py))
}

/// Interpolate coarse grid to fine grid using bilinear interpolation
///
/// # Arguments
/// * `coarse` - Coarse-resolution 2D array
/// * `fine_shape` - Target (rows, cols) for fine grid
///
/// # Returns
/// Fine-resolution array with interpolated values
#[pyfunction]
fn interpolate_to_fine_grid<'py>(
    py: Python<'py>,
    coarse: PyReadonlyArray2<f64>,
    fine_rows: usize,
    fine_cols: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let coarse = coarse.as_array();
    let (coarse_rows, coarse_cols) = (coarse.shape()[0], coarse.shape()[1]);

    let mut fine = Array2::<f64>::zeros((fine_rows, fine_cols));

    // Bilinear interpolation
    for i in 0..fine_rows {
        for j in 0..fine_cols {
            // Map to coarse coordinates (continuous)
            let ci = (i as f64 + 0.5) * coarse_rows as f64 / fine_rows as f64 - 0.5;
            let cj = (j as f64 + 0.5) * coarse_cols as f64 / fine_cols as f64 - 0.5;

            // Bilinear interpolation
            let ci0 = (ci.floor() as usize).min(coarse_rows - 1);
            let ci1 = (ci0 + 1).min(coarse_rows - 1);
            let cj0 = (cj.floor() as usize).min(coarse_cols - 1);
            let cj1 = (cj0 + 1).min(coarse_cols - 1);

            let wi = ci - ci0 as f64;
            let wj = cj - cj0 as f64;

            let wi = wi.clamp(0.0, 1.0);
            let wj = wj.clamp(0.0, 1.0);

            fine[[i, j]] = (1.0 - wi) * (1.0 - wj) * coarse[[ci0, cj0]]
                + (1.0 - wi) * wj * coarse[[ci0, cj1]]
                + wi * (1.0 - wj) * coarse[[ci1, cj0]]
                + wi * wj * coarse[[ci1, cj1]];
        }
    }

    Ok(fine.into_pyarray_bound(py))
}

/// Compute Laplacian difference matrix products for smoothness regularization
///
/// Computes D^T @ D @ x for both x and y directions where D is the first-order
/// difference matrix. This is used in the smoothness cost function:
/// J_smooth = 0.5 * x^T @ (Dx^T @ Dx + Dy^T @ Dy) @ x / sigma^2
#[pyfunction]
fn apply_laplacian<'py>(
    py: Python<'py>,
    data: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let data = data.as_array();
    let (rows, cols) = (data.shape()[0], data.shape()[1]);

    let mut result = Array2::<f64>::zeros((rows, cols));

    // Apply D^T @ D in x direction (row differences)
    for i in 0..rows {
        for j in 0..cols {
            let mut val = 0.0;

            // Central coefficient: 2 for interior, 1 for edges
            let central_coef = if i == 0 || i == rows - 1 { 1.0 } else { 2.0 };
            val += central_coef * data[[i, j]];

            // Left neighbor
            if i > 0 {
                val -= data[[i - 1, j]];
            }

            // Right neighbor
            if i < rows - 1 {
                val -= data[[i + 1, j]];
            }

            result[[i, j]] += val;
        }
    }

    // Apply D^T @ D in y direction (column differences)
    for i in 0..rows {
        for j in 0..cols {
            let mut val = 0.0;

            // Central coefficient
            let central_coef = if j == 0 || j == cols - 1 { 1.0 } else { 2.0 };
            val += central_coef * data[[i, j]];

            // Top neighbor
            if j > 0 {
                val -= data[[i, j - 1]];
            }

            // Bottom neighbor
            if j < cols - 1 {
                val -= data[[i, j + 1]];
            }

            result[[i, j]] += val;
        }
    }

    Ok(result.into_pyarray_bound(py))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_remap_identity() {
        // Same size should preserve values
        let data = Array2::from_shape_vec((4, 4), (0..16).map(|x| x as f64).collect()).unwrap();

        // This would require Python context, so we test the logic manually
        let (high_rows, high_cols) = (4, 4);
        let (coarse_rows, coarse_cols) = (4, 4);

        let mut coarse = Array2::<f64>::zeros((coarse_rows, coarse_cols));
        let mut counts = Array2::<f64>::zeros((coarse_rows, coarse_cols));

        for i in 0..high_rows {
            for j in 0..high_cols {
                let ci = (i * coarse_rows) / high_rows;
                let cj = (j * coarse_cols) / high_cols;
                coarse[[ci, cj]] += data[[i, j]];
                counts[[ci, cj]] += 1.0;
            }
        }

        for i in 0..coarse_rows {
            for j in 0..coarse_cols {
                if counts[[i, j]] > 0.0 {
                    coarse[[i, j]] /= counts[[i, j]];
                }
            }
        }

        // Check values match original
        for i in 0..4 {
            for j in 0..4 {
                assert!((coarse[[i, j]] - data[[i, j]]).abs() < 1e-10);
            }
        }
    }
}
