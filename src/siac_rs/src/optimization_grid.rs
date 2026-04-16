//! Pure grid-remap, interpolation, and Laplacian helpers extracted from
//! `optimization.rs`.  These three functions have no dependency on the
//! multigrid cost/refinement machinery and live here so the main file can
//! focus on solver logic.

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;

/// Remap high-resolution data to coarse grid by averaging.
///
/// Each high-res pixel is mapped center-based to the coarse cell it falls
/// into (to avoid aliasing); the coarse value is the arithmetic mean of the
/// finite high-res pixels in its cell.  NaN pixels are ignored; coarse cells
/// receiving no finite pixels are left at zero.
///
/// # Memory
/// Returns a new NumPy array of `coarse_rows * coarse_cols * 8` bytes (f64).
#[pyfunction]
pub fn remap_to_coarse_grid<'py>(
    py: Python<'py>,
    data: PyReadonlyArray2<f64>,
    coarse_rows: usize,
    coarse_cols: usize,
) -> PyResult<&'py PyArray2<f64>> {
    let data = data.as_array();
    let (high_rows, high_cols) = (data.shape()[0], data.shape()[1]);

    let mut coarse = Array2::<f64>::zeros((coarse_rows, coarse_cols));
    let mut counts = Array2::<f64>::zeros((coarse_rows, coarse_cols));

    for i in 0..high_rows {
        for j in 0..high_cols {
            let ci = ((2 * i + 1) * coarse_rows / (2 * high_rows)).min(coarse_rows - 1);
            let cj = ((2 * j + 1) * coarse_cols / (2 * high_cols)).min(coarse_cols - 1);

            let val = data[[i, j]];
            if val.is_finite() {
                coarse[[ci, cj]] += val;
                counts[[ci, cj]] += 1.0;
            }
        }
    }

    for i in 0..coarse_rows {
        for j in 0..coarse_cols {
            if counts[[i, j]] > 0.0 {
                coarse[[i, j]] /= counts[[i, j]];
            }
        }
    }

    Ok(coarse.into_pyarray(py))
}

/// Interpolate coarse grid to fine grid using bilinear interpolation.
///
/// # Memory
/// Returns a new NumPy array of `fine_rows * fine_cols * 8` bytes (f64).
#[pyfunction]
pub fn interpolate_to_fine_grid<'py>(
    py: Python<'py>,
    coarse: PyReadonlyArray2<f64>,
    fine_rows: usize,
    fine_cols: usize,
) -> PyResult<&'py PyArray2<f64>> {
    let coarse = coarse.as_array();
    let (coarse_rows, coarse_cols) = (coarse.shape()[0], coarse.shape()[1]);

    let mut fine = Array2::<f64>::zeros((fine_rows, fine_cols));

    for i in 0..fine_rows {
        for j in 0..fine_cols {
            let ci = (i as f64 + 0.5) * coarse_rows as f64 / fine_rows as f64 - 0.5;
            let cj = (j as f64 + 0.5) * coarse_cols as f64 / fine_cols as f64 - 0.5;

            let ci0 = (ci.floor() as usize).min(coarse_rows - 1);
            let ci1 = (ci0 + 1).min(coarse_rows - 1);
            let cj0 = (cj.floor() as usize).min(coarse_cols - 1);
            let cj1 = (cj0 + 1).min(coarse_cols - 1);

            let wi = (ci - ci0 as f64).clamp(0.0, 1.0);
            let wj = (cj - cj0 as f64).clamp(0.0, 1.0);

            fine[[i, j]] = (1.0 - wi) * (1.0 - wj) * coarse[[ci0, cj0]]
                + (1.0 - wi) * wj * coarse[[ci0, cj1]]
                + wi * (1.0 - wj) * coarse[[ci1, cj0]]
                + wi * wj * coarse[[ci1, cj1]];
        }
    }

    Ok(fine.into_pyarray(py))
}

/// Compute Laplacian difference matrix products for smoothness regularization.
///
/// Computes `D^T @ D @ x` for both x and y directions where D is the
/// first-order difference matrix.  Used in the smoothness cost function:
/// `J_smooth = 0.5 * x^T @ (Dx^T @ Dx + Dy^T @ Dy) @ x / sigma^2`.
#[pyfunction]
pub fn apply_laplacian<'py>(
    py: Python<'py>,
    data: PyReadonlyArray2<f64>,
) -> PyResult<&'py PyArray2<f64>> {
    let data = data.as_array();
    let (rows, cols) = (data.shape()[0], data.shape()[1]);

    let mut result = Array2::<f64>::zeros((rows, cols));

    // D^T @ D along rows (first-order difference in x).
    for i in 0..rows {
        for j in 0..cols {
            let central_coef = if i == 0 || i == rows - 1 { 1.0 } else { 2.0 };
            let mut val = central_coef * data[[i, j]];
            if i > 0 {
                val -= data[[i - 1, j]];
            }
            if i < rows - 1 {
                val -= data[[i + 1, j]];
            }
            result[[i, j]] += val;
        }
    }

    // D^T @ D along columns (first-order difference in y).
    for i in 0..rows {
        for j in 0..cols {
            let central_coef = if j == 0 || j == cols - 1 { 1.0 } else { 2.0 };
            let mut val = central_coef * data[[i, j]];
            if j > 0 {
                val -= data[[i, j - 1]];
            }
            if j < cols - 1 {
                val -= data[[i, j + 1]];
            }
            result[[i, j]] += val;
        }
    }

    Ok(result.into_pyarray(py))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use ndarray::{Array2, ArrayView2};

    // Local inline copy of the algorithm bodies so the tests don't require a
    // Python runtime; they validate the numeric behaviour of the kernels.

    fn remap_local(data: ArrayView2<f64>, coarse_rows: usize, coarse_cols: usize) -> Array2<f64> {
        let (high_rows, high_cols) = (data.shape()[0], data.shape()[1]);
        let mut coarse = Array2::<f64>::zeros((coarse_rows, coarse_cols));
        let mut counts = Array2::<f64>::zeros((coarse_rows, coarse_cols));
        for i in 0..high_rows {
            for j in 0..high_cols {
                let ci = ((2 * i + 1) * coarse_rows / (2 * high_rows)).min(coarse_rows - 1);
                let cj = ((2 * j + 1) * coarse_cols / (2 * high_cols)).min(coarse_cols - 1);
                let val = data[[i, j]];
                if val.is_finite() {
                    coarse[[ci, cj]] += val;
                    counts[[ci, cj]] += 1.0;
                }
            }
        }
        for i in 0..coarse_rows {
            for j in 0..coarse_cols {
                if counts[[i, j]] > 0.0 {
                    coarse[[i, j]] /= counts[[i, j]];
                }
            }
        }
        coarse
    }

    #[test]
    fn remap_identity_when_shapes_match() {
        let data = Array2::from_shape_fn((4, 4), |(i, j)| (i * 4 + j) as f64);
        let out = remap_local(data.view(), 4, 4);
        for i in 0..4 {
            for j in 0..4 {
                assert_relative_eq!(out[[i, j]], data[[i, j]], epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn remap_ignores_nan_pixels() {
        let mut data = Array2::<f64>::ones((4, 4));
        data[[0, 0]] = f64::NAN;
        let out = remap_local(data.view(), 2, 2);
        // Every coarse cell must be 1.0: either all ones, or three ones and
        // one NaN (which is ignored).
        for v in out.iter() {
            assert_relative_eq!(*v, 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn remap_empty_cell_defaults_to_zero() {
        // All NaNs -> coarse cells have no finite inputs, so remain zero.
        let data = Array2::from_elem((4, 4), f64::NAN);
        let out = remap_local(data.view(), 2, 2);
        for v in out.iter() {
            assert_eq!(*v, 0.0);
        }
    }
}
