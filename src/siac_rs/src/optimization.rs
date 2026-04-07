//! Multi-grid Optimization Utilities
//!
//! Helper functions for the multi-grid L-BFGS-B solver used in aerosol retrieval.

use ndarray::{s, Array2, Array4, ArrayView1, ArrayView2, ArrayView3, ArrayView4, Axis};
use numpy::{
    IntoPyArray, PyArray2, PyArray3, PyArray4, PyReadonlyArray1, PyReadonlyArray2,
    PyReadonlyArray3, PyReadonlyArray4,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FixedParameter {
    None,
    Aot,
    Tcwv,
}

fn parse_fixed_parameter(value: &str) -> PyResult<FixedParameter> {
    match value.trim().to_ascii_lowercase().as_str() {
        "" | "none" | "false" | "no" => Ok(FixedParameter::None),
        "aot" => Ok(FixedParameter::Aot),
        "tcwv" => Ok(FixedParameter::Tcwv),
        _ => Err(PyValueError::new_err(
            "fixed_parameter must be one of 'none', 'aot', or 'tcwv'",
        )),
    }
}

fn validate_fixed_axes(
    fixed_parameter: FixedParameter,
    n_aot: usize,
    n_tcwv: usize,
) -> PyResult<()> {
    match fixed_parameter {
        FixedParameter::Aot if n_aot != 1 => Err(PyValueError::new_err(
            "fixed_parameter='aot' requires aot_axis length 1",
        )),
        FixedParameter::Tcwv if n_tcwv != 1 => Err(PyValueError::new_err(
            "fixed_parameter='tcwv' requires tcwv_axis length 1",
        )),
        _ => Ok(()),
    }
}

/// Remap high-resolution data to coarse grid by averaging
///
/// # Arguments
/// * `data` - High-resolution 2D array
/// * `coarse_shape` - Target (rows, cols) for coarse grid
///
/// # Returns
/// Coarse-resolution array with averaged values
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

    // Map each high-res pixel to coarse grid cell (center-based to avoid aliasing)
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

    // Average
    for i in 0..coarse_rows {
        for j in 0..coarse_cols {
            if counts[[i, j]] > 0.0 {
                coarse[[i, j]] /= counts[[i, j]];
            }
        }
    }

    Ok(coarse.into_pyarray(py))
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
pub fn interpolate_to_fine_grid<'py>(
    py: Python<'py>,
    coarse: PyReadonlyArray2<f64>,
    fine_rows: usize,
    fine_cols: usize,
) -> PyResult<&'py PyArray2<f64>> {
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

    Ok(fine.into_pyarray(py))
}

/// Compute Laplacian difference matrix products for smoothness regularization
///
/// Computes D^T @ D @ x for both x and y directions where D is the first-order
/// difference matrix. This is used in the smoothness cost function:
/// J_smooth = 0.5 * x^T @ (Dx^T @ Dx + Dy^T @ Dy) @ x / sigma^2
#[pyfunction]
pub fn apply_laplacian<'py>(
    py: Python<'py>,
    data: PyReadonlyArray2<f64>,
) -> PyResult<&'py PyArray2<f64>> {
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

    Ok(result.into_pyarray(py))
}

/// Evaluate one grid-search candidate cost map from RT coefficients.
///
/// Computes per-pixel:
/// J = J_obs + J_prior
///
/// with observation term aggregated over bands using provided spectral weights.
#[pyfunction]
pub fn evaluate_grid_search_candidate_cost<'py>(
    py: Python<'py>,
    toa: PyReadonlyArray3<f32>,
    xap: PyReadonlyArray3<f32>,
    xbp: PyReadonlyArray3<f32>,
    xcp: PyReadonlyArray3<f32>,
    boa_prior: PyReadonlyArray3<f32>,
    boa_unc: PyReadonlyArray3<f32>,
    band_weights: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    aot_val: f32,
    tcwv_val: f32,
    aot_prior: PyReadonlyArray2<f32>,
    tcwv_prior: PyReadonlyArray2<f32>,
    aot_prior_unc: PyReadonlyArray2<f32>,
    tcwv_prior_unc: PyReadonlyArray2<f32>,
) -> PyResult<&'py PyArray2<f32>> {
    let toa = toa.as_array();
    let xap = xap.as_array();
    let xbp = xbp.as_array();
    let xcp = xcp.as_array();
    let boa_prior = boa_prior.as_array();
    let boa_unc = boa_unc.as_array();
    let band_weights = band_weights.as_array();
    let valid = valid_mask.as_array();
    let aot_prior = aot_prior.as_array();
    let tcwv_prior = tcwv_prior.as_array();
    let aot_prior_unc = aot_prior_unc.as_array();
    let tcwv_prior_unc = tcwv_prior_unc.as_array();

    let toa_shape = toa.shape();
    if toa_shape.len() != 3 {
        return Err(PyValueError::new_err("toa must be 3D (band, y, x)"));
    }
    let (n_band, ny, nx) = (toa_shape[0], toa_shape[1], toa_shape[2]);

    for (name, arr) in [
        ("xap", xap.shape()),
        ("xbp", xbp.shape()),
        ("xcp", xcp.shape()),
        ("boa_prior", boa_prior.shape()),
        ("boa_unc", boa_unc.shape()),
    ] {
        if arr != [n_band, ny, nx] {
            return Err(PyValueError::new_err(format!(
                "{name} must have shape (band, y, x) matching toa"
            )));
        }
    }
    if band_weights.len() != n_band {
        return Err(PyValueError::new_err(
            "band_weights length must match toa.shape[0]",
        ));
    }
    if valid.shape() != [ny, nx]
        || aot_prior.shape() != [ny, nx]
        || tcwv_prior.shape() != [ny, nx]
        || aot_prior_unc.shape() != [ny, nx]
        || tcwv_prior_unc.shape() != [ny, nx]
    {
        return Err(PyValueError::new_err(
            "2D inputs must all match shape (y, x)",
        ));
    }

    let mut cost = Array2::<f32>::zeros((ny, nx));

    // Parallel observation cost: compute per-row slices in parallel.
    let row_costs: Vec<Vec<f32>> = (0..ny)
        .into_par_iter()
        .map(|iy| {
            let mut row = vec![0.0f32; nx];
            for ix in 0..nx {
                if !valid[[iy, ix]] {
                    continue;
                }
                let mut pixel_cost = 0.0_f64;
                for ib in 0..n_band {
                    let band_w = band_weights[ib] as f64;
                    let unc = boa_unc[[ib, iy, ix]] as f64;
                    if !unc.is_finite() || unc <= 0.0 {
                        continue;
                    }

                    let toa_v = toa[[ib, iy, ix]] as f64;
                    let xap_v = xap[[ib, iy, ix]] as f64;
                    let xbp_v = xbp[[ib, iy, ix]] as f64;
                    let xcp_v = xcp[[ib, iy, ix]] as f64;
                    let prior_v = boa_prior[[ib, iy, ix]] as f64;
                    if !toa_v.is_finite()
                        || !xap_v.is_finite()
                        || !xbp_v.is_finite()
                        || !xcp_v.is_finite()
                        || !prior_v.is_finite()
                    {
                        continue;
                    }

                    let y = xap_v * toa_v - xbp_v;
                    let denom = 1.0 + xcp_v * y;
                    if !denom.is_finite() || denom.abs() < 1e-12 {
                        continue;
                    }

                    let boa_model = y / denom;
                    let diff = boa_model - prior_v;
                    if !diff.is_finite() {
                        continue;
                    }

                    let w = band_w / (unc * unc).max(1e-12);
                    if !w.is_finite() {
                        continue;
                    }
                    let add = 0.5 * w * diff * diff;
                    if add.is_finite() {
                        pixel_cost += add;
                    }
                }
                row[ix] = pixel_cost as f32;
            }
            row
        })
        .collect();

    for (iy, row) in row_costs.into_iter().enumerate() {
        for (ix, val) in row.into_iter().enumerate() {
            cost[[iy, ix]] = val;
        }
    }

    let aot_val = aot_val as f64;
    let tcwv_val = tcwv_val as f64;

    // Parallel prior cost: compute per-row in parallel.
    let prior_costs: Vec<Vec<f32>> = (0..ny)
        .into_par_iter()
        .map(|iy| {
            let mut row = vec![0.0f32; nx];
            for ix in 0..nx {
                if !valid[[iy, ix]] {
                    continue;
                }

                let aot_p = aot_prior[[iy, ix]] as f64;
                let tcwv_p = tcwv_prior[[iy, ix]] as f64;
                let aot_u = (aot_prior_unc[[iy, ix]] as f64).abs();
                let tcwv_u = (tcwv_prior_unc[[iy, ix]] as f64).abs();
                if !aot_p.is_finite()
                    || !tcwv_p.is_finite()
                    || !aot_u.is_finite()
                    || !tcwv_u.is_finite()
                {
                    continue;
                }

                let prior = 0.5
                    * (((aot_val - aot_p) * (aot_val - aot_p)) / (aot_u * aot_u).max(1e-12)
                        + ((tcwv_val - tcwv_p) * (tcwv_val - tcwv_p))
                            / (tcwv_u * tcwv_u).max(1e-12));
                if prior.is_finite() {
                    row[ix] = prior as f32;
                }
            }
            row
        })
        .collect();

    for (iy, row) in prior_costs.into_iter().enumerate() {
        for (ix, val) in row.into_iter().enumerate() {
            cost[[iy, ix]] += val;
        }
    }

    Ok(cost.into_pyarray(py))
}

/// Evaluate full AOT/TCWV candidate-cost cube with Rust-controlled loops.
///
/// The `coeff_provider(aot, tcwv)` callback must return a tuple of
/// `(xap, xbp, xcp)` arrays, each shaped `(band, y, x)` and `float32`.
#[pyfunction]
#[pyo3(signature = (
    coeff_provider,
    aot_axis,
    tcwv_axis,
    toa,
    boa_prior,
    boa_unc,
    band_weights,
    valid_mask,
    aot_prior,
    tcwv_prior,
    aot_prior_unc,
    tcwv_prior_unc,
    fixed_parameter = "none"
))]
pub fn evaluate_grid_search_cost_cube_with_provider<'py>(
    py: Python<'py>,
    coeff_provider: PyObject,
    aot_axis: PyReadonlyArray1<f32>,
    tcwv_axis: PyReadonlyArray1<f32>,
    toa: PyReadonlyArray3<f32>,
    boa_prior: PyReadonlyArray3<f32>,
    boa_unc: PyReadonlyArray3<f32>,
    band_weights: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    aot_prior: PyReadonlyArray2<f32>,
    tcwv_prior: PyReadonlyArray2<f32>,
    aot_prior_unc: PyReadonlyArray2<f32>,
    tcwv_prior_unc: PyReadonlyArray2<f32>,
    fixed_parameter: &str,
) -> PyResult<&'py PyArray4<f32>> {
    let fixed_parameter = parse_fixed_parameter(fixed_parameter)?;
    let (costs, _obs_counts) = compute_grid_search_cost_cube(
        py,
        coeff_provider,
        aot_axis,
        tcwv_axis,
        toa,
        boa_prior,
        boa_unc,
        band_weights,
        valid_mask,
        aot_prior,
        tcwv_prior,
        aot_prior_unc,
        tcwv_prior_unc,
        fixed_parameter,
    )?;

    Ok(costs.into_pyarray(py))
}

#[pyfunction]
#[pyo3(signature = (
    coeff_provider,
    aot_axis,
    tcwv_axis,
    toa,
    boa_prior,
    boa_unc,
    band_weights,
    valid_mask,
    aot_prior,
    tcwv_prior,
    aot_prior_unc,
    tcwv_prior_unc,
    fixed_parameter = "none"
))]
pub fn evaluate_grid_search_cost_cube_with_provider_qa<'py>(
    py: Python<'py>,
    coeff_provider: PyObject,
    aot_axis: PyReadonlyArray1<f32>,
    tcwv_axis: PyReadonlyArray1<f32>,
    toa: PyReadonlyArray3<f32>,
    boa_prior: PyReadonlyArray3<f32>,
    boa_unc: PyReadonlyArray3<f32>,
    band_weights: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    aot_prior: PyReadonlyArray2<f32>,
    tcwv_prior: PyReadonlyArray2<f32>,
    aot_prior_unc: PyReadonlyArray2<f32>,
    tcwv_prior_unc: PyReadonlyArray2<f32>,
    fixed_parameter: &str,
) -> PyResult<(&'py PyArray4<f32>, &'py PyArray4<u16>)> {
    let fixed_parameter = parse_fixed_parameter(fixed_parameter)?;
    let (costs, obs_counts) = compute_grid_search_cost_cube(
        py,
        coeff_provider,
        aot_axis,
        tcwv_axis,
        toa,
        boa_prior,
        boa_unc,
        band_weights,
        valid_mask,
        aot_prior,
        tcwv_prior,
        aot_prior_unc,
        tcwv_prior_unc,
        fixed_parameter,
    )?;

    Ok((costs.into_pyarray(py), obs_counts.into_pyarray(py)))
}

#[inline]
fn evaluate_candidate_pixel_cost(
    iy: usize,
    ix: usize,
    coeff_y: usize,
    coeff_x: usize,
    n_band: usize,
    aot_val_f64: f64,
    tcwv_val_f64: f64,
    toa: ArrayView3<'_, f32>,
    boa_prior: ArrayView3<'_, f32>,
    boa_unc: ArrayView3<'_, f32>,
    band_weights: ArrayView1<'_, f32>,
    valid: ArrayView2<'_, bool>,
    aot_prior: ArrayView2<'_, f32>,
    tcwv_prior: ArrayView2<'_, f32>,
    aot_prior_unc: ArrayView2<'_, f32>,
    tcwv_prior_unc: ArrayView2<'_, f32>,
    xap: ArrayView3<'_, f32>,
    xbp: ArrayView3<'_, f32>,
    xcp: ArrayView3<'_, f32>,
    fixed_parameter: FixedParameter,
) -> (f32, u16) {
    if !valid[[iy, ix]] {
        return (0.0, 0);
    }

    let mut pixel_cost = 0.0_f64;
    let mut pixel_obs_count = 0u16;

    for ib in 0..n_band {
        let band_w = band_weights[ib] as f64;
        let unc = boa_unc[[ib, iy, ix]] as f64;
        if !unc.is_finite() || unc <= 0.0 {
            continue;
        }

        let toa_v = toa[[ib, iy, ix]] as f64;
        let xap_v = xap[[ib, coeff_y, coeff_x]] as f64;
        let xbp_v = xbp[[ib, coeff_y, coeff_x]] as f64;
        let xcp_v = xcp[[ib, coeff_y, coeff_x]] as f64;
        let prior_v = boa_prior[[ib, iy, ix]] as f64;
        if !toa_v.is_finite()
            || !xap_v.is_finite()
            || !xbp_v.is_finite()
            || !xcp_v.is_finite()
            || !prior_v.is_finite()
        {
            continue;
        }

        let y = xap_v * toa_v - xbp_v;
        let denom = 1.0 + xcp_v * y;
        if !denom.is_finite() || denom.abs() < 1e-12 {
            continue;
        }

        let boa_model = y / denom;
        let diff = boa_model - prior_v;
        if !diff.is_finite() {
            continue;
        }

        let w = band_w / (unc * unc).max(1e-12);
        if !w.is_finite() || w <= 0.0 {
            continue;
        }
        let add = 0.5 * w * diff * diff;
        if add.is_finite() {
            pixel_cost += add;
            pixel_obs_count = pixel_obs_count.saturating_add(1);
        }
    }

    let aot_p = aot_prior[[iy, ix]] as f64;
    let tcwv_p = tcwv_prior[[iy, ix]] as f64;
    let aot_u = (aot_prior_unc[[iy, ix]] as f64).abs();
    let tcwv_u = (tcwv_prior_unc[[iy, ix]] as f64).abs();
    match fixed_parameter {
        FixedParameter::None => {
            if aot_p.is_finite() && tcwv_p.is_finite() && aot_u.is_finite() && tcwv_u.is_finite() {
                let prior = 0.5
                    * (((aot_val_f64 - aot_p) * (aot_val_f64 - aot_p))
                        / (aot_u * aot_u).max(1e-12)
                        + ((tcwv_val_f64 - tcwv_p) * (tcwv_val_f64 - tcwv_p))
                            / (tcwv_u * tcwv_u).max(1e-12));
                if prior.is_finite() {
                    pixel_cost += prior;
                }
            }
        }
        FixedParameter::Aot => {
            if tcwv_p.is_finite() && tcwv_u.is_finite() {
                let prior = 0.5 * ((tcwv_val_f64 - tcwv_p) * (tcwv_val_f64 - tcwv_p))
                    / (tcwv_u * tcwv_u).max(1e-12);
                if prior.is_finite() {
                    pixel_cost += prior;
                }
            }
        }
        FixedParameter::Tcwv => {
            if aot_p.is_finite() && aot_u.is_finite() {
                let prior = 0.5 * ((aot_val_f64 - aot_p) * (aot_val_f64 - aot_p))
                    / (aot_u * aot_u).max(1e-12);
                if prior.is_finite() {
                    pixel_cost += prior;
                }
            }
        }
    }

    (pixel_cost as f32, pixel_obs_count)
}

#[pyfunction]
#[pyo3(signature = (
    coeff_provider,
    aot_axis,
    tcwv_axis,
    toa,
    boa_prior,
    boa_unc,
    band_weights,
    valid_mask,
    aot_prior,
    tcwv_prior,
    aot_prior_unc,
    tcwv_prior_unc,
    block_size,
    fixed_parameter = "none"
))]
pub fn evaluate_block_grid_search_cost_cube_with_provider_qa<'py>(
    py: Python<'py>,
    coeff_provider: PyObject,
    aot_axis: PyReadonlyArray1<f32>,
    tcwv_axis: PyReadonlyArray1<f32>,
    toa: PyReadonlyArray3<f32>,
    boa_prior: PyReadonlyArray3<f32>,
    boa_unc: PyReadonlyArray3<f32>,
    band_weights: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    aot_prior: PyReadonlyArray2<f32>,
    tcwv_prior: PyReadonlyArray2<f32>,
    aot_prior_unc: PyReadonlyArray2<f32>,
    tcwv_prior_unc: PyReadonlyArray2<f32>,
    block_size: usize,
    fixed_parameter: &str,
) -> PyResult<(&'py PyArray4<f32>, &'py PyArray4<u16>, &'py PyArray2<bool>)> {
    if block_size == 0 {
        return Err(PyValueError::new_err("block_size must be positive"));
    }
    let fixed_parameter = parse_fixed_parameter(fixed_parameter)?;

    let aot_axis = aot_axis.as_array();
    let tcwv_axis = tcwv_axis.as_array();
    let toa = toa.as_array();
    let boa_prior = boa_prior.as_array();
    let boa_unc = boa_unc.as_array();
    let band_weights = band_weights.as_array();
    let valid = valid_mask.as_array();
    let aot_prior = aot_prior.as_array();
    let tcwv_prior = tcwv_prior.as_array();
    let aot_prior_unc = aot_prior_unc.as_array();
    let tcwv_prior_unc = tcwv_prior_unc.as_array();

    if aot_axis.is_empty() || tcwv_axis.is_empty() {
        return Err(PyValueError::new_err(
            "aot_axis and tcwv_axis must be non-empty",
        ));
    }
    validate_fixed_axes(fixed_parameter, aot_axis.len(), tcwv_axis.len())?;

    let toa_shape = toa.shape();
    if toa_shape.len() != 3 {
        return Err(PyValueError::new_err("toa must be 3D (band, y, x)"));
    }
    let (n_band, ny, nx) = (toa_shape[0], toa_shape[1], toa_shape[2]);

    for (name, arr) in [
        ("boa_prior", boa_prior.shape()),
        ("boa_unc", boa_unc.shape()),
    ] {
        if arr != [n_band, ny, nx] {
            return Err(PyValueError::new_err(format!(
                "{name} must have shape (band, y, x) matching toa"
            )));
        }
    }
    if band_weights.len() != n_band {
        return Err(PyValueError::new_err(
            "band_weights length must match toa.shape[0]",
        ));
    }
    if valid.shape() != [ny, nx]
        || aot_prior.shape() != [ny, nx]
        || tcwv_prior.shape() != [ny, nx]
        || aot_prior_unc.shape() != [ny, nx]
        || tcwv_prior_unc.shape() != [ny, nx]
    {
        return Err(PyValueError::new_err(
            "2D inputs must all match shape (y, x)",
        ));
    }

    let block_rows = (ny + block_size - 1) / block_size;
    let block_cols = (nx + block_size - 1) / block_size;
    let mut costs = Array4::<f32>::zeros((aot_axis.len(), tcwv_axis.len(), block_rows, block_cols));
    let mut obs_counts =
        Array4::<u16>::zeros((aot_axis.len(), tcwv_axis.len(), block_rows, block_cols));
    let mut block_valid = Array2::<bool>::from_elem((block_rows, block_cols), false);
    for iy in 0..ny {
        for ix in 0..nx {
            if valid[[iy, ix]] {
                block_valid[[iy / block_size, ix / block_size]] = true;
            }
        }
    }

    for (ia, aot_val) in aot_axis.iter().enumerate() {
        for (it, tcwv_val) in tcwv_axis.iter().enumerate() {
            let returned = coeff_provider.call1(py, (*aot_val, *tcwv_val))?;
            let (xap_obj, xbp_obj, xcp_obj): (&PyArray3<f32>, &PyArray3<f32>, &PyArray3<f32>) =
                returned.extract(py)?;

            let xap_ro = xap_obj.readonly();
            let xbp_ro = xbp_obj.readonly();
            let xcp_ro = xcp_obj.readonly();
            let xap = xap_ro.as_array();
            let xbp = xbp_ro.as_array();
            let xcp = xcp_ro.as_array();

            let full_coeff_shape = [n_band, ny, nx];
            let block_coeff_shape = [n_band, block_rows, block_cols];
            let full_coeffs = xap.shape() == full_coeff_shape
                && xbp.shape() == full_coeff_shape
                && xcp.shape() == full_coeff_shape;
            let block_coeffs = xap.shape() == block_coeff_shape
                && xbp.shape() == block_coeff_shape
                && xcp.shape() == block_coeff_shape;
            if !full_coeffs && !block_coeffs {
                return Err(PyValueError::new_err(
                    "coeff_provider returned coefficients with invalid shape; expected all coefficient arrays to match toa (band, y, x) or the block grid",
                ));
            }
            let coeffs_are_block_grid = block_coeffs;

            let aot_val_f64 = *aot_val as f64;
            let tcwv_val_f64 = *tcwv_val as f64;
            let mut cost_slice = costs.slice_mut(s![ia, it, .., ..]);
            let mut count_slice = obs_counts.slice_mut(s![ia, it, .., ..]);
            cost_slice
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .zip(count_slice.axis_iter_mut(Axis(0)).into_par_iter())
                .enumerate()
                .for_each(|(by, (mut cost_row, mut count_row))| {
                    let y0 = by * block_size;
                    let y1 = ((by + 1) * block_size).min(ny);
                    for bx in 0..block_cols {
                        let x0 = bx * block_size;
                        let x1 = ((bx + 1) * block_size).min(nx);
                        let mut block_cost = 0.0_f64;
                        let mut block_obs_count = 0u16;

                        for iy in y0..y1 {
                            for ix in x0..x1 {
                                let coeff_y = if coeffs_are_block_grid { by } else { iy };
                                let coeff_x = if coeffs_are_block_grid { bx } else { ix };
                                let (pixel_cost, pixel_obs_count) = evaluate_candidate_pixel_cost(
                                    iy,
                                    ix,
                                    coeff_y,
                                    coeff_x,
                                    n_band,
                                    aot_val_f64,
                                    tcwv_val_f64,
                                    toa,
                                    boa_prior,
                                    boa_unc,
                                    band_weights,
                                    valid,
                                    aot_prior,
                                    tcwv_prior,
                                    aot_prior_unc,
                                    tcwv_prior_unc,
                                    xap,
                                    xbp,
                                    xcp,
                                    fixed_parameter,
                                );
                                block_cost += pixel_cost as f64;
                                block_obs_count = block_obs_count.saturating_add(pixel_obs_count);
                            }
                        }

                        cost_row[bx] = block_cost as f32;
                        count_row[bx] = block_obs_count;
                    }
                });
        }
    }

    Ok((
        costs.into_pyarray(py),
        obs_counts.into_pyarray(py),
        block_valid.into_pyarray(py),
    ))
}

fn compute_grid_search_cost_cube<'py>(
    py: Python<'py>,
    coeff_provider: PyObject,
    aot_axis: PyReadonlyArray1<f32>,
    tcwv_axis: PyReadonlyArray1<f32>,
    toa: PyReadonlyArray3<f32>,
    boa_prior: PyReadonlyArray3<f32>,
    boa_unc: PyReadonlyArray3<f32>,
    band_weights: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    aot_prior: PyReadonlyArray2<f32>,
    tcwv_prior: PyReadonlyArray2<f32>,
    aot_prior_unc: PyReadonlyArray2<f32>,
    tcwv_prior_unc: PyReadonlyArray2<f32>,
    fixed_parameter: FixedParameter,
) -> PyResult<(Array4<f32>, Array4<u16>)> {
    let aot_axis = aot_axis.as_array();
    let tcwv_axis = tcwv_axis.as_array();
    let toa = toa.as_array();
    let boa_prior = boa_prior.as_array();
    let boa_unc = boa_unc.as_array();
    let band_weights = band_weights.as_array();
    let valid = valid_mask.as_array();
    let aot_prior = aot_prior.as_array();
    let tcwv_prior = tcwv_prior.as_array();
    let aot_prior_unc = aot_prior_unc.as_array();
    let tcwv_prior_unc = tcwv_prior_unc.as_array();

    if aot_axis.is_empty() || tcwv_axis.is_empty() {
        return Err(PyValueError::new_err(
            "aot_axis and tcwv_axis must be non-empty",
        ));
    }
    validate_fixed_axes(fixed_parameter, aot_axis.len(), tcwv_axis.len())?;

    let toa_shape = toa.shape();
    if toa_shape.len() != 3 {
        return Err(PyValueError::new_err("toa must be 3D (band, y, x)"));
    }
    let (n_band, ny, nx) = (toa_shape[0], toa_shape[1], toa_shape[2]);

    for (name, arr) in [
        ("boa_prior", boa_prior.shape()),
        ("boa_unc", boa_unc.shape()),
    ] {
        if arr != [n_band, ny, nx] {
            return Err(PyValueError::new_err(format!(
                "{name} must have shape (band, y, x) matching toa"
            )));
        }
    }
    if band_weights.len() != n_band {
        return Err(PyValueError::new_err(
            "band_weights length must match toa.shape[0]",
        ));
    }
    if valid.shape() != [ny, nx]
        || aot_prior.shape() != [ny, nx]
        || tcwv_prior.shape() != [ny, nx]
        || aot_prior_unc.shape() != [ny, nx]
        || tcwv_prior_unc.shape() != [ny, nx]
    {
        return Err(PyValueError::new_err(
            "2D inputs must all match shape (y, x)",
        ));
    }

    let mut costs = Array4::<f32>::zeros((aot_axis.len(), tcwv_axis.len(), ny, nx));
    let mut obs_counts = Array4::<u16>::zeros((aot_axis.len(), tcwv_axis.len(), ny, nx));

    for (ia, aot_val) in aot_axis.iter().enumerate() {
        for (it, tcwv_val) in tcwv_axis.iter().enumerate() {
            let returned = coeff_provider.call1(py, (*aot_val, *tcwv_val))?;
            let (xap_obj, xbp_obj, xcp_obj): (&PyArray3<f32>, &PyArray3<f32>, &PyArray3<f32>) =
                returned.extract(py)?;

            let xap_ro = xap_obj.readonly();
            let xbp_ro = xbp_obj.readonly();
            let xcp_ro = xcp_obj.readonly();
            let xap = xap_ro.as_array();
            let xbp = xbp_ro.as_array();
            let xcp = xcp_ro.as_array();

            for (name, arr) in [
                ("xap", xap.shape()),
                ("xbp", xbp.shape()),
                ("xcp", xcp.shape()),
            ] {
                if arr != [n_band, ny, nx] {
                    return Err(PyValueError::new_err(format!(
                        "coeff_provider returned {name} with invalid shape; expected (band, y, x) matching toa"
                    )));
                }
            }

            // Compute observation + prior cost per-row in parallel for this (aot, tcwv) slice.
            let aot_val_f64 = *aot_val as f64;
            let tcwv_val_f64 = *tcwv_val as f64;
            let mut cost_slice = costs.slice_mut(s![ia, it, .., ..]);
            let mut count_slice = obs_counts.slice_mut(s![ia, it, .., ..]);
            cost_slice
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .zip(count_slice.axis_iter_mut(Axis(0)).into_par_iter())
                .enumerate()
                .for_each(|(iy, (mut cost_row, mut count_row))| {
                    for ix in 0..nx {
                        let (pixel_cost, pixel_obs_count) = evaluate_candidate_pixel_cost(
                            iy,
                            ix,
                            iy,
                            ix,
                            n_band,
                            aot_val_f64,
                            tcwv_val_f64,
                            toa,
                            boa_prior,
                            boa_unc,
                            band_weights,
                            valid,
                            aot_prior,
                            tcwv_prior,
                            aot_prior_unc,
                            tcwv_prior_unc,
                            xap,
                            xbp,
                            xcp,
                            fixed_parameter,
                        );
                        cost_row[ix] = pixel_cost;
                        count_row[ix] = pixel_obs_count;
                    }
                });
        }
    }

    Ok((costs, obs_counts))
}

/// Refine per-pixel AOT/TCWV estimates from a discrete cost cube.
///
/// The routine finds the minimum over (aot, tcwv) candidates for each pixel,
/// then performs a local quadratic Newton step using 3x3 finite differences
/// around the minimum (when interior). Uncertainty is approximated from the
/// inverse Hessian diagonal.
#[pyfunction]
#[pyo3(signature = (costs, aot_axis, tcwv_axis, valid_mask, fixed_parameter = "none"))]
pub fn quadratic_refine_grid_search<'py>(
    py: Python<'py>,
    costs: PyReadonlyArray4<f32>,
    aot_axis: PyReadonlyArray1<f32>,
    tcwv_axis: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    fixed_parameter: &str,
) -> PyResult<(
    &'py PyArray2<f32>,
    &'py PyArray2<f32>,
    &'py PyArray2<f32>,
    &'py PyArray2<f32>,
)> {
    let (
        aot_best,
        tcwv_best,
        aot_unc,
        tcwv_unc,
        _invalid,
        _boundary,
        _lower_aot_boundary,
        _zero_obs,
    ) = refine_grid_search_with_qa(
        costs.as_array(),
        None,
        aot_axis.as_array(),
        tcwv_axis.as_array(),
        valid_mask.as_array(),
        parse_fixed_parameter(fixed_parameter)?,
    )?;

    Ok((
        aot_best.into_pyarray(py),
        tcwv_best.into_pyarray(py),
        aot_unc.into_pyarray(py),
        tcwv_unc.into_pyarray(py),
    ))
}

#[pyfunction]
#[pyo3(signature = (costs, obs_counts, aot_axis, tcwv_axis, valid_mask, fixed_parameter = "none"))]
pub fn quadratic_refine_grid_search_qa<'py>(
    py: Python<'py>,
    costs: PyReadonlyArray4<f32>,
    obs_counts: PyReadonlyArray4<u16>,
    aot_axis: PyReadonlyArray1<f32>,
    tcwv_axis: PyReadonlyArray1<f32>,
    valid_mask: PyReadonlyArray2<bool>,
    fixed_parameter: &str,
) -> PyResult<(
    &'py PyArray2<f32>,
    &'py PyArray2<f32>,
    &'py PyArray2<f32>,
    &'py PyArray2<f32>,
    &'py PyArray2<bool>,
    &'py PyArray2<bool>,
    &'py PyArray2<bool>,
    &'py PyArray2<bool>,
)> {
    let (
        aot_best,
        tcwv_best,
        aot_unc,
        tcwv_unc,
        invalid_mask,
        boundary_mask,
        lower_aot_boundary_mask,
        zero_obs_mask,
    ) = refine_grid_search_with_qa(
        costs.as_array(),
        Some(obs_counts.as_array()),
        aot_axis.as_array(),
        tcwv_axis.as_array(),
        valid_mask.as_array(),
        parse_fixed_parameter(fixed_parameter)?,
    )?;

    Ok((
        aot_best.into_pyarray(py),
        tcwv_best.into_pyarray(py),
        aot_unc.into_pyarray(py),
        tcwv_unc.into_pyarray(py),
        invalid_mask.into_pyarray(py),
        boundary_mask.into_pyarray(py),
        lower_aot_boundary_mask.into_pyarray(py),
        zero_obs_mask.into_pyarray(py),
    ))
}

fn refine_grid_search_with_qa(
    costs: ArrayView4<'_, f32>,
    obs_counts: Option<ArrayView4<'_, u16>>,
    aot_axis: ArrayView1<'_, f32>,
    tcwv_axis: ArrayView1<'_, f32>,
    valid: ArrayView2<'_, bool>,
    fixed_parameter: FixedParameter,
) -> PyResult<(
    Array2<f32>,
    Array2<f32>,
    Array2<f32>,
    Array2<f32>,
    Array2<bool>,
    Array2<bool>,
    Array2<bool>,
    Array2<bool>,
)> {
    let shape = costs.shape();
    let (n_aot, n_tcwv, ny, nx) = (shape[0], shape[1], shape[2], shape[3]);

    if aot_axis.len() != n_aot {
        return Err(PyValueError::new_err(
            "aot_axis length must match costs.shape[0]",
        ));
    }
    if tcwv_axis.len() != n_tcwv {
        return Err(PyValueError::new_err(
            "tcwv_axis length must match costs.shape[1]",
        ));
    }
    validate_fixed_axes(fixed_parameter, n_aot, n_tcwv)?;
    if valid.shape() != [ny, nx] {
        return Err(PyValueError::new_err(
            "valid_mask shape must match costs.shape[2:4]",
        ));
    }
    if let Some(obs_counts) = obs_counts {
        if obs_counts.shape() != [n_aot, n_tcwv, ny, nx] {
            return Err(PyValueError::new_err(
                "obs_counts shape must match costs.shape",
            ));
        }
    }

    let mut aot_best = Array2::<f32>::zeros((ny, nx));
    let mut tcwv_best = Array2::<f32>::zeros((ny, nx));
    let mut aot_unc = Array2::<f32>::zeros((ny, nx));
    let mut tcwv_unc = Array2::<f32>::zeros((ny, nx));
    let mut invalid_mask = Array2::<bool>::from_elem((ny, nx), false);
    let mut boundary_mask = Array2::<bool>::from_elem((ny, nx), false);
    let mut lower_aot_boundary_mask = Array2::<bool>::from_elem((ny, nx), false);
    let mut zero_obs_mask = Array2::<bool>::from_elem((ny, nx), false);

    let da = if n_aot > 1 {
        (aot_axis[1] - aot_axis[0]).abs().max(1e-6)
    } else {
        0.05
    };
    let dt = if n_tcwv > 1 {
        (tcwv_axis[1] - tcwv_axis[0]).abs().max(1e-6)
    } else {
        0.2
    };
    aot_unc.fill(da);
    tcwv_unc.fill(dt);

    for iy in 0..ny {
        for ix in 0..nx {
            if !valid[[iy, ix]] {
                continue;
            }

            // Discrete minimum.
            let mut best_ia = 0usize;
            let mut best_it = 0usize;
            let mut best_cost = f32::INFINITY;
            for ia in 0..n_aot {
                for it in 0..n_tcwv {
                    let c = costs[[ia, it, iy, ix]];
                    if c.is_finite() && c < best_cost {
                        best_cost = c;
                        best_ia = ia;
                        best_it = it;
                    }
                }
            }

            aot_best[[iy, ix]] = aot_axis[best_ia];
            tcwv_best[[iy, ix]] = tcwv_axis[best_it];

            if !best_cost.is_finite() {
                invalid_mask[[iy, ix]] = true;
                continue;
            }

            if let Some(obs_counts) = obs_counts {
                if obs_counts[[best_ia, best_it, iy, ix]] == 0 {
                    invalid_mask[[iy, ix]] = true;
                    zero_obs_mask[[iy, ix]] = true;
                    continue;
                }
            }

            if fixed_parameter == FixedParameter::Tcwv {
                if best_ia == 0 || best_ia + 1 >= n_aot {
                    boundary_mask[[iy, ix]] = true;
                    if best_ia == 0 {
                        lower_aot_boundary_mask[[iy, ix]] = true;
                    }
                    continue;
                }

                let ia = best_ia;
                let f00 = costs[[ia, 0, iy, ix]] as f64;
                let fxm = costs[[ia - 1, 0, iy, ix]] as f64;
                let fxp = costs[[ia + 1, 0, iy, ix]] as f64;
                if !f00.is_finite() || !fxm.is_finite() || !fxp.is_finite() {
                    continue;
                }

                let da64 = da as f64;
                let dfdx = (fxp - fxm) / (2.0 * da64);
                let d2fdx2 = (fxp - 2.0 * f00 + fxm) / (da64 * da64);
                if d2fdx2 > 1e-6 {
                    let delta_a = dfdx / d2fdx2;
                    let mut a_fit = (aot_axis[ia] as f64) - delta_a;
                    let a_lo = aot_axis[ia - 1] as f64;
                    let a_hi = aot_axis[ia + 1] as f64;
                    a_fit = a_fit.clamp(a_lo.min(a_hi), a_lo.max(a_hi));
                    aot_best[[iy, ix]] = a_fit as f32;
                    aot_unc[[iy, ix]] = (1.0 / d2fdx2).max(1e-12).sqrt() as f32;
                }
                continue;
            }

            if fixed_parameter == FixedParameter::Aot {
                if best_it == 0 || best_it + 1 >= n_tcwv {
                    boundary_mask[[iy, ix]] = true;
                    continue;
                }

                let it = best_it;
                let f00 = costs[[0, it, iy, ix]] as f64;
                let fym = costs[[0, it - 1, iy, ix]] as f64;
                let fyp = costs[[0, it + 1, iy, ix]] as f64;
                if !f00.is_finite() || !fym.is_finite() || !fyp.is_finite() {
                    continue;
                }

                let dt64 = dt as f64;
                let dfdy = (fyp - fym) / (2.0 * dt64);
                let d2fdy2 = (fyp - 2.0 * f00 + fym) / (dt64 * dt64);
                if d2fdy2 > 1e-6 {
                    let delta_t = dfdy / d2fdy2;
                    let mut t_fit = (tcwv_axis[it] as f64) - delta_t;
                    let t_lo = tcwv_axis[it - 1] as f64;
                    let t_hi = tcwv_axis[it + 1] as f64;
                    t_fit = t_fit.clamp(t_lo.min(t_hi), t_lo.max(t_hi));
                    tcwv_best[[iy, ix]] = t_fit as f32;
                    tcwv_unc[[iy, ix]] = (1.0 / d2fdy2).max(1e-12).sqrt() as f32;
                }
                continue;
            }

            // Need interior point for finite-difference Hessian.
            if best_ia == 0 || best_ia + 1 >= n_aot || best_it == 0 || best_it + 1 >= n_tcwv {
                boundary_mask[[iy, ix]] = true;
                if best_ia == 0 {
                    lower_aot_boundary_mask[[iy, ix]] = true;
                }
                continue;
            }

            let ia = best_ia;
            let it = best_it;

            let f00 = costs[[ia, it, iy, ix]] as f64;
            let fxm = costs[[ia - 1, it, iy, ix]] as f64;
            let fxp = costs[[ia + 1, it, iy, ix]] as f64;
            let fym = costs[[ia, it - 1, iy, ix]] as f64;
            let fyp = costs[[ia, it + 1, iy, ix]] as f64;
            let fmm = costs[[ia - 1, it - 1, iy, ix]] as f64;
            let fmp = costs[[ia - 1, it + 1, iy, ix]] as f64;
            let fpm = costs[[ia + 1, it - 1, iy, ix]] as f64;
            let fpp = costs[[ia + 1, it + 1, iy, ix]] as f64;

            if !f00.is_finite()
                || !fxm.is_finite()
                || !fxp.is_finite()
                || !fym.is_finite()
                || !fyp.is_finite()
                || !fmm.is_finite()
                || !fmp.is_finite()
                || !fpm.is_finite()
                || !fpp.is_finite()
            {
                continue;
            }

            let da64 = da as f64;
            let dt64 = dt as f64;
            let dfdx = (fxp - fxm) / (2.0 * da64);
            let dfdy = (fyp - fym) / (2.0 * dt64);
            let d2fdx2 = (fxp - 2.0 * f00 + fxm) / (da64 * da64);
            let d2fdy2 = (fyp - 2.0 * f00 + fym) / (dt64 * dt64);
            let d2fdxdy = (fpp - fpm - fmp + fmm) / (4.0 * da64 * dt64);

            let a = d2fdx2;
            let b = d2fdxdy;
            let c = d2fdy2;
            let det = a * c - b * b;
            // Use 1e-6 threshold: data originates from f32, so tighter bounds
            // cause spurious NaN/inf from ill-conditioned Hessians.
            if !(a > 1e-6 && c > 1e-6 && det > 1e-6) {
                continue;
            }

            // Newton update x* = x0 - H^{-1} g for [aot, tcwv].
            let delta_a = (c * dfdx - b * dfdy) / det;
            let delta_t = (-b * dfdx + a * dfdy) / det;
            let mut a_fit = (aot_axis[ia] as f64) - delta_a;
            let mut t_fit = (tcwv_axis[it] as f64) - delta_t;

            let a_lo = aot_axis[ia - 1] as f64;
            let a_hi = aot_axis[ia + 1] as f64;
            let t_lo = tcwv_axis[it - 1] as f64;
            let t_hi = tcwv_axis[it + 1] as f64;
            a_fit = a_fit.clamp(a_lo.min(a_hi), a_lo.max(a_hi));
            t_fit = t_fit.clamp(t_lo.min(t_hi), t_lo.max(t_hi));

            aot_best[[iy, ix]] = a_fit as f32;
            tcwv_best[[iy, ix]] = t_fit as f32;

            // Diagonal of inverse Hessian as variance proxy.
            let var_a = (c / det).max(1e-12);
            let var_t = (a / det).max(1e-12);
            aot_unc[[iy, ix]] = var_a.sqrt() as f32;
            tcwv_unc[[iy, ix]] = var_t.sqrt() as f32;
        }
    }

    Ok((
        aot_best,
        tcwv_best,
        aot_unc,
        tcwv_unc,
        invalid_mask,
        boundary_mask,
        lower_aot_boundary_mask,
        zero_obs_mask,
    ))
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
                let ci = ((2 * i + 1) * coarse_rows / (2 * high_rows)).min(coarse_rows - 1);
                let cj = ((2 * j + 1) * coarse_cols / (2 * high_cols)).min(coarse_cols - 1);
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
