//! Surface-driven aerosol solver kernel.
//!
//! The surface-driven retrieval inverts AOT directly from the surface-prior
//! mismatch: each AOT node's TOA->BOA correction is scored against the surface
//! prior, giving a per-node observation cost. That `(n_aot, y, x)` cost cube is
//! built (cheaply, vectorised) on the Python side; this kernel does the heavy,
//! loop-bound part — spatially median-pool each node slice, add a Gaussian
//! prior (CAMS/MAIAC) backstop, and pick the per-pixel arg-min over the AOT
//! axis. The windowed median + per-pixel sweep is what rayon parallelises here.
//!
//! It is a bit-exact (f64) replication of the reference numpy `_resolve_map`:
//!   pooled = rolling-median(window, center, min_periods) over each node;
//!   finite = all nodes pooled-finite;
//!   aod    = aot_axis[argmin_k(pooled_k + (aot_k - prior)^2 / unc^2)];
//!   jmin   = min_k pooled_k                 (obs-only, for cost reporting).
//! The pooling window matches pandas/xarray `center=True`: for window `W` the
//! window covering output `i` is `[i-(W-1)/2 .. i+W/2]`, clamped to bounds, and
//! a cell is valid only if it holds >= `pool_min_count` finite samples. TCWV is
//! held at the prior for the sweep (carried through by the Python solver).

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Median of the finite values in `buf` (NaN if none). Sorts in place.
/// Even counts average the two central values (numpy/xarray convention).
fn finite_median(buf: &mut Vec<f64>) -> f64 {
    buf.retain(|v| v.is_finite());
    if buf.is_empty() {
        return f64::NAN;
    }
    buf.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = buf.len();
    if n % 2 == 1 {
        buf[n / 2]
    } else {
        0.5 * (buf[n / 2 - 1] + buf[n / 2])
    }
}

/// Spatially median-pool a per-node surface-cost cube, add a Gaussian AOT-prior
/// backstop, and arg-min over the AOT axis per pixel. Bit-exact replication of
/// the reference numpy resolve (see module docs).
///
/// `cost_cube` is `(n_aot, y, x)` (NaN where a node has no usable observation).
/// `aot_prior`/`aot_prior_unc` are the per-pixel backstop centre and sigma; a
/// non-finite or non-positive sigma disables the backstop for that pixel.
/// `pool_window` is the full rolling window edge length (1 = no pooling). A node
/// contributes to a pixel only if its window holds >= `pool_min_count` finite
/// samples; a pixel yields a finite AOT only if *every* node is pooled-finite
/// there (mirrors numpy `isfinite(pooled).all(0)`).
///
/// Returns `(aot, aot_unc, jmin)`, each `(y, x)`; invalid pixels are NaN.
/// `aot_unc` is the half-width of the cost basin within `min_total + 0.5`,
/// floored at 0.02. `jmin` is the per-pixel min pooled obs-cost over nodes.
#[pyfunction]
#[pyo3(signature = (cost_cube, aot_axis, aot_prior, aot_prior_unc, valid_mask, pool_window, pool_min_count))]
pub fn surface_driven_pool_argmin<'py>(
    py: Python<'py>,
    cost_cube: PyReadonlyArray3<f64>,
    aot_axis: PyReadonlyArray1<f64>,
    aot_prior: PyReadonlyArray2<f64>,
    aot_prior_unc: PyReadonlyArray2<f64>,
    valid_mask: PyReadonlyArray2<bool>,
    pool_window: i64,
    pool_min_count: i64,
) -> PyResult<(&'py PyArray2<f64>, &'py PyArray2<f64>, &'py PyArray2<f64>)> {
    let cost = cost_cube.as_array();
    let aot_axis = aot_axis.as_array();
    let aot_prior = aot_prior.as_array();
    let aot_prior_unc = aot_prior_unc.as_array();
    let valid = valid_mask.as_array();

    let cube_shape = cost.shape();
    if cube_shape.len() != 3 {
        return Err(PyValueError::new_err("cost_cube must be 3D (n_aot, y, x)"));
    }
    let (n_aot, ny, nx) = (cube_shape[0], cube_shape[1], cube_shape[2]);
    if aot_axis.len() != n_aot {
        return Err(PyValueError::new_err(
            "aot_axis length must match cost_cube.shape[0]",
        ));
    }
    if n_aot == 0 {
        return Err(PyValueError::new_err("aot_axis must be non-empty"));
    }
    if valid.shape() != [ny, nx]
        || aot_prior.shape() != [ny, nx]
        || aot_prior_unc.shape() != [ny, nx]
    {
        return Err(PyValueError::new_err(
            "valid_mask, aot_prior, aot_prior_unc must all be (y, x)",
        ));
    }

    // pandas/xarray center=True window offsets for window length W.
    let w = pool_window.max(1) as usize;
    let lo = (w - 1) / 2;
    let hi = w / 2;
    let min_count = pool_min_count.max(1) as usize;

    let mut aot_out = Array2::<f64>::from_elem((ny, nx), f64::NAN);
    let mut unc_out = Array2::<f64>::from_elem((ny, nx), f64::NAN);
    let mut jmin_out = Array2::<f64>::from_elem((ny, nx), f64::NAN);

    py.allow_threads(|| {
        // 1. Pool each node slice (rolling median, center, min_periods). One
        //    (y, x) slab per AOT node; parallelised over nodes.
        let pooled: Vec<Array2<f64>> = (0..n_aot)
            .into_par_iter()
            .map(|k| {
                let mut slab = Array2::<f64>::from_elem((ny, nx), f64::NAN);
                let mut win: Vec<f64> = Vec::with_capacity(w * w);
                for iy in 0..ny {
                    let y0 = iy.saturating_sub(lo);
                    let y1 = (iy + hi + 1).min(ny);
                    for ix in 0..nx {
                        let x0 = ix.saturating_sub(lo);
                        let x1 = (ix + hi + 1).min(nx);
                        win.clear();
                        for jy in y0..y1 {
                            for jx in x0..x1 {
                                let v = cost[[k, jy, jx]];
                                if v.is_finite() {
                                    win.push(v);
                                }
                            }
                        }
                        if win.len() >= min_count {
                            slab[[iy, ix]] = finite_median(&mut win);
                        }
                    }
                }
                slab
            })
            .collect();

        // 2. Per-pixel: all-nodes-finite gate, arg-min(pooled + backstop), jmin.
        let rows: Vec<(usize, Vec<f64>, Vec<f64>, Vec<f64>)> = (0..ny)
            .into_par_iter()
            .map(|iy| {
                let mut aot_row = vec![f64::NAN; nx];
                let mut unc_row = vec![f64::NAN; nx];
                let mut jmin_row = vec![f64::NAN; nx];
                let mut totals = vec![f64::NAN; n_aot];
                for ix in 0..nx {
                    if !valid[[iy, ix]] {
                        continue;
                    }
                    let aot_p = aot_prior[[iy, ix]];
                    let aot_u = aot_prior_unc[[iy, ix]].abs();
                    let backstop_ok = aot_p.is_finite() && aot_u.is_finite() && aot_u > 0.0;

                    let mut all_finite = true;
                    let mut jmin = f64::INFINITY;
                    let mut best = f64::INFINITY;
                    let mut best_k = usize::MAX;
                    for k in 0..n_aot {
                        let p = pooled[k][[iy, ix]];
                        if !p.is_finite() {
                            all_finite = false;
                            totals[k] = f64::NAN;
                            continue;
                        }
                        if p < jmin {
                            jmin = p;
                        }
                        let mut total = p;
                        if backstop_ok {
                            let d = aot_axis[k] - aot_p;
                            total += d * d / (aot_u * aot_u).max(1e-12);
                        }
                        totals[k] = total;
                        if total < best {
                            best = total;
                            best_k = k;
                        }
                    }
                    // numpy medians aod_post only over all-nodes-finite pixels.
                    if !all_finite || best_k == usize::MAX {
                        continue;
                    }
                    aot_row[ix] = aot_axis[best_k];
                    jmin_row[ix] = jmin;
                    let thr = best + 0.5;
                    let mut blo = aot_axis[best_k];
                    let mut bhi = blo;
                    for k in 0..n_aot {
                        if totals[k].is_finite() && totals[k] <= thr {
                            let a = aot_axis[k];
                            blo = blo.min(a);
                            bhi = bhi.max(a);
                        }
                    }
                    unc_row[ix] = (0.5 * (bhi - blo)).max(0.02);
                }
                (iy, aot_row, unc_row, jmin_row)
            })
            .collect();

        for (iy, aot_row, unc_row, jmin_row) in rows {
            for ix in 0..nx {
                aot_out[[iy, ix]] = aot_row[ix];
                unc_out[[iy, ix]] = unc_row[ix];
                jmin_out[[iy, ix]] = jmin_row[ix];
            }
        }
    });

    Ok((
        aot_out.into_pyarray(py),
        unc_out.into_pyarray(py),
        jmin_out.into_pyarray(py),
    ))
}
