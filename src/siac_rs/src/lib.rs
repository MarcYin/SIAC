//! SIAC Rust Extensions
//!
//! High-performance implementations of computationally intensive SIAC operations:
//! - BRDF kernel calculations (Ross-Thick, Li-Sparse)
//! - PSF convolution for scale matching
//! - Neural network emulator forward pass
//! - Multi-grid optimization utilities

use pyo3::prelude::*;

mod emulator;
mod kernels;
mod optimization;
mod psf;

/// SIAC Rust extension module
#[pymodule]
fn siac_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // BRDF kernels
    m.add_class::<kernels::RossThickLiSparse>()?;

    // PSF convolution
    m.add_class::<psf::PSFConvolver>()?;

    // Neural network emulator
    m.add_class::<emulator::TwoLayerNN>()?;

    // Optimization utilities
    m.add_function(wrap_pyfunction!(optimization::remap_to_coarse_grid, m)?)?;
    m.add_function(wrap_pyfunction!(optimization::interpolate_to_fine_grid, m)?)?;

    Ok(())
}
