#![allow(non_local_definitions)]
//! BRDF Kernel Calculations
//!
//! Implements Ross-Thick and Li-Sparse BRDF kernels used in MODIS MCD43 products.

use ndarray::{Array2, Zip};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

/// Ross-Thick Li-Sparse kernel calculator
///
/// Implements the MODIS BRDF model kernels for surface reflectance modeling.
/// The model is: ρ = f0 + f1 * K_vol + f2 * K_geo
///
/// where K_vol is the Ross-Thick volumetric kernel and K_geo is the Li-Sparse
/// geometric kernel.
#[pyclass]
pub struct RossThickLiSparse {
    /// Height-to-breadth ratio for Li kernel (default: 2.0 for MODIS)
    hb: f64,
    /// Crown relative height for Li kernel (default: 1.0 for MODIS sparse)
    br: f64,
}

#[pymethods]
impl RossThickLiSparse {
    /// Create a new kernel calculator
    ///
    /// # Arguments
    /// * `hb` - Height-to-breadth ratio (default 2.0 for MODIS)
    /// * `br` - Crown relative height (default 1.0 for sparse vegetation)
    #[new]
    #[pyo3(signature = (hb=2.0, br=1.0))]
    fn new(hb: f64, br: f64) -> Self {
        Self { hb, br }
    }

    /// Compute Ross-Thick and Li-Sparse kernel values
    ///
    /// # Arguments
    /// * `vza` - View zenith angle in radians
    /// * `sza` - Solar zenith angle in radians
    /// * `raa` - Relative azimuth angle in radians
    ///
    /// # Returns
    /// Tuple of (Ross kernel array, Li kernel array)
    fn compute<'py>(
        &self,
        py: Python<'py>,
        vza: PyReadonlyArray2<f64>,
        sza: PyReadonlyArray2<f64>,
        raa: PyReadonlyArray2<f64>,
    ) -> PyResult<(&'py PyArray2<f64>, &'py PyArray2<f64>)> {
        let vza = vza.as_array();
        let sza = sza.as_array();
        let raa = raa.as_array();

        let shape = (vza.shape()[0], vza.shape()[1]);
        let mut ross = Array2::<f64>::zeros(shape);
        let mut li = Array2::<f64>::zeros(shape);

        // Release the GIL while running the rayon-parallel kernel computation
        // so other Python threads can make progress. All buffers are Rust-owned;
        // we only touch Python again to construct the return arrays below.
        py.allow_threads(|| {
            Zip::from(&mut ross)
                .and(&mut li)
                .and(&vza)
                .and(&sza)
                .and(&raa)
                .par_for_each(|r, l, &v, &s, &a| {
                    let (ross_val, li_val) = self.compute_single(v, s, a);
                    *r = ross_val;
                    *l = li_val;
                });
        });

        Ok((ross.into_pyarray(py), li.into_pyarray(py)))
    }
}

impl RossThickLiSparse {
    /// Compute kernel values for a single pixel
    fn compute_single(&self, vza: f64, sza: f64, raa: f64) -> (f64, f64) {
        let cos_sza = sza.cos();
        let cos_vza = vza.cos();
        let sin_sza = sza.sin();
        let sin_vza = vza.sin();
        let raa_abs = raa.abs();
        let cos_raa = raa_abs.cos();
        let sin_raa = raa_abs.sin();

        // Phase angle
        let cos_phase = (cos_sza * cos_vza + sin_sza * sin_vza * cos_raa).clamp(-1.0, 1.0);
        let phase = cos_phase.acos();

        // Ross-Thick kernel
        let ross = self.ross_thick(cos_sza, cos_vza, cos_phase, phase);

        // Li-Sparse kernel
        let li = self.li_sparse(cos_sza, cos_vza, sin_sza, sin_vza, cos_raa, sin_raa);

        (ross, li)
    }

    /// Ross-Thick volumetric scattering kernel
    fn ross_thick(&self, cos_sza: f64, cos_vza: f64, cos_phase: f64, phase: f64) -> f64 {
        let denom = cos_sza + cos_vza;
        if denom.abs() < 1e-10 {
            return 0.0;
        }

        // Main Ross-Thick formula (Lucht 2000 / Wanner 1995), normalized to 0
        // at nadir via the -π/4 constant. A previous -π/2 left a constant -π/4
        // offset that darkened every BRDF-derived reflectance by ~0.785*f_vol.
        ((FRAC_PI_2 - phase) * cos_phase + phase.sin()) / denom - FRAC_PI_4
    }

    /// Li-Sparse geometric scattering kernel
    fn li_sparse(
        &self,
        cos_sza: f64,
        cos_vza: f64,
        sin_sza: f64,
        sin_vza: f64,
        cos_raa: f64,
        sin_raa: f64,
    ) -> f64 {
        // Prime angles (scaled by br for sparse vegetation)
        let tan_sza = sin_sza / cos_sza.max(1e-10);
        let tan_vza = sin_vza / cos_vza.max(1e-10);

        let tan_sza_prime = self.br * tan_sza;
        let tan_vza_prime = self.br * tan_vza;

        let sza_prime = tan_sza_prime.atan();
        let vza_prime = tan_vza_prime.atan();

        let cos_sza_prime = sza_prime.cos();
        let cos_vza_prime = vza_prime.cos();
        let sin_sza_prime = sza_prime.sin();
        let sin_vza_prime = vza_prime.sin();

        // Prime phase angle
        let cos_phase_prime = (cos_sza_prime * cos_vza_prime
            + sin_sza_prime * sin_vza_prime * cos_raa)
            .clamp(-1.0, 1.0);

        // Distance term
        let d2 = tan_sza_prime.powi(2) + tan_vza_prime.powi(2)
            - 2.0 * tan_sza_prime * tan_vza_prime * cos_raa;
        // sec values
        let sec_sza_prime = 1.0 / cos_sza_prime.max(1e-10);
        let sec_vza_prime = 1.0 / cos_vza_prime.max(1e-10);

        // Overlap function
        let cost = self.hb
            * ((d2 + (tan_sza_prime * tan_vza_prime * sin_raa).powi(2)).sqrt()
                / (sec_sza_prime + sec_vza_prime));
        let cost = cost.clamp(-1.0, 1.0);
        let t = cost.acos();
        let overlap = (1.0 / PI) * (t - t.sin() * t.cos()) * (sec_sza_prime + sec_vza_prime);

        // Final Li-Sparse kernel
        overlap - sec_sza_prime - sec_vza_prime
            + 0.5 * (1.0 + cos_phase_prime) * sec_sza_prime * sec_vza_prime
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_nadir_nadir_known_answer() {
        // At vza=sza=0 (sun and sensor both at nadir), the standard MCD43
        // RossThick/LiSparse kernels are normalized to 0 (the BRDF there is the
        // pure isotropic term f_iso). raa is irrelevant at nadir.
        let kernel = RossThickLiSparse::new(2.0, 1.0);
        let (ross, li) = kernel.compute_single(0.0, 0.0, 0.0);
        assert_relative_eq!(ross, 0.0, epsilon = 1e-10);
        assert_relative_eq!(li, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_raa_symmetry() {
        // Both kernels depend only on |raa| (cos(raa) and sin(raa_abs));
        // flipping the sign of raa must leave the kernel values unchanged.
        let kernel = RossThickLiSparse::new(2.0, 1.0);
        let (ross_pos, li_pos) = kernel.compute_single(0.3, 0.5, 0.5);
        let (ross_neg, li_neg) = kernel.compute_single(0.3, 0.5, -0.5);
        assert_relative_eq!(ross_pos, ross_neg, epsilon = 1e-10);
        assert_relative_eq!(li_pos, li_neg, epsilon = 1e-10);
    }

    #[test]
    fn test_outputs_finite_over_geometry_grid() {
        // Guard against NaN/Inf for representative view/solar geometries.
        let kernel = RossThickLiSparse::new(2.0, 1.0);
        for sza_deg in [0.0f64, 15.0, 30.0, 45.0, 60.0, 75.0] {
            for vza_deg in [0.0f64, 15.0, 30.0, 45.0, 60.0] {
                for raa_deg in [0.0f64, 45.0, 90.0, 135.0, 180.0] {
                    let (ross, li) = kernel.compute_single(
                        vza_deg.to_radians(),
                        sza_deg.to_radians(),
                        raa_deg.to_radians(),
                    );
                    assert!(
                        ross.is_finite(),
                        "ross not finite at sza={sza_deg} vza={vza_deg} raa={raa_deg}"
                    );
                    assert!(
                        li.is_finite(),
                        "li not finite at sza={sza_deg} vza={vza_deg} raa={raa_deg}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_br_scaling_affects_li_only() {
        // The br (crown-height) parameter only enters the Li-Sparse kernel.
        // Ross-Thick is independent of br; Li-Sparse must change.
        let k_ref = RossThickLiSparse::new(2.0, 1.0);
        let k_tall = RossThickLiSparse::new(2.0, 2.5);
        let (ross_ref, li_ref) = k_ref.compute_single(0.4, 0.6, 0.5);
        let (ross_tall, li_tall) = k_tall.compute_single(0.4, 0.6, 0.5);
        assert_relative_eq!(ross_ref, ross_tall, epsilon = 1e-12);
        assert!(
            (li_ref - li_tall).abs() > 1e-6,
            "Li kernel should depend on br"
        );
    }
}
