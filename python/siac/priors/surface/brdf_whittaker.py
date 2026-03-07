"""Route-A temporal BRDF smoothing using the Rust Whittaker smoother."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac._rust import whittaker_smooth_cube
from siac.core.types import BRDFKernelWeights, GeometryAngles, SurfacePrior
from siac.priors.brdf.kernels import BRDFKernels, compute_reflectance
from siac.priors.surface.kernel_model import KernelModelDeriver

if TYPE_CHECKING:
    from datetime import datetime

logger = logging.getLogger(__name__)


class BRDFWhittakerDeriver(KernelModelDeriver):
    """Derive a sensing-date surface prior from a temporal BRDF stack."""

    def __init__(
        self,
        *,
        temporal_lambda: float = 10.0,
        psf_sigma_x: float = 29.75,
        psf_sigma_y: float = 39.0,
        apply_psf: bool = True,
    ) -> None:
        super().__init__(
            psf_sigma_x=psf_sigma_x,
            psf_sigma_y=psf_sigma_y,
            apply_psf=apply_psf,
        )
        if temporal_lambda <= 0:
            raise ValueError("temporal_lambda must be > 0")
        self.temporal_lambda = float(temporal_lambda)
        self._kernels = BRDFKernels(hb=2.0, br=1.0)

    def compute_surface_prior(
        self,
        brdf_weights: BRDFKernelWeights,
        geometry: GeometryAngles,
        *,
        obs_time: datetime | None = None,
    ) -> SurfacePrior:
        if "time" not in brdf_weights.f0.dims:
            raise ValueError("Whittaker BRDF derivation requires BRDF weights with a 'time' dimension")

        ref = brdf_weights.f0.isel(time=0, drop=True)
        k_vol, k_geo = self._kernels.compute(geometry.vza, geometry.sza, geometry.raa)
        k_vol, k_geo = self._align_kernels_to_brdf_grid(k_vol, k_geo, ref)

        reflectance = compute_reflectance(
            brdf_weights.f0,
            brdf_weights.f1,
            brdf_weights.f2,
            k_vol,
            k_geo,
        ).transpose("time", "band", "y", "x")
        reflectance_unc = np.sqrt(
            brdf_weights.f0_unc**2
            + (k_vol * brdf_weights.f1_unc) ** 2
            + (k_geo * brdf_weights.f2_unc) ** 2
        ).transpose("time", "band", "y", "x")

        reflectance_values = np.asarray(reflectance.values, dtype=np.float32)
        unc_values = np.asarray(reflectance_unc.values, dtype=np.float32)
        logger.info(
            "Route-A Whittaker BRDF prior: %d days, %d band%s.",
            reflectance_values.shape[0],
            reflectance_values.shape[1],
            "" if reflectance_values.shape[1] == 1 else "s",
        )
        weights = np.where(
            np.isfinite(reflectance_values) & np.isfinite(unc_values) & (unc_values > 0.0),
            1.0 / np.maximum(unc_values**2, 1.0e-6),
            0.0,
        ).astype(np.float32)
        max_weight = np.max(weights, axis=0, keepdims=True)
        normalized_weights = np.zeros_like(weights, dtype=np.float32)
        np.divide(weights, max_weight, out=normalized_weights, where=max_weight > 0.0)
        weights = normalized_weights

        smoothed = whittaker_smooth_cube(
            np.ascontiguousarray(reflectance_values, dtype=np.float32),
            np.ascontiguousarray(weights, dtype=np.float32),
            self.temporal_lambda,
        )
        smoothed = np.asarray(smoothed, dtype=np.float32)

        target_index = self._target_index(brdf_weights.f0.coords["time"].values, obs_time)
        has_data = np.any(weights > 0.0, axis=0)
        prior = smoothed[target_index]

        target_unc = unc_values[target_index]
        fallback_unc = np.full(target_unc.shape, np.nan, dtype=np.float32)
        finite_any = np.any(np.isfinite(unc_values), axis=0)
        if np.any(finite_any):
            fallback_unc[finite_any] = np.nanmin(unc_values[:, finite_any], axis=0)
        target_unc = np.where(np.isfinite(target_unc), target_unc, fallback_unc)

        residual = np.where(weights > 0.0, reflectance_values - smoothed, np.nan)
        weighted_var = np.sum(np.where(np.isfinite(residual), weights * residual**2, 0.0), axis=0)
        weight_sum = np.sum(weights, axis=0)
        residual_ratio = np.full(weighted_var.shape, np.nan, dtype=np.float32)
        np.divide(weighted_var, weight_sum, out=residual_ratio, where=weight_sum > 0.0)
        residual_std = np.sqrt(residual_ratio).astype(np.float32)
        prior_unc = np.sqrt(np.maximum(target_unc, 0.0) ** 2 + np.nan_to_num(residual_std, nan=0.0) ** 2)

        prior = np.where(has_data, prior, 0.20).astype(np.float32)
        prior_unc = np.where(has_data, prior_unc, 0.08).astype(np.float32)

        coords = {
            "band": brdf_weights.f0.coords["band"],
            "y": brdf_weights.f0.coords["y"],
            "x": brdf_weights.f0.coords["x"],
        }
        boa = xr.DataArray(prior, dims=["band", "y", "x"], coords=coords)
        boa_unc = xr.DataArray(prior_unc, dims=["band", "y", "x"], coords=coords)

        if self.apply_psf and self.psf_sigma_x > 0 and self.psf_sigma_y > 0:
            boa = self._apply_psf(boa, self.psf_sigma_x, self.psf_sigma_y)
            boa_unc = self._apply_psf(boa_unc, self.psf_sigma_x, self.psf_sigma_y)

        mask = xr.DataArray(
            np.all(np.isfinite(boa.values), axis=0) & np.all(np.isfinite(boa_unc.values), axis=0),
            dims=["y", "x"],
            coords={"y": boa.coords["y"], "x": boa.coords["x"]},
        )
        return SurfacePrior(
            boa=boa,
            boa_unc=boa_unc,
            kernels=None,
            mask=mask,
        )

    @staticmethod
    def _target_index(time_values: np.ndarray, obs_time: datetime | None) -> int:
        if time_values.size == 0:
            raise ValueError("Temporal BRDF stack has no time samples")
        if obs_time is None:
            return time_values.size // 2
        obs_np = np.datetime64(obs_time.date())
        deltas = np.abs(time_values.astype("datetime64[D]") - obs_np)
        return int(np.argmin(deltas))
