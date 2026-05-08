"""Route-A temporal BRDF smoothing using the Rust Whittaker smoother."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt

from siac._rust_compat import whittaker_smooth_cube
from siac.algorithms.brdf.kernels import BRDFKernels, compute_reflectance
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.algorithms.surface.spectral_mapping import map_multispectral_reflectance
from siac.geo._spatial import copy_spatial_metadata_like
from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import datetime

    from siac.algorithms.surface.spectral_mapping import SpectralMappingConfig
    from siac.domain import SensorBand

logger = logging.getLogger(__name__)
Float32Array: TypeAlias = npt.NDArray[np.float32]


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
        psf_params: tuple[float, float] | None = None,
        *,
        source_bands: Sequence[SensorBand] | None = None,
        target_bands: Sequence[SensorBand] | None = None,
        spectral_library: SpectralMappingConfig | None = None,
        spectral_k_neighbors: int = 5,
        **kwargs: Any,
    ) -> SurfacePrior:
        obs_time = self._pop_obs_time(kwargs)
        if "time" not in brdf_weights.f0.dims:
            raise ValueError(
                "Whittaker BRDF derivation requires BRDF weights with a 'time' dimension"
            )
        sigma_x, sigma_y = self._resolve_psf_sigmas(psf_params)

        ref = brdf_weights.f0.isel(time=0, drop=True)
        k_vol, k_geo = self._kernels.compute(geometry.vza, geometry.sza, geometry.raa)
        k_vol = self._require_data_array(k_vol, name="k_vol")
        k_geo = self._require_data_array(k_geo, name="k_geo")
        k_vol, k_geo = self._align_kernels_to_brdf_grid(k_vol, k_geo, ref)

        reflectance = self._require_data_array(
            compute_reflectance(
                brdf_weights.f0,
                brdf_weights.f1,
                brdf_weights.f2,
                k_vol,
                k_geo,
            ),
            name="reflectance",
        ).transpose("time", "band", "y", "x")
        reflectance_unc = self._require_data_array(
            brdf_weights.compute_reflectance_uncertainty(k_vol, k_geo),
            name="reflectance_unc",
        ).transpose("time", "band", "y", "x")

        reflectance_values = np.asarray(reflectance.values, dtype=np.float32)
        unc_values = np.asarray(reflectance_unc.values, dtype=np.float32)
        logger.info(
            "Route-A Whittaker BRDF prior: %d days, %d band%s.",
            reflectance_values.shape[0],
            reflectance_values.shape[1],
            "" if reflectance_values.shape[1] == 1 else "s",
        )
        weights = self._normalized_temporal_weights(reflectance_values, unc_values)

        smoothed = whittaker_smooth_cube(
            np.ascontiguousarray(reflectance_values, dtype=np.float32),
            np.ascontiguousarray(weights, dtype=np.float32),
            self.temporal_lambda,
        )
        smoothed = np.asarray(smoothed, dtype=np.float32)

        target_index = self._target_index(brdf_weights.f0.coords["time"].values, obs_time)
        has_data = np.any(weights > 0.0, axis=0)
        prior = smoothed[target_index]

        prior_unc = self._prior_uncertainty(
            reflectance_values,
            unc_values,
            smoothed,
            weights,
            target_index,
        )

        # Fallback values for no-data pixels — see ``siac.constants`` for
        # rationale. Previously hard-coded ``0.20`` / ``0.08`` (REVIEW.md §3.5).
        from siac.constants import DEFAULT_NO_DATA_BOA, DEFAULT_NO_DATA_BOA_UNC

        prior = np.where(has_data, prior, DEFAULT_NO_DATA_BOA).astype(np.float32)
        prior_unc = np.where(has_data, prior_unc, DEFAULT_NO_DATA_BOA_UNC).astype(np.float32)

        coords = {
            "band": brdf_weights.f0.coords["band"],
            "y": brdf_weights.f0.coords["y"],
            "x": brdf_weights.f0.coords["x"],
        }
        boa = copy_spatial_metadata_like(
            xr.DataArray(prior, dims=["band", "y", "x"], coords=coords),
            ref,
        )
        boa_unc = copy_spatial_metadata_like(
            xr.DataArray(prior_unc, dims=["band", "y", "x"], coords=coords),
            ref,
        )

        if source_bands and target_bands:
            boa, boa_unc = map_multispectral_reflectance(
                boa,
                source_bands=source_bands,
                target_bands=target_bands,
                source_uncertainty=boa_unc,
                spectral_library=spectral_library,
                k_neighbors=spectral_k_neighbors,
            )

        if self.apply_psf and sigma_x > 0 and sigma_y > 0:
            boa = copy_spatial_metadata_like(self._apply_psf(boa, sigma_x, sigma_y), boa)
            boa_unc = copy_spatial_metadata_like(
                self._apply_psf(boa_unc, sigma_x, sigma_y), boa_unc
            )

        mask = copy_spatial_metadata_like(
            xr.DataArray(
                np.all(np.isfinite(boa.values), axis=0)
                & np.all(np.isfinite(boa_unc.values), axis=0),
                dims=["y", "x"],
                coords={"y": boa.coords["y"], "x": boa.coords["x"]},
            ),
            boa,
        )
        return SurfacePrior(
            boa=boa,
            boa_unc=boa_unc,
            kernels=None,
            mask=mask,
        )

    @staticmethod
    def _pop_obs_time(kwargs: dict[str, Any]) -> datetime | None:
        obs_time_value = kwargs.pop("obs_time", None)
        if kwargs:
            unexpected = ", ".join(sorted(kwargs))
            raise TypeError(f"Unexpected keyword argument(s): {unexpected}")
        return cast("datetime | None", obs_time_value)

    def _resolve_psf_sigmas(
        self,
        psf_params: tuple[float, float] | None,
    ) -> tuple[float, float]:
        if psf_params is None:
            return self.psf_sigma_x, self.psf_sigma_y
        return psf_params

    @staticmethod
    def _normalized_temporal_weights(
        reflectance_values: Float32Array,
        unc_values: Float32Array,
    ) -> Float32Array:
        weights = np.where(
            np.isfinite(reflectance_values) & np.isfinite(unc_values) & (unc_values > 0.0),
            1.0 / np.maximum(unc_values**2, 1.0e-6),
            0.0,
        ).astype(np.float32)
        max_weight = np.max(weights, axis=0, keepdims=True)
        normalized_weights = np.zeros_like(weights, dtype=np.float32)
        np.divide(weights, max_weight, out=normalized_weights, where=max_weight > 0.0)
        return normalized_weights

    @staticmethod
    def _target_uncertainty(
        unc_values: Float32Array,
        target_index: int,
    ) -> Float32Array:
        target_unc = np.asarray(unc_values[target_index], dtype=np.float32)
        fallback_unc = np.full(target_unc.shape, np.nan, dtype=np.float32)
        finite_any = np.any(np.isfinite(unc_values), axis=0)
        if np.any(finite_any):
            fallback_unc[finite_any] = np.asarray(
                np.nanmin(unc_values[:, finite_any], axis=0),
                dtype=np.float32,
            )
        return np.asarray(
            np.where(np.isfinite(target_unc), target_unc, fallback_unc), dtype=np.float32
        )

    @classmethod
    def _prior_uncertainty(
        cls,
        reflectance_values: Float32Array,
        unc_values: Float32Array,
        smoothed: Float32Array,
        weights: Float32Array,
        target_index: int,
    ) -> Float32Array:
        target_unc = cls._target_uncertainty(unc_values, target_index)
        residual = np.asarray(
            np.where(weights > 0.0, reflectance_values - smoothed, np.nan),
            dtype=np.float32,
        )
        weighted_var = np.asarray(
            np.sum(np.where(np.isfinite(residual), weights * residual**2, 0.0), axis=0),
            dtype=np.float32,
        )
        weight_sum = np.asarray(np.sum(weights, axis=0), dtype=np.float32)
        residual_ratio = np.full(weighted_var.shape, np.nan, dtype=np.float32)
        np.divide(weighted_var, weight_sum, out=residual_ratio, where=weight_sum > 0.0)
        residual_std = np.sqrt(residual_ratio).astype(np.float32)
        return np.asarray(
            np.sqrt(np.maximum(target_unc, 0.0) ** 2 + np.nan_to_num(residual_std, nan=0.0) ** 2),
            dtype=np.float32,
        )

    @staticmethod
    def _require_data_array(
        value: xr.DataArray | np.ndarray[Any, Any], *, name: str
    ) -> xr.DataArray:
        if not isinstance(value, xr.DataArray):
            raise TypeError(f"{name} must be an xarray.DataArray")
        return value

    @staticmethod
    def _target_index(time_values: np.ndarray, obs_time: datetime | None) -> int:
        if time_values.size == 0:
            raise ValueError("Temporal BRDF stack has no time samples")
        if obs_time is None:
            return time_values.size // 2
        obs_np = np.datetime64(obs_time.date())
        deltas = np.abs(time_values.astype("datetime64[D]") - obs_np)
        return int(np.argmin(deltas))
