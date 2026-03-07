"""Spectral mapping between differing source and target sensor band sets."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from sklearn.neighbors import NearestNeighbors

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.core.types import SensorBand


_VISIBLE_UPPER_NM = 700.0
_DEFAULT_K_NEIGHBORS = 5
_DISTANCE_EPS = 1.0e-6
_UNCERTAINTY_FLOOR = 0.005


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


@dataclass(frozen=True)
class HyperspectralLibrary:
    """Hyperspectral library sampled on a common wavelength axis."""

    wavelengths_nm: np.ndarray
    spectra: np.ndarray
    sample_ids: tuple[str, ...]
    source_id: str = "siac-default-spectral-library"
    source_version: str = "1"

    def __post_init__(self) -> None:
        wavelengths = np.asarray(self.wavelengths_nm, dtype=np.float32)
        spectra = np.asarray(self.spectra, dtype=np.float32)
        if wavelengths.ndim != 1 or wavelengths.size < 2:
            raise ValueError("wavelengths_nm must be a 1-D array with at least two samples")
        if spectra.ndim != 2 or spectra.shape[1] != wavelengths.size:
            raise ValueError("spectra must have shape (n_samples, n_wavelengths)")
        if spectra.shape[0] < 1:
            raise ValueError("spectra must contain at least one sample")
        if len(self.sample_ids) != spectra.shape[0]:
            raise ValueError("sample_ids must match the number of spectra")
        if not np.all(np.isfinite(wavelengths)) or np.any(np.diff(wavelengths) <= 0.0):
            raise ValueError("wavelengths_nm must be finite and strictly increasing")
        if not np.all(np.isfinite(spectra)):
            raise ValueError("spectra must be finite")
        if np.any((spectra < 0.0) | (spectra > 1.5)):
            raise ValueError("spectra must be bounded to a physically plausible reflectance range")
        object.__setattr__(self, "wavelengths_nm", wavelengths)
        object.__setattr__(self, "spectra", spectra)


def _sigmoid(wavelengths_nm: np.ndarray, center_nm: float, width_nm: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(wavelengths_nm - center_nm) / width_nm))


def _gaussian(wavelengths_nm: np.ndarray, center_nm: float, sigma_nm: float) -> np.ndarray:
    return np.exp(-0.5 * np.square((wavelengths_nm - center_nm) / sigma_nm))


@lru_cache(maxsize=1)
def build_default_spectral_library() -> HyperspectralLibrary:
    """Build the default shipped spectral library used by Appendix-D mapping."""
    wavelengths = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)

    def _vegetation() -> np.ndarray:
        vis = (
            0.035
            + 0.045 * _gaussian(wavelengths, 550.0, 35.0)
            - 0.030 * _gaussian(wavelengths, 450.0, 22.0)
            - 0.040 * _gaussian(wavelengths, 675.0, 22.0)
        )
        red_edge = 0.42 * _sigmoid(wavelengths, 715.0, 18.0)
        water_abs = (
            1.0
            - 0.10 * _gaussian(wavelengths, 970.0, 25.0)
            - 0.08 * _gaussian(wavelengths, 1200.0, 30.0)
            - 0.38 * _gaussian(wavelengths, 1450.0, 45.0)
            - 0.48 * _gaussian(wavelengths, 1940.0, 70.0)
            - 0.12 * _gaussian(wavelengths, 2200.0, 50.0)
        )
        swir = 0.04 * _gaussian(wavelengths, 1650.0, 180.0)
        return np.clip((vis + red_edge + swir) * water_abs, 0.0, 0.95)

    def _dry_vegetation() -> np.ndarray:
        vis = (
            0.08
            + 0.025 * _gaussian(wavelengths, 560.0, 45.0)
            - 0.020 * _gaussian(wavelengths, 670.0, 28.0)
        )
        red_edge = 0.24 * _sigmoid(wavelengths, 720.0, 22.0)
        swir = 0.10 + 0.06 * _gaussian(wavelengths, 1700.0, 220.0)
        water_abs = (
            1.0
            - 0.08 * _gaussian(wavelengths, 970.0, 25.0)
            - 0.04 * _gaussian(wavelengths, 1200.0, 30.0)
            - 0.16 * _gaussian(wavelengths, 1450.0, 40.0)
            - 0.22 * _gaussian(wavelengths, 1940.0, 60.0)
        )
        return np.clip((vis + red_edge + swir) * water_abs, 0.0, 0.9)

    def _soil(brightness: float, slope: float) -> np.ndarray:
        baseline = brightness + slope * (wavelengths - 400.0)
        features = (
            1.0
            - 0.03 * _gaussian(wavelengths, 1400.0, 50.0)
            - 0.04 * _gaussian(wavelengths, 1900.0, 80.0)
            - 0.02 * _gaussian(wavelengths, 2200.0, 55.0)
        )
        return np.clip(baseline * features, 0.0, 0.8)

    def _water() -> np.ndarray:
        baseline = 0.025 * np.exp(-(wavelengths - 400.0) / 280.0)
        absorption = (
            1.0
            - 0.60 * _gaussian(wavelengths, 740.0, 45.0)
            - 0.85 * _gaussian(wavelengths, 970.0, 55.0)
            - 0.95 * _gaussian(wavelengths, 1200.0, 70.0)
            - 0.98 * _gaussian(wavelengths, 1450.0, 90.0)
            - 0.99 * _gaussian(wavelengths, 1940.0, 110.0)
        )
        return np.clip(baseline * absorption, 0.0, 0.08)

    def _snow() -> np.ndarray:
        baseline = 0.85 - 0.10 * (wavelengths - 400.0) / 2100.0
        absorption = (
            1.0
            - 0.08 * _gaussian(wavelengths, 1030.0, 35.0)
            - 0.15 * _gaussian(wavelengths, 1250.0, 40.0)
            - 0.55 * _gaussian(wavelengths, 1500.0, 70.0)
            - 0.75 * _gaussian(wavelengths, 2000.0, 110.0)
        )
        return np.clip(baseline * absorption, 0.0, 0.95)

    endmembers = np.stack(
        [
            _vegetation(),
            _dry_vegetation(),
            _soil(0.05, 8.0e-5),
            _soil(0.12, 1.4e-4),
            _water(),
            _snow(),
        ],
        axis=0,
    ).astype(np.float32)

    rng = np.random.default_rng(42)
    spectra: list[np.ndarray] = []
    sample_ids: list[str] = []
    for idx in range(384):
        weights = rng.dirichlet(np.array([2.2, 1.8, 1.3, 1.2, 0.35, 0.25], dtype=np.float32))
        spectrum = np.sum(weights[:, np.newaxis] * endmembers, axis=0)
        moisture = float(rng.uniform(0.0, 0.35))
        moisture_mask = (
            1.0
            - moisture * 0.15 * _gaussian(wavelengths, 970.0, 24.0)
            - moisture * 0.20 * _gaussian(wavelengths, 1200.0, 28.0)
            - moisture * 0.35 * _gaussian(wavelengths, 1450.0, 42.0)
            - moisture * 0.45 * _gaussian(wavelengths, 1940.0, 65.0)
            - moisture * 0.12 * _gaussian(wavelengths, 2200.0, 52.0)
        )
        red_edge_shift = float(rng.uniform(-18.0, 18.0))
        red_edge_delta = float(rng.uniform(-0.04, 0.05))
        slope = float(rng.uniform(-0.05, 0.08))
        jitter = rng.normal(0.0, 0.002, size=wavelengths.size).astype(np.float32)
        spectrum = spectrum * moisture_mask
        spectrum += red_edge_delta * _sigmoid(wavelengths, 710.0 + red_edge_shift, 17.0)
        spectrum *= 1.0 + slope * (wavelengths - wavelengths.mean()) / 1000.0
        spectra.append(np.clip(spectrum + jitter, 0.0, 0.98).astype(np.float32))
        sample_ids.append(f"sample_{idx:03d}")

    return HyperspectralLibrary(
        wavelengths_nm=wavelengths,
        spectra=np.stack(spectra, axis=0),
        sample_ids=tuple(sample_ids),
    )


def _normalized_band_response(band: SensorBand, wavelengths_nm: np.ndarray) -> np.ndarray:
    response = np.asarray(band.effective_response(wavelengths_nm), dtype=np.float32)
    area = _trapezoid(response, wavelengths_nm)
    if not np.isfinite(area) or area <= 0.0:
        raise ValueError(f"Band {band.name!r} has zero support on the requested wavelength grid")
    return response / area


def _project_library_to_bands(
    library: HyperspectralLibrary,
    bands: Sequence[SensorBand],
) -> np.ndarray:
    matrix = np.stack(
        [_normalized_band_response(band, library.wavelengths_nm) for band in bands],
        axis=0,
    )
    return np.asarray(library.spectra @ matrix.T, dtype=np.float32)


def _classify_band_region(band: SensorBand) -> str:
    return "visible" if band.center_wavelength < _VISIBLE_UPPER_NM else "infrared"


def needs_spectral_mapping(
    source_bands: Sequence[SensorBand],
    target_bands: Sequence[SensorBand],
    *,
    response_tolerance: float = 5.0e-3,
) -> bool:
    """Return True when source and target band sets differ materially."""
    if len(source_bands) != len(target_bands):
        return True

    for source_band, target_band in zip(source_bands, target_bands, strict=True):
        if source_band.name != target_band.name:
            return True
        common_min = max(
            source_band.center_wavelength - 4.0 * source_band.bandwidth,
            target_band.center_wavelength - 4.0 * target_band.bandwidth,
            350.0,
        )
        common_max = min(
            source_band.center_wavelength + 4.0 * source_band.bandwidth,
            target_band.center_wavelength + 4.0 * target_band.bandwidth,
            2500.0,
        )
        if common_max <= common_min:
            return True
        grid = np.arange(common_min, common_max + 1.0, 1.0, dtype=np.float32)
        source_resp = _normalized_band_response(source_band, grid)
        target_resp = _normalized_band_response(target_band, grid)
        if float(np.max(np.abs(source_resp - target_resp))) > response_tolerance:
            return True
    return False


def convolve_hyperspectral_reflectance(
    reflectance: xr.DataArray,
    wavelengths_nm: np.ndarray,
    target_bands: Sequence[SensorBand],
) -> xr.DataArray:
    """Project a hyperspectral reflectance cube onto multispectral target bands."""
    if "wavelength" not in reflectance.dims:
        raise ValueError("reflectance must have a 'wavelength' dimension")
    source = reflectance.transpose(*[dim for dim in reflectance.dims if dim != "wavelength"], "wavelength")
    if source.sizes["wavelength"] != len(wavelengths_nm):
        raise ValueError("wavelength dimension does not match wavelengths_nm")

    matrix = np.stack(
        [_normalized_band_response(band, np.asarray(wavelengths_nm, dtype=np.float32)) for band in target_bands],
        axis=0,
    ).astype(np.float32)
    values = np.asarray(source.values, dtype=np.float32)
    flat = values.reshape(-1, values.shape[-1])
    projected = flat @ matrix.T
    out_shape = values.shape[:-1] + (len(target_bands),)
    projected = projected.reshape(out_shape)

    other_dims = [dim for dim in source.dims if dim != "wavelength"]
    coords = {dim: source.coords[dim] for dim in other_dims if dim in source.coords}
    coords["band"] = [band.name for band in target_bands]
    return xr.DataArray(
        np.moveaxis(projected, -1, 0),
        dims=["band", *other_dims],
        coords=coords,
    )


class SpectralMapper:
    """Appendix-D-style mapper from one multispectral basis to another."""

    def __init__(
        self,
        source_bands: Sequence[SensorBand],
        target_bands: Sequence[SensorBand],
        *,
        spectral_library: HyperspectralLibrary | None = None,
        k_neighbors: int = _DEFAULT_K_NEIGHBORS,
    ) -> None:
        if k_neighbors < 1:
            raise ValueError("k_neighbors must be >= 1")
        self.source_bands = tuple(source_bands)
        self.target_bands = tuple(target_bands)
        self.spectral_library = spectral_library or build_default_spectral_library()
        self.k_neighbors = int(k_neighbors)
        self._identity = not needs_spectral_mapping(self.source_bands, self.target_bands)

        self._source_projection = _project_library_to_bands(self.spectral_library, self.source_bands)
        self._target_projection = _project_library_to_bands(self.spectral_library, self.target_bands)
        self._source_name_to_index = {band.name: idx for idx, band in enumerate(self.source_bands)}
        self._target_band_names = [band.name for band in self.target_bands]
        self._visible_source_idx = np.array(
            [idx for idx, band in enumerate(self.source_bands) if _classify_band_region(band) == "visible"],
            dtype=np.int32,
        )
        self._infrared_source_idx = np.array(
            [idx for idx, band in enumerate(self.source_bands) if _classify_band_region(band) == "infrared"],
            dtype=np.int32,
        )
        self._visible_target_idx = np.array(
            [idx for idx, band in enumerate(self.target_bands) if _classify_band_region(band) == "visible"],
            dtype=np.int32,
        )
        self._infrared_target_idx = np.array(
            [idx for idx, band in enumerate(self.target_bands) if _classify_band_region(band) == "infrared"],
            dtype=np.int32,
        )
        self._nn_visible = self._make_nn(self._visible_source_idx)
        self._nn_infrared = self._make_nn(self._infrared_source_idx)

    def _make_nn(self, source_idx: np.ndarray) -> NearestNeighbors | None:
        if source_idx.size == 0:
            return None
        features = self._source_projection[:, source_idx]
        if features.shape[0] == 0:
            return None
        model = NearestNeighbors(
            n_neighbors=min(self.k_neighbors, features.shape[0]),
            algorithm="auto",
            metric="euclidean",
        )
        model.fit(features)
        return model

    def map(
        self,
        source_reflectance: xr.DataArray,
        *,
        source_uncertainty: xr.DataArray | None = None,
    ) -> tuple[xr.DataArray, xr.DataArray]:
        """Map source-basis multispectral reflectance to target bands."""
        source_data = self._align_source_data(source_reflectance)
        original_dims = tuple(source_data.dims)
        if source_uncertainty is not None:
            source_unc = self._align_source_data(source_uncertainty)
        else:
            source_unc = None

        if self._identity:
            identity = source_data.assign_coords(band=self._target_band_names)
            if source_unc is not None:
                unc = source_unc.assign_coords(band=self._target_band_names)
            else:
                unc = xr.zeros_like(identity, dtype=np.float32) + _UNCERTAINTY_FLOOR
            return (
                identity.transpose(*original_dims).astype(np.float32),
                unc.transpose(*original_dims).astype(np.float32),
            )

        transpose_dims = ["band", *[dim for dim in source_data.dims if dim != "band"]]
        source_data = source_data.transpose(*transpose_dims)
        values = np.asarray(source_data.values, dtype=np.float32)
        spatial_dims = transpose_dims[1:]
        flat = values.reshape(values.shape[0], -1).T

        if source_unc is not None:
            unc_values = np.asarray(source_unc.transpose(*transpose_dims).values, dtype=np.float32)
            unc_flat = unc_values.reshape(unc_values.shape[0], -1).T
        else:
            unc_flat = None

        target_flat = np.full((flat.shape[0], len(self.target_bands)), np.nan, dtype=np.float32)
        target_unc_flat = np.full_like(target_flat, np.nan)

        self._map_region(
            flat,
            target_flat,
            target_unc_flat,
            source_idx=self._visible_source_idx,
            target_idx=self._visible_target_idx,
            nn=self._nn_visible,
            unc_flat=unc_flat,
        )
        self._map_region(
            flat,
            target_flat,
            target_unc_flat,
            source_idx=self._infrared_source_idx,
            target_idx=self._infrared_target_idx,
            nn=self._nn_infrared,
            unc_flat=unc_flat,
        )

        shape = (len(self.target_bands), *[source_data.sizes[dim] for dim in spatial_dims])
        reflectance = target_flat.T.reshape(shape)
        uncertainty = target_unc_flat.T.reshape(shape)
        coords = {"band": self._target_band_names}
        coords.update({dim: source_data.coords[dim] for dim in spatial_dims if dim in source_data.coords})
        reflectance_da = xr.DataArray(reflectance, dims=transpose_dims, coords=coords).astype(np.float32)
        uncertainty_da = xr.DataArray(uncertainty, dims=transpose_dims, coords=coords).astype(np.float32)
        return reflectance_da.transpose(*original_dims), uncertainty_da.transpose(*original_dims)

    def _align_source_data(self, data: xr.DataArray) -> xr.DataArray:
        if "band" not in data.dims:
            raise ValueError("source_reflectance must have a 'band' dimension")
        missing = [band.name for band in self.source_bands if band.name not in set(data.coords["band"].values.tolist())]
        if missing:
            raise KeyError(f"source_reflectance is missing source bands: {missing}")
        return data.sel(band=[band.name for band in self.source_bands])

    def _map_region(
        self,
        flat: np.ndarray,
        target_flat: np.ndarray,
        target_unc_flat: np.ndarray,
        *,
        source_idx: np.ndarray,
        target_idx: np.ndarray,
        nn: NearestNeighbors | None,
        unc_flat: np.ndarray | None,
    ) -> None:
        if target_idx.size == 0:
            return
        if source_idx.size == 0 or nn is None:
            return

        query = flat[:, source_idx]
        valid = np.all(np.isfinite(query), axis=1)
        if not np.any(valid):
            return

        distances, neighbors = nn.kneighbors(query[valid], return_distance=True)
        neighbor_source = self._source_projection[neighbors][:, :, source_idx]
        neighbor_target = self._target_projection[neighbors][:, :, target_idx]

        weights = np.zeros_like(distances, dtype=np.float32)
        zero_mask = distances <= _DISTANCE_EPS
        if np.any(zero_mask):
            zero_rows = np.any(zero_mask, axis=1)
            weights[zero_rows] = zero_mask[zero_rows].astype(np.float32)
        nonzero_rows = ~np.any(zero_mask, axis=1)
        if np.any(nonzero_rows):
            inv = 1.0 / np.maximum(distances[nonzero_rows], _DISTANCE_EPS)
            weights[nonzero_rows] = inv
        weights_sum = np.sum(weights, axis=1, keepdims=True)
        weights = np.divide(weights, weights_sum, out=np.zeros_like(weights), where=weights_sum > 0.0)

        estimate = np.einsum("nk,nkb->nb", weights, neighbor_target, optimize=True).astype(np.float32)
        recon_source = np.einsum("nk,nks->ns", weights, neighbor_source, optimize=True).astype(np.float32)
        source_rmse = np.sqrt(
            np.mean(np.square(recon_source - query[valid]), axis=1, dtype=np.float32)
        ).astype(np.float32)
        target_spread = np.sqrt(
            np.einsum(
                "nk,nkb->nb",
                weights,
                np.square(neighbor_target - estimate[:, np.newaxis, :], dtype=np.float32),
                optimize=True,
            )
        ).astype(np.float32)
        if unc_flat is not None:
            input_unc = np.sqrt(np.mean(np.square(unc_flat[valid][:, source_idx]), axis=1, dtype=np.float32))
        else:
            input_unc = np.zeros(source_rmse.shape, dtype=np.float32)
        total_unc = np.sqrt(
            np.maximum(target_spread, _UNCERTAINTY_FLOOR) ** 2
            + source_rmse[:, np.newaxis] ** 2
            + input_unc[:, np.newaxis] ** 2
        ).astype(np.float32)

        target_flat[np.flatnonzero(valid)[:, np.newaxis], target_idx[np.newaxis, :]] = estimate
        target_unc_flat[np.flatnonzero(valid)[:, np.newaxis], target_idx[np.newaxis, :]] = total_unc


def map_multispectral_reflectance(
    source_reflectance: xr.DataArray,
    *,
    source_bands: Sequence[SensorBand],
    target_bands: Sequence[SensorBand],
    source_uncertainty: xr.DataArray | None = None,
    spectral_library: HyperspectralLibrary | None = None,
    k_neighbors: int = _DEFAULT_K_NEIGHBORS,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Convenience wrapper around :class:`SpectralMapper`."""
    mapper = SpectralMapper(
        source_bands,
        target_bands,
        spectral_library=spectral_library,
        k_neighbors=k_neighbors,
    )
    return mapper.map(source_reflectance, source_uncertainty=source_uncertainty)


__all__ = [
    "HyperspectralLibrary",
    "SpectralMapper",
    "build_default_spectral_library",
    "convolve_hyperspectral_reflectance",
    "map_multispectral_reflectance",
    "needs_spectral_mapping",
]
