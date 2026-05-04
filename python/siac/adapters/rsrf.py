"""Thin adapter layer over the external RSRF package."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rsrf
import rsrf.registry as rsrf_registry
from rsrf.convolve import convolution_weights as _convolution_weights

from siac.catalog import get_sensor_config
from siac.domain import SensorBand, SensorConfig
from siac.domain.spectral import RelativeSpectralResponse


def _normalize_rsrf_root(rsrf_root: Path | str | None) -> Path | None:
    if rsrf_root is None:
        return _default_rsrf_root()
    return Path(rsrf_root).expanduser().resolve()


def _default_rsrf_root() -> Path | None:
    root = rsrf_registry.discover_repo_root(None)
    if _has_runtime_rsrf_assets(root):
        return root
    runtime_root = rsrf_registry._runtime_release_root()
    if _has_runtime_rsrf_assets(runtime_root):
        return runtime_root
    return root


def _has_runtime_rsrf_assets(root: Path) -> bool:
    data_root = root / "data"
    return (
        (data_root / "registry" / "sensors.parquet").exists()
        and (data_root / "registry" / "bands.parquet").exists()
        and (data_root / "canonical").is_dir()
    )


def _curve_from_response_definition(
    response_definition: object,
    *,
    band_id: str,
    source_variant: str | None = None,
) -> object:
    resolved = rsrf.coerce_response_definition(
        response_definition,
        band_id=band_id,
        source_variant=source_variant or "custom",
    )
    if hasattr(resolved, "wavelength_nm") and hasattr(resolved, "response"):
        return resolved
    return rsrf.realize_curve(resolved, source_variant=source_variant)


def band_response_definition_input(band: SensorBand) -> object:
    """Return an rsrf-compatible response-definition input for a SIAC band."""
    band_id = str(band.rsrf_band_id or band.name)
    if band.rsrf_wavelengths_nm is not None and band.rsrf_response is not None:
        return {
            "kind": "sampled",
            "band_id": band_id,
            "wavelength_nm": np.asarray(band.rsrf_wavelengths_nm, dtype=float).tolist(),
            "response": np.asarray(band.rsrf_response, dtype=float).tolist(),
        }
    return {
        "kind": "gaussian",
        "band_id": band_id,
        "band_name": str(band.name),
        "band_index": int(band.band_index),
        "center_wavelength_nm": float(band.center_wavelength),
        "fwhm_nm": float(band.bandwidth),
    }


def coerce_band_rsrf(
    band: SensorBand,
    *,
    sensor_id: str,
    satellite_id: str,
) -> RelativeSpectralResponse:
    """Realize any SIAC band through rsrf and return a tabulated response."""
    curve = _curve_from_response_definition(
        band_response_definition_input(band),
        band_id=str(band.rsrf_band_id or band.name),
        source_variant=band.rsrf_representation_variant or None,
    )
    source_variant = getattr(curve, "source_variant", None)
    return RelativeSpectralResponse.from_tabulated(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        band_name=band.name,
        wavelengths_nm=np.asarray(curve.wavelength_nm, dtype=np.float32),
        response=np.asarray(curve.response, dtype=np.float32),
        source_id="RSRF",
        source_version=None if source_variant is None else str(source_variant),
    )


def band_convolution_weights(
    band: SensorBand,
    spectrum_wavelength_nm: np.ndarray,
) -> np.ndarray:
    """Return rsrf-owned normalized band weights on a target spectral grid."""
    weights = _convolution_weights(
        np.asarray(spectrum_wavelength_nm, dtype=float),
        band_response_definition_input(band),
    )
    return np.asarray(weights, dtype=np.float32)


def load_band_rsrf(
    band: SensorBand,
    *,
    sensor_id: str,
    satellite_id: str,
    rsrf_root: Path | str | None = None,
) -> RelativeSpectralResponse:
    """Load a canonical RSRF for a SIAC band from the external RSRF catalog."""
    if band.rsrf_sensor_unit_id is None:
        raise ValueError(f"Band {band.name!r} does not define an RSRF sensor unit id")

    root = _normalize_rsrf_root(rsrf_root)
    response_definition = rsrf.load_response_definition(
        band.rsrf_sensor_unit_id,
        band.rsrf_band_id or band.name,
        band.rsrf_representation_variant or "band_average",
        root=root,
    )
    curve = _curve_from_response_definition(
        response_definition,
        band_id=str(band.rsrf_band_id or band.name),
        source_variant=band.rsrf_representation_variant or None,
    )
    source_variant = getattr(curve, "source_variant", None)

    return RelativeSpectralResponse.from_tabulated(
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        band_name=band.name,
        wavelengths_nm=np.asarray(curve.wavelength_nm, dtype=np.float32),
        response=np.asarray(curve.response, dtype=np.float32),
        source_id="RSRF",
        source_version=None if source_variant is None else str(source_variant),
    )


def load_sensor_config_with_rsrf(
    sensor_id: str,
    satellite_id: str,
    *,
    rsrf_root: Path | str | None = None,
    base_config: SensorConfig | None = None,
) -> SensorConfig:
    """Attach RSRF curves to a built-in sensor config."""
    resolved_base = base_config or get_sensor_config(sensor_id, satellite_id)
    bands: list[SensorBand] = []

    for band in resolved_base.bands:
        band_rsrf = load_band_rsrf(
            band,
            sensor_id=resolved_base.sensor_id,
            satellite_id=resolved_base.satellite_id,
            rsrf_root=rsrf_root,
        )
        bands.append(
            SensorBand(
                name=band.name,
                center_wavelength=float(
                    band_rsrf.effective_wavelength_nm or band.center_wavelength
                ),
                bandwidth=float(band_rsrf.fwhm_nm or band.bandwidth),
                resolution=band.resolution,
                band_index=band.band_index,
                rsrf_wavelengths_nm=band_rsrf.wavelengths_nm.copy(),
                rsrf_response=band_rsrf.response.copy(),
                rsrf_sensor_unit_id=band.rsrf_sensor_unit_id,
                rsrf_representation_variant=band.rsrf_representation_variant,
                rsrf_band_id=band.rsrf_band_id,
            )
        )

    return SensorConfig(
        sensor_id=resolved_base.sensor_id,
        satellite_id=resolved_base.satellite_id,
        bands=tuple(bands),
        default_ref_scale=resolved_base.default_ref_scale,
        default_ref_offset=resolved_base.default_ref_offset,
    )
