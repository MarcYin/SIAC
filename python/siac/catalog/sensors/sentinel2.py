"""Sentinel-2 MSI sensor catalog entries.

Sentinel-2A, -2B, and -2C share the same band layout (13 bands, identical
indices, resolutions, and bandwidths). Only a handful of band centers shift
between satellites, reflecting the per-satellite spectral characterisation
published by ESA. The ``_make_s2_config`` helper below centralises the shared
band specification and keeps the per-satellite differences in compact override
maps.

Reference: ESA "Sentinel-2 Spectral Response Functions (S2-SRF)" reports —
band-centre / FWHM tables in S2-PDGS-MPC-L2A-IODD differ slightly between S2A
and S2B/S2C. The numbers below preserve the values that were previously
hard-coded; this refactor only de-duplicates them.
"""

from __future__ import annotations

from typing import Any

from siac.catalog.sensors._common import rsrf_band
from siac.domain.sensors import SensorConfig

# Common (S2A) band layout. Each entry is the kwargs we pass to ``rsrf_band``
# minus ``sensor_unit_id`` (which is per-satellite). Band centers and FWHMs
# default to S2A values; per-satellite overrides only need to specify the
# fields that differ.
_S2_COMMON_BANDS: tuple[dict[str, Any], ...] = (
    {"name": "B01", "center_wavelength": 443.0, "bandwidth": 20.0, "resolution": 60.0, "band_index": 0},
    {"name": "B02", "center_wavelength": 490.0, "bandwidth": 65.0, "resolution": 10.0, "band_index": 1},
    {"name": "B03", "center_wavelength": 560.0, "bandwidth": 35.0, "resolution": 10.0, "band_index": 2},
    {"name": "B04", "center_wavelength": 665.0, "bandwidth": 30.0, "resolution": 10.0, "band_index": 3},
    {"name": "B05", "center_wavelength": 705.0, "bandwidth": 15.0, "resolution": 20.0, "band_index": 4},
    {"name": "B06", "center_wavelength": 740.0, "bandwidth": 15.0, "resolution": 20.0, "band_index": 5},
    {"name": "B07", "center_wavelength": 783.0, "bandwidth": 20.0, "resolution": 20.0, "band_index": 6},
    {"name": "B08", "center_wavelength": 842.0, "bandwidth": 115.0, "resolution": 10.0, "band_index": 7},
    {"name": "B8A", "center_wavelength": 865.0, "bandwidth": 20.0, "resolution": 20.0, "band_index": 8},
    {"name": "B09", "center_wavelength": 945.0, "bandwidth": 20.0, "resolution": 60.0, "band_index": 9},
    {"name": "B10", "center_wavelength": 1375.0, "bandwidth": 30.0, "resolution": 60.0, "band_index": 10},
    {"name": "B11", "center_wavelength": 1610.0, "bandwidth": 90.0, "resolution": 20.0, "band_index": 11},
    {"name": "B12", "center_wavelength": 2190.0, "bandwidth": 180.0, "resolution": 20.0, "band_index": 12},
)


def _make_s2_config(
    satellite_id: str,
    band_overrides: dict[str, dict[str, Any]] | None = None,
) -> SensorConfig:
    """Build a Sentinel-2 ``SensorConfig`` for a given satellite.

    ``band_overrides`` is a mapping of ``{band_name: {field: value}}`` patching
    the per-satellite differences (typically ``center_wavelength``). The base
    layout is taken from ``_S2_COMMON_BANDS``; the resulting RSRF
    ``sensor_unit_id`` is derived from ``satellite_id`` (e.g. ``"S2B"`` →
    ``"sentinel-2b_msi"``).

    Aerosol retrieval defaults to ``("B02", "B04")`` for every Sentinel-2
    satellite (visible blue + red anchor); see
    ``SensorConfig.default_aerosol_solver_bands``.
    """
    band_overrides = band_overrides or {}
    sensor_unit_id = f"sentinel-{satellite_id[1:].lower()}_msi"

    bands = []
    for spec in _S2_COMMON_BANDS:
        merged: dict[str, Any] = dict(spec)
        merged.update(band_overrides.get(merged["name"], {}))
        bands.append(rsrf_band(sensor_unit_id=sensor_unit_id, **merged))

    return SensorConfig(
        sensor_id="MSI",
        satellite_id=satellite_id,
        bands=tuple(bands),
        default_ref_scale=1.0 / 10000.0,
        default_ref_offset=0.0,
        # MSI-preferred aerosol bands: blue + red anchors. Moved into the
        # catalog (REVIEW.md §3.2) so the generic SensorConfig does not
        # hard-code "MSI" branching.
        aerosol_solver_band_names=("B02", "B04"),
    )


# S2A is the reference layout (zero overrides).
SENTINEL2A_CONFIG = _make_s2_config("S2A")

# S2B band centers as previously defined; per-satellite spectral
# characterisation (ESA S2-SRF / S2-PDGS-MPC-L2A-IODD).
SENTINEL2B_CONFIG = _make_s2_config(
    "S2B",
    band_overrides={
        "B01": {"center_wavelength": 442.0},
        "B02": {"center_wavelength": 492.0},
        "B03": {"center_wavelength": 559.0},
        "B05": {"center_wavelength": 704.0},
        "B06": {"center_wavelength": 739.0},
        "B07": {"center_wavelength": 780.0},
        "B08": {"center_wavelength": 833.0},
        "B8A": {"center_wavelength": 864.0},
        "B09": {"center_wavelength": 943.0},
        "B10": {"center_wavelength": 1377.0},
        "B12": {"center_wavelength": 2186.0},
    },
)

# S2C inherits the S2A band centers verbatim in the previously hard-coded
# config; no overrides needed.
SENTINEL2C_CONFIG = _make_s2_config("S2C")
