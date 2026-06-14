"""Offline tests for the GEE-sourced MCD43 BRDF provider.

These build synthetic MCD43A1-style GeoTIFFs (the layout edown produces) and a
fake ``edown`` executable, so the provider's download orchestration, GeoTIFF
reading, QA decoding, and BRDFKernelWeights assembly are all exercised without
Earth Engine or network access.
"""

from __future__ import annotations

import os
import stat
import textwrap
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from siac.adapters.brdf.mcd43_gee import MCD43GEEProvider
from siac.adapters.earthdata_common import build_target_template
from siac.domain.sensors import SensorBand

# MODIS sinusoidal (sphere) — the CRS edown writes for MCD43A1.
_SINU = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs +type=crs"
_PARAM_FILL = 32767
_BANDS_TWO = ("Band1", "Band2")


def _band_descriptions(labels: tuple[str, ...]) -> list[str]:
    names: list[str] = []
    for label in labels:
        names += [f"BRDF_Albedo_Parameters_{label}_{k}" for k in ("iso", "vol", "geo")]
        names.append(f"BRDF_Albedo_Band_Mandatory_Quality_{label}")
    return names


# AOI shared by the read/assemble tests; the synthetic tile is placed to cover it.
_TEST_AOI = (2.0, 48.5, 2.05, 48.55)


def _write_synthetic_tiff(path: Path, labels: tuple[str, ...] = _BANDS_TWO) -> None:
    """Write a small MCD43A1-style multiband sinusoidal GeoTIFF.

    Parameter bands are int16 scaled by 1000 (iso≈0.20, vol≈0.05, geo≈0.02);
    QA bands are 0 (best quality). A nodata corner exercises masking. The tile is
    placed in MODIS sinusoidal space so it covers ``_TEST_AOI`` with margin.
    """
    from pyproj import Transformer

    descriptions = _band_descriptions(labels)
    res = 463.31
    to_sinu = Transformer.from_crs("EPSG:4326", _SINU, always_xy=True)
    xs, ys = to_sinu.transform(
        [_TEST_AOI[0], _TEST_AOI[2], _TEST_AOI[0], _TEST_AOI[2]],
        [_TEST_AOI[1], _TEST_AOI[1], _TEST_AOI[3], _TEST_AOI[3]],
    )
    margin = 10 * res
    origin_x = min(xs) - margin
    origin_y = max(ys) + margin
    width = int((max(xs) - min(xs)) / res) + 20
    height = int((max(ys) - min(ys)) / res) + 20
    transform = from_origin(origin_x, origin_y, res, res)
    data = np.zeros((len(descriptions), height, width), dtype=np.int16)
    for i, name in enumerate(descriptions):
        if name.endswith("_iso"):
            data[i] = 200
        elif name.endswith("_vol"):
            data[i] = 50
        elif name.endswith("_geo"):
            data[i] = 20
        else:  # Mandatory_Quality
            data[i] = 0
    # One nodata pixel in every parameter band.
    for i, name in enumerate(descriptions):
        if "Parameters" in name:
            data[i, 0, 0] = _PARAM_FILL

    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": len(descriptions),
        "dtype": "int16",
        "crs": _SINU,
        "transform": transform,
        "nodata": _PARAM_FILL,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data)
        for i, name in enumerate(descriptions, start=1):
            dst.set_band_description(i, name)


def _s2_bands() -> list[SensorBand]:
    return [
        SensorBand(
            name="B04", center_wavelength=665.0, bandwidth=30.0, resolution=10.0, band_index=0
        ),
        SensorBand(
            name="B8A", center_wavelength=865.0, bandwidth=20.0, resolution=20.0, band_index=1
        ),
    ]


def _target_template():
    # WGS84 AOI that the synthetic sinusoidal tile is built to cover.
    return _TEST_AOI, build_target_template(_TEST_AOI, "EPSG:4326", 0.01)


def test_resolve_requested_bands_matches_nearest_modis_band(tmp_path: Path) -> None:
    provider = MCD43GEEProvider(cache_dir=tmp_path)
    resolved = provider._resolve_requested_bands(_s2_bands())
    assert [(coord, definition.label) for coord, definition in resolved] == [
        ("B04", "Band1"),
        ("B8A", "Band2"),
    ]


def test_resolve_requested_bands_rejects_non_sensorband(tmp_path: Path) -> None:
    provider = MCD43GEEProvider(cache_dir=tmp_path)
    with pytest.raises(TypeError):
        provider._resolve_requested_bands(["B04"])  # type: ignore[list-item]


def test_source_bands_use_terra_modis_unit(tmp_path: Path) -> None:
    # The spectral-library RSRF lookup only accepts known source-sensor ids;
    # MCD43 must report "terra_modis" (matching the Earthaccess provider).
    provider = MCD43GEEProvider(cache_dir=tmp_path)
    units = {band.rsrf_sensor_unit_id for band in provider.source_bands}
    variants = {band.rsrf_representation_variant for band in provider.source_bands}
    assert units == {"terra_modis"}
    assert variants == {"band_average"}


def test_batch_validates_equal_lengths(tmp_path: Path) -> None:
    provider = MCD43GEEProvider(cache_dir=tmp_path)
    with pytest.raises(ValueError, match="same length"):
        provider.get_temporal_brdf_parameters_batch(
            bounds=_TEST_AOI,
            crs="EPSG:4326",
            obs_times=[datetime(2024, 2, 17)],
            target_resolution=0.01,
            bands=_s2_bands(),
            temporal_windows=[4, 4],
            sample_date_sets=[None],
        )


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("MODIS_061_MCD43A1_2024_02_16.tif", np.datetime64("2024-02-16", "D")),
        ("MODIS_061_MCD43A1_2024_12_31.tif", np.datetime64("2024-12-31", "D")),
        ("no_date_here.tif", None),
    ],
)
def test_date_from_tiff_name(name: str, expected) -> None:
    assert MCD43GEEProvider._date_from_tiff_name(name) == expected


def test_read_and_assemble_single_composite(tmp_path: Path) -> None:
    provider = MCD43GEEProvider(cache_dir=tmp_path)
    tiff = tmp_path / "MODIS_061_MCD43A1_2024_02_17.tif"
    _write_synthetic_tiff(tiff)
    bounds, template = _target_template()
    requested = provider._resolve_requested_bands(_s2_bands())

    layers = provider._read_tiff_to_band_layers(tiff, requested, template)
    assert set(layers) == {"B04", "B8A"}
    f0 = layers["B04"].sel(parameter="f0")
    # iso 200 * 0.001 == 0.20 where valid.
    assert np.isclose(np.nanmax(f0.values), 0.20, atol=1e-3)

    per_date = {np.datetime64("2024-02-17", "D"): layers}
    composite = provider._composite_nearest_valid(per_date, datetime(2024, 2, 17))
    weights = provider._weights_from_band_layers(composite, requested, template)
    assert weights.f0.dims == ("band", "y", "x")
    assert weights.f0.sizes["band"] == 2
    # Nadir reflectance == f0 ≈ 0.20.
    refl = weights.compute_reflectance(weights.f0 * 0, weights.f0 * 0)
    assert np.isclose(np.nanmax(refl.values), 0.20, atol=1e-3)
    # Uncertainty is finite where QA decoded.
    assert np.isfinite(weights.f0_unc.values).any()


def test_temporal_stack_has_time_dim(tmp_path: Path) -> None:
    provider = MCD43GEEProvider(cache_dir=tmp_path)
    bounds, template = _target_template()
    requested = provider._resolve_requested_bands(_s2_bands())
    per_date = {}
    for day in ("2024-02-16", "2024-02-17"):
        tiff = tmp_path / f"MODIS_061_MCD43A1_{day.replace('-', '_')}.tif"
        _write_synthetic_tiff(tiff)
        per_date[np.datetime64(day, "D")] = provider._read_tiff_to_band_layers(
            tiff, requested, template
        )
    time_axis = provider._time_axis(datetime(2024, 2, 17), 16)
    weights = provider._temporal_weights_from_layers(per_date, requested, template, time_axis)
    assert weights.f0.dims == ("time", "band", "y", "x")
    assert weights.f0.sizes["time"] == len(time_axis)
    populated = np.isfinite(weights.f0.values).any(axis=(1, 2, 3))
    assert populated.sum() == 2  # only the two supplied dates carry data


def _write_fake_edown(path: Path, payload_dir: Path) -> None:
    """A stand-in edown that writes a synthetic tiff into the output root."""
    script = textwrap.dedent(
        f"""
        #!/usr/bin/env python3
        import sys
        from pathlib import Path
        sys.path.insert(0, {str(Path(__file__).parent)!r})
        from test_mcd43_gee import _write_synthetic_tiff
        args = sys.argv
        out_root = Path(args[args.index("--output-root") + 1])
        images = out_root / "images" / "MODIS_061_MCD43A1"
        images.mkdir(parents=True, exist_ok=True)
        _write_synthetic_tiff(images / "MODIS_061_MCD43A1_2024_02_17.tif")
        print("fake edown wrote 1 image")
        """
    ).strip()
    path.write_text(script)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IRUSR)


def test_get_brdf_parameters_end_to_end_with_fake_edown(tmp_path: Path) -> None:
    fake_edown = tmp_path / "fake_edown"
    _write_fake_edown(fake_edown, tmp_path)
    provider = MCD43GEEProvider(
        cache_dir=tmp_path / "cache",
        edown_executable=f"{os.sys.executable} {fake_edown}".split()[0],
    )
    # The provider invokes [executable, "download", ...]; route through python.
    provider.edown_executable = str(fake_edown)

    bounds, _ = _target_template()
    weights = provider.get_brdf_parameters(
        bounds=bounds,
        crs="EPSG:4326",
        obs_time=datetime(2024, 2, 17),
        target_resolution=0.01,
        bands=_s2_bands(),
        temporal_window=4,
    )
    assert weights.f0.dims == ("band", "y", "x")
    assert np.isfinite(weights.f0.values).any()
    # Second call hits the populated cache (no re-download needed).
    again = provider.get_brdf_parameters(
        bounds=bounds,
        crs="EPSG:4326",
        obs_time=datetime(2024, 2, 17),
        target_resolution=0.01,
        bands=_s2_bands(),
        temporal_window=4,
    )
    assert np.isfinite(again.f0.values).any()


def test_qa_uncertainty_scales_with_reflectance() -> None:
    """QA uncertainty is reflectance-relative when a reflectance is supplied."""
    from siac.adapters.brdf._mcd43_qa import (
        _QA_REFLECTANCE_UNCERTAINTY_FLOOR as FLOOR,
    )
    from siac.adapters.brdf._mcd43_qa import (
        _QA_REFLECTANCE_UNCERTAINTY_FRACTION as FRAC,
    )
    from siac.adapters.brdf._mcd43_qa import (
        qa_values_to_uncertainty,
    )

    qa = np.array([0.0, 0.0, 1.0], dtype=np.float32)

    # No reflectance -> legacy absolute best-QA constant (earthaccess path).
    legacy = qa_values_to_uncertainty(qa)
    np.testing.assert_allclose(
        legacy, np.array([0.015, 0.015, 0.015 * 2.0**1.6], dtype=np.float32), rtol=1e-6
    )

    refl = np.array([0.03, 0.30, 0.03], dtype=np.float32)  # dark, bright, dark
    rel = qa_values_to_uncertainty(qa, reflectance=refl)
    dark0 = float(np.sqrt(FLOOR**2 + (FRAC * 0.03) ** 2))
    bright0 = float(np.sqrt(FLOOR**2 + (FRAC * 0.30) ** 2))
    # Dark target at best QA is ~5x tighter than the flat 0.015.
    assert rel[0] == pytest.approx(dark0, rel=1e-5)
    assert rel[0] < 0.004
    # Bright target stays near the legacy ~0.015.
    assert rel[1] == pytest.approx(bright0, rel=1e-5)
    assert 0.014 < rel[1] < 0.016
    # QA degradation still multiplies the (now reflectance-scaled) best-QA value.
    assert rel[2] == pytest.approx(dark0 * 2.0**1.6, rel=1e-5)

    # Missing reflectance falls back to the absolute best-QA value.
    fallback = qa_values_to_uncertainty(
        np.array([0.0], dtype=np.float32), reflectance=np.array([np.nan], dtype=np.float32)
    )
    assert fallback[0] == pytest.approx(0.015, rel=1e-6)
