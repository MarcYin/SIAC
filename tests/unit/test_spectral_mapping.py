from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.spectral_mapping import (
    HyperspectralLibrary,
    SpectralMapper,
    convolve_hyperspectral_reflectance,
    map_multispectral_reflectance,
    needs_spectral_mapping,
)
from siac.domain import BRDFKernelWeights, GeometryAngles, SensorBand


def _library() -> HyperspectralLibrary:
    wavelengths = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)
    veg = np.clip(
        0.03
        + 0.04 * np.exp(-0.5 * ((wavelengths - 550.0) / 35.0) ** 2)
        - 0.03 * np.exp(-0.5 * ((wavelengths - 675.0) / 20.0) ** 2)
        + 0.45 / (1.0 + np.exp(-(wavelengths - 715.0) / 18.0)),
        0.0,
        0.9,
    )
    veg *= (
        1.0
        - 0.10 * np.exp(-0.5 * ((wavelengths - 970.0) / 25.0) ** 2)
        - 0.40 * np.exp(-0.5 * ((wavelengths - 1940.0) / 70.0) ** 2)
    )
    soil = np.clip(
        (0.09 + 1.2e-4 * (wavelengths - 400.0))
        * (1.0 - 0.03 * np.exp(-0.5 * ((wavelengths - 1900.0) / 80.0) ** 2)),
        0.0,
        0.7,
    )
    water = np.clip(0.02 * np.exp(-(wavelengths - 400.0) / 280.0), 0.0, 0.05)
    spectra = np.stack([veg.astype(np.float32), soil.astype(np.float32), water.astype(np.float32)], axis=0)
    return HyperspectralLibrary(
        wavelengths_nm=wavelengths,
        spectra=spectra,
        sample_ids=("veg", "soil", "water"),
        source_id="unit-test-library",
    )


def _source_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand("Band3", 469.0, 20.0, 500.0, 0),
        SensorBand("Band4", 555.0, 20.0, 500.0, 1),
        SensorBand("Band2", 858.5, 35.0, 500.0, 2),
        SensorBand("Band6", 1640.0, 24.0, 500.0, 3),
        SensorBand("Band7", 2130.0, 50.0, 500.0, 4),
    )


def _target_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B03", 560.0, 35.0, 10.0, 1),
        SensorBand("B08", 842.0, 115.0, 10.0, 2),
        SensorBand("B11", 1610.0, 90.0, 20.0, 3),
        SensorBand("B12", 2190.0, 180.0, 20.0, 4),
    )


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.4, dtype=np.float32), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.3, dtype=np.float32), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.7, dtype=np.float32), dims=["y", "x"]),
    )


def test_identity_mapping_keeps_reflectance_and_uncertainty() -> None:
    bands = _target_bands()
    coords = {"band": [band.name for band in bands], "y": [0], "x": [0]}
    source = xr.DataArray(
        np.array([[[0.12]], [[0.18]], [[0.40]], [[0.30]], [[0.22]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords=coords,
    )
    unc = xr.full_like(source, 0.03)

    mapped, mapped_unc = map_multispectral_reflectance(
        source,
        source_bands=bands,
        target_bands=bands,
        source_uncertainty=unc,
        spectral_library=_library(),
    )

    np.testing.assert_allclose(mapped.values, source.values, atol=1e-6)
    np.testing.assert_allclose(mapped_unc.values, unc.values, atol=1e-6)
    assert needs_spectral_mapping(bands, bands) is False


def test_multispectral_mapping_recovers_target_projection_for_library_sample() -> None:
    library = _library()
    spectrum = xr.DataArray(
        library.spectra[0].reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": library.wavelengths_nm},
    )
    source = convolve_hyperspectral_reflectance(spectrum, library.wavelengths_nm, _source_bands())
    target = convolve_hyperspectral_reflectance(spectrum, library.wavelengths_nm, _target_bands())

    mapped, mapped_unc = map_multispectral_reflectance(
        source,
        source_bands=_source_bands(),
        target_bands=_target_bands(),
        spectral_library=library,
        k_neighbors=1,
    )

    np.testing.assert_allclose(mapped.values, target.values, atol=1e-4)
    assert np.all(mapped_unc.values >= 0.0)
    assert needs_spectral_mapping(_source_bands(), _target_bands()) is True


def test_convolve_hyperspectral_reflectance_projects_to_target_bands() -> None:
    library = _library()
    cube = xr.DataArray(
        np.stack([library.spectra[1], library.spectra[2]], axis=0).reshape(2, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0, 1], "x": [0], "wavelength": library.wavelengths_nm},
    )
    projected = convolve_hyperspectral_reflectance(cube, library.wavelengths_nm, _target_bands())

    assert projected.dims == ("band", "y", "x")
    assert tuple(projected.coords["band"].values.tolist()) == tuple(band.name for band in _target_bands())
    assert projected.shape == (5, 2, 1)


def test_whittaker_route_maps_source_basis_to_target_basis() -> None:
    library = _library()
    sample = xr.DataArray(
        library.spectra[0].reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": library.wavelengths_nm},
    )
    expected = convolve_hyperspectral_reflectance(sample, library.wavelengths_nm, _target_bands())
    source = convolve_hyperspectral_reflectance(sample, library.wavelengths_nm, _source_bands())

    times = np.array(
        [(datetime(2024, 7, 15) + timedelta(days=offset)).date().isoformat() for offset in (-1, 0, 1)],
        dtype="datetime64[D]",
    )
    source_values = np.repeat(source.values[np.newaxis, ...], len(times), axis=0)
    source_unc = np.full_like(source_values, 0.02, dtype=np.float32)
    coords = {"time": times, "band": source.coords["band"], "y": [0], "x": [0]}
    weights = BRDFKernelWeights(
        f0=xr.DataArray(source_values, dims=["time", "band", "y", "x"], coords=coords),
        f1=xr.DataArray(np.zeros_like(source_values), dims=["time", "band", "y", "x"], coords=coords),
        f2=xr.DataArray(np.zeros_like(source_values), dims=["time", "band", "y", "x"], coords=coords),
        f0_unc=xr.DataArray(source_unc, dims=["time", "band", "y", "x"], coords=coords),
        f1_unc=xr.DataArray(source_unc, dims=["time", "band", "y", "x"], coords=coords),
        f2_unc=xr.DataArray(source_unc, dims=["time", "band", "y", "x"], coords=coords),
    )

    prior = BRDFWhittakerDeriver(temporal_lambda=10.0, apply_psf=False).compute_surface_prior(
        weights,
        _geometry((1, 1)),
        obs_time=datetime(2024, 7, 15),
        source_bands=_source_bands(),
        target_bands=_target_bands(),
        spectral_library=library,
        spectral_k_neighbors=1,
    )

    assert tuple(prior.boa.coords["band"].values.tolist()) == tuple(band.name for band in _target_bands())
    np.testing.assert_allclose(prior.boa.values, expected.values, atol=2.0e-2)
    assert np.all(prior.boa_unc.values > 0.0)


def test_spectral_mapper_handles_time_dimension() -> None:
    library = _library()
    sample0 = xr.DataArray(
        library.spectra[0].reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": library.wavelengths_nm},
    )
    sample1 = xr.DataArray(
        library.spectra[1].reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": library.wavelengths_nm},
    )
    source0 = convolve_hyperspectral_reflectance(sample0, library.wavelengths_nm, _source_bands())
    source1 = convolve_hyperspectral_reflectance(sample1, library.wavelengths_nm, _source_bands())
    source = xr.concat([source0, source1], dim=xr.IndexVariable("time", np.array(["2024-07-01", "2024-07-08"], dtype="datetime64[D]")))
    unc = xr.full_like(source, 0.02)

    mapper = SpectralMapper(_source_bands(), _target_bands(), spectral_library=library, k_neighbors=1)
    mapped, mapped_unc = mapper.map(source, source_uncertainty=unc)

    assert mapped.dims == ("time", "band", "y", "x")
    assert mapped_unc.shape == mapped.shape
    assert tuple(mapped.coords["band"].values.tolist()) == tuple(band.name for band in _target_bands())


def test_mapping_uses_rsrf_resolution_when_requested(monkeypatch) -> None:
    calls: list[tuple[str, str, str, str | None]] = []
    monkeypatch.setenv("RSRF_ROOT", "/tmp/fake-rsrf")

    class _FakeCurve:
        def __init__(self) -> None:
            self.wavelength_nm = np.array([480.0, 490.0, 500.0], dtype=np.float32)
            self.response = np.array([0.0, 1.0, 0.0], dtype=np.float32)

    def _fake_load_response_definition(
        sensor_unit_id: str,
        band_id: str,
        representation_variant: str,
        *,
        root=None,
    ) -> _FakeCurve:
        calls.append((sensor_unit_id, band_id, representation_variant, None if root is None else str(root)))
        return _FakeCurve()

    monkeypatch.setattr(
        "siac.algorithms.surface.spectral_mapping.rsrf.load_response_definition",
        _fake_load_response_definition,
    )

    source_band = SensorBand(
        "B02",
        490.0,
        65.0,
        10.0,
        0,
        rsrf_sensor_unit_id="sentinel-2a_msi",
        rsrf_representation_variant="band_average",
        rsrf_band_id="B02",
    )
    target_band = SensorBand(
        "B03",
        560.0,
        35.0,
        10.0,
        0,
        rsrf_sensor_unit_id="sentinel-2a_msi",
        rsrf_representation_variant="band_average",
        rsrf_band_id="B03",
    )
    coords = {"band": [source_band.name], "y": [0], "x": [0]}
    source = xr.DataArray(np.array([[[0.2]]], dtype=np.float32), dims=["band", "y", "x"], coords=coords)
    unc = xr.full_like(source, 0.01)

    mapped, mapped_unc = map_multispectral_reflectance(
        source,
        source_bands=(source_band,),
        target_bands=(target_band,),
        source_uncertainty=unc,
        spectral_library=_library(),
    )

    assert mapped.shape == (1, 1, 1)
    assert mapped_unc.shape == mapped.shape
    assert calls == [
        ("sentinel-2a_msi", "B02", "band_average", str(Path("/tmp/fake-rsrf").resolve())),
        ("sentinel-2a_msi", "B03", "band_average", str(Path("/tmp/fake-rsrf").resolve())),
    ]


def test_mapping_requires_explicit_spectral_library(monkeypatch) -> None:
    monkeypatch.delenv("SIAC_SPECTRAL_LIBRARY_ROOT", raising=False)

    coords = {"band": [band.name for band in _source_bands()], "y": [0], "x": [0]}
    source = xr.DataArray(
        np.array([[[0.12]], [[0.18]], [[0.40]], [[0.30]], [[0.22]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords=coords,
    )

    with np.testing.assert_raises_regex(ValueError, "explicit SIAC spectral library"):
        map_multispectral_reflectance(
            source,
            source_bands=_source_bands(),
            target_bands=_target_bands(),
            spectral_library=None,
        )
