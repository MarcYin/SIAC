from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
import xarray as xr

import siac.algorithms.surface.brdf_whittaker as brdf_whittaker_mod
import siac.algorithms.surface.kernel_model as kernel_model_mod
import siac.algorithms.surface.spectral_mapping as spectral_mapping_mod
from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.spectral_mapping import (
    SpectralMapper,
    SpectralMappingConfig,
    convolve_hyperspectral_reflectance,
    map_multispectral_reflectance,
    needs_spectral_mapping,
)
from siac.domain import SensorBand
from siac.runtime import BRDFKernelWeights, GeometryAngles


@pytest.fixture(autouse=True)
def _stub_kernel_model_kernels(monkeypatch: pytest.MonkeyPatch) -> None:
    class _FakeBRDFKernels:
        def __init__(self, hb: float = 2.0, br: float = 1.0) -> None:
            self.hb = hb
            self.br = br

        def compute(self, vza, sza, raa):  # noqa: ANN001
            del sza, raa
            if isinstance(vza, xr.DataArray):
                zeros = xr.zeros_like(vza, dtype=np.float32)
                return zeros, zeros.copy(deep=False)
            shape = np.asarray(vza).shape
            zeros = np.zeros(shape, dtype=np.float32)
            return zeros, zeros.copy()

    monkeypatch.setattr(kernel_model_mod, "BRDFKernels", _FakeBRDFKernels)
    monkeypatch.setattr(brdf_whittaker_mod, "BRDFKernels", _FakeBRDFKernels)

    def _fake_whittaker_smooth_cube(cube, weights, temporal_lambda):  # noqa: ANN001
        del weights, temporal_lambda
        values = np.asarray(cube, dtype=np.float32).copy()
        smoothed = values.copy()
        time_axis = np.arange(values.shape[0], dtype=np.float32)
        for band_index in range(values.shape[1]):
            for y_index in range(values.shape[2]):
                for x_index in range(values.shape[3]):
                    series = values[:, band_index, y_index, x_index]
                    finite = np.isfinite(series)
                    if not np.any(finite):
                        continue
                    if finite.sum() == 1:
                        smoothed[:, band_index, y_index, x_index] = series[finite][0]
                        continue
                    smoothed[:, band_index, y_index, x_index] = np.interp(
                        time_axis,
                        time_axis[finite],
                        series[finite],
                    ).astype(np.float32)
        return smoothed

    monkeypatch.setattr(brdf_whittaker_mod, "whittaker_smooth_cube", _fake_whittaker_smooth_cube)


def _full_spectrum() -> tuple[np.ndarray, np.ndarray]:
    wavelengths = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)
    spectrum = np.clip(
        0.03
        + 0.04 * np.exp(-0.5 * ((wavelengths - 550.0) / 35.0) ** 2)
        - 0.03 * np.exp(-0.5 * ((wavelengths - 675.0) / 20.0) ** 2)
        + 0.45 / (1.0 + np.exp(-(wavelengths - 715.0) / 18.0)),
        0.0,
        0.9,
    )
    spectrum *= (
        1.0
        - 0.10 * np.exp(-0.5 * ((wavelengths - 970.0) / 25.0) ** 2)
        - 0.40 * np.exp(-0.5 * ((wavelengths - 1940.0) / 70.0) ** 2)
    )
    return wavelengths, np.asarray(spectrum, dtype=np.float32)


def _source_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand("Band3", 469.0, 20.0, 500.0, 0, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B3"),
        SensorBand("Band4", 555.0, 20.0, 500.0, 1, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B4"),
        SensorBand("Band2", 858.5, 35.0, 500.0, 2, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B2"),
        SensorBand("Band6", 1640.0, 24.0, 500.0, 3, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B6"),
        SensorBand("Band7", 2130.0, 50.0, 500.0, 4, rsrf_sensor_unit_id="terra_modis", rsrf_band_id="B7"),
    )


def _sampled_band(
    name: str,
    center_wavelength: float,
    bandwidth: float,
    resolution: float,
    band_index: int,
    *,
    sensor_unit_id: str,
    rsrf_band_id: str,
) -> SensorBand:
    wavelengths = np.arange(
        center_wavelength - 3.0 * bandwidth,
        center_wavelength + 3.0 * bandwidth + 1.0,
        1.0,
        dtype=np.float32,
    )
    sigma = bandwidth / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    response = np.exp(-0.5 * ((wavelengths - center_wavelength) / sigma) ** 2).astype(np.float32)
    return SensorBand(
        name,
        center_wavelength,
        bandwidth,
        resolution,
        band_index,
        rsrf_wavelengths_nm=wavelengths,
        rsrf_response=response,
        rsrf_sensor_unit_id=sensor_unit_id,
        rsrf_band_id=rsrf_band_id,
    )


def _target_bands() -> tuple[SensorBand, ...]:
    return (
        _sampled_band("B02", 490.0, 65.0, 10.0, 0, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B02"),
        _sampled_band("B03", 560.0, 35.0, 10.0, 1, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B03"),
        _sampled_band("B08", 842.0, 115.0, 10.0, 2, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B08"),
        _sampled_band("B11", 1610.0, 90.0, 20.0, 3, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B11"),
        _sampled_band("B12", 2190.0, 180.0, 20.0, 4, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B12"),
    )


def _unsupported_source_bands() -> tuple[SensorBand, ...]:
    return (
        SensorBand("M2", 445.0, 18.0, 500.0, 0, rsrf_sensor_unit_id="snpp_viirs", rsrf_band_id="M2"),
        SensorBand("M4", 555.0, 20.0, 500.0, 1, rsrf_sensor_unit_id="snpp_viirs", rsrf_band_id="M4"),
        SensorBand("M7", 865.0, 39.0, 500.0, 2, rsrf_sensor_unit_id="snpp_viirs", rsrf_band_id="M7"),
    )


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.4, dtype=np.float32), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.3, dtype=np.float32), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.7, dtype=np.float32), dims=["y", "x"]),
    )


def _with_geo(data: xr.DataArray) -> xr.DataArray:
    height = int(data.sizes["y"])
    width = int(data.sizes["x"])
    xmin = 399960.0
    ymax = 4700040.0
    resolution = 500.0
    xmax = xmin + width * resolution
    ymin = ymax - height * resolution
    transform = rasterio.transform.from_bounds(xmin, ymin, xmax, ymax, width, height)
    x = np.linspace(xmin + resolution / 2.0, xmax - resolution / 2.0, width, dtype=np.float64)
    y = np.linspace(ymax - resolution / 2.0, ymin + resolution / 2.0, height, dtype=np.float64)
    out = data.assign_coords({"y": y, "x": x}).rio.set_spatial_dims(x_dim="x", y_dim="y")
    return out.rio.write_crs("EPSG:32615").rio.write_transform(transform)


class _FakeSchema:
    def __init__(self, band_ids: tuple[str, ...]) -> None:
        self._band_ids = band_ids

    def band_ids(self) -> tuple[str, ...]:
        return self._band_ids


class _FakePublishedMapper:
    def __init__(
        self,
        *,
        prepared_root: Path,
        source_sensors: tuple[str, ...],
        target_rows_by_sensor: dict[str, np.ndarray] | None,
        calls: dict[str, object],
    ) -> None:
        self.prepared_root = prepared_root
        self.manifest = SimpleNamespace(source_sensors=source_sensors)
        self._calls = calls
        self._target_rows_by_sensor = {
            sensor_id: np.asarray(values, dtype=np.float64)
            for sensor_id, values in (target_rows_by_sensor or {}).items()
        }

    def get_sensor_schema(self, sensor_id: str) -> _FakeSchema:
        self._calls.setdefault("schemas", []).append(sensor_id)
        if sensor_id == "terra_modis":
            return _FakeSchema(("blue", "green", "red", "nir", "swir1", "swir2"))
        if sensor_id == "sentinel-2a_msi":
            return _FakeSchema(("ultra_blue", "blue", "green", "red", "nir", "swir1", "swir2"))
        if sensor_id == "snpp_viirs":
            return _FakeSchema(("blue", "green", "red", "nir", "swir1", "swir2"))
        raise KeyError(sensor_id)

    def map_reflectance_batch_arrays_ndarray(self, **kwargs):  # noqa: ANN003
        self._calls["batch"] = kwargs
        rows = int(np.asarray(kwargs["reflectance_rows"]).shape[0])
        target_sensor = str(kwargs.get("target_sensor"))
        target_rows = self._target_rows_by_sensor.get(target_sensor)
        if target_rows is None:
            raise KeyError(f"missing fake target rows for {target_sensor}")
        return SimpleNamespace(
            reflectance=np.repeat(target_rows[np.newaxis, :], rows, axis=0),
            source_fit_rmse=np.full(rows, 0.01, dtype=np.float32),
            output_columns=self.get_sensor_schema(target_sensor).band_ids(),
        )


def _install_fake_runtime(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    target_rows_by_sensor: dict[str, np.ndarray] | None = None,
    source_sensors: tuple[str, ...] = ("terra_modis",),
) -> dict[str, object]:
    calls: dict[str, object] = {}
    prepared_root = tmp_path / "published-runtime"
    prepared_root.mkdir()
    monkeypatch.setattr(spectral_mapping_mod, "resolve_prepared_library_root", lambda: prepared_root)

    def _factory(prepared_root_arg, verify_checksums=False):  # noqa: ANN001
        calls["init"] = {
            "prepared_root": Path(prepared_root_arg),
            "verify_checksums": verify_checksums,
        }
        return _FakePublishedMapper(
            prepared_root=Path(prepared_root_arg),
            source_sensors=source_sensors,
            target_rows_by_sensor=target_rows_by_sensor,
            calls=calls,
        )

    monkeypatch.setattr(spectral_mapping_mod, "PackageSpectralMapper", _factory)
    return calls


def _published_target_sensor_bands() -> tuple[SensorBand, ...]:
    return (
        _sampled_band("B01", 443.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B01"),
        _sampled_band("B02", 490.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B02"),
        _sampled_band("B03", 560.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B03"),
        _sampled_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B04"),
        _sampled_band("B08", 842.0, 115.0, 10.0, 4, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B08"),
        _sampled_band("B11", 1610.0, 90.0, 20.0, 5, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B11"),
        _sampled_band("B12", 2190.0, 180.0, 20.0, 6, sensor_unit_id="sentinel-2a_msi", rsrf_band_id="B12"),
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
        spectral_library=SpectralMappingConfig(),
    )

    np.testing.assert_allclose(mapped.values, source.values, atol=1e-6)
    np.testing.assert_allclose(mapped_unc.values, unc.values, atol=1e-6)
    assert needs_spectral_mapping(bands, bands) is False


def test_multispectral_mapping_projects_requested_target_bands_from_full_spectrum(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    wavelengths, spectrum = _full_spectrum()
    sample = xr.DataArray(
        spectrum.reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": wavelengths},
    )
    published_target = convolve_hyperspectral_reflectance(sample, wavelengths, _published_target_sensor_bands())
    calls = _install_fake_runtime(
        monkeypatch,
        tmp_path,
        target_rows_by_sensor={"sentinel-2a_msi": np.asarray(published_target.values[:, 0, 0], dtype=np.float64)},
    )
    source = convolve_hyperspectral_reflectance(sample, wavelengths, _source_bands())
    expected = convolve_hyperspectral_reflectance(sample, wavelengths, _target_bands())

    mapped, mapped_unc = map_multispectral_reflectance(
        source,
        source_bands=_source_bands(),
        target_bands=_target_bands(),
        spectral_library=SpectralMappingConfig(knn_backend="numpy"),
        k_neighbors=1,
    )

    np.testing.assert_allclose(mapped.values, expected.values, atol=1e-4)
    assert np.all(mapped_unc.values >= 0.0)
    assert calls["init"]["prepared_root"] == tmp_path / "published-runtime"
    assert calls["init"]["verify_checksums"] is False
    assert calls["batch"]["source_sensor"] == "terra_modis"
    assert calls["batch"]["output_mode"] == "target_sensor"
    assert calls["batch"]["target_sensor"] == "sentinel-2a_msi"
    np.testing.assert_allclose(
        np.asarray(calls["batch"]["reflectance_rows"][0], dtype=np.float64),
        np.array([source.values[0, 0, 0], source.values[1, 0, 0], np.nan, source.values[2, 0, 0], source.values[3, 0, 0], source.values[4, 0, 0]], dtype=np.float64),
        equal_nan=True,
    )
    assert needs_spectral_mapping(_source_bands(), _target_bands()) is True


def test_whittaker_route_uses_published_runtime_mapping(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    wavelengths, spectrum = _full_spectrum()
    sample = xr.DataArray(
        spectrum.reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": wavelengths},
    )
    published_target = convolve_hyperspectral_reflectance(sample, wavelengths, _published_target_sensor_bands())
    _install_fake_runtime(
        monkeypatch,
        tmp_path,
        target_rows_by_sensor={"sentinel-2a_msi": np.asarray(published_target.values[:, 0, 0], dtype=np.float64)},
    )
    expected = convolve_hyperspectral_reflectance(sample, wavelengths, _target_bands())
    source = convolve_hyperspectral_reflectance(sample, wavelengths, _source_bands())

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
        spectral_library=SpectralMappingConfig(),
        spectral_k_neighbors=1,
    )

    assert tuple(prior.boa.coords["band"].values.tolist()) == tuple(band.name for band in _target_bands())
    np.testing.assert_allclose(prior.boa.values, expected.values, atol=2.0e-2)
    assert np.all(prior.boa_unc.values > 0.0)


def test_spectral_mapper_handles_time_dimension(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    wavelengths, spectrum = _full_spectrum()
    sample = xr.DataArray(
        spectrum.reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": wavelengths},
    )
    published_target = convolve_hyperspectral_reflectance(sample, wavelengths, _published_target_sensor_bands())
    _install_fake_runtime(
        monkeypatch,
        tmp_path,
        target_rows_by_sensor={"sentinel-2a_msi": np.asarray(published_target.values[:, 0, 0], dtype=np.float64)},
    )
    source0 = convolve_hyperspectral_reflectance(sample, wavelengths, _source_bands())
    source1 = convolve_hyperspectral_reflectance(sample * 0.8, wavelengths, _source_bands())
    source = xr.concat(
        [source0, source1],
        dim=xr.IndexVariable("time", np.array(["2024-07-01", "2024-07-08"], dtype="datetime64[D]")),
    )
    source = _with_geo(source)
    unc = _with_geo(xr.full_like(source, 0.02))

    mapper = SpectralMapper(_source_bands(), _target_bands(), spectral_library=SpectralMappingConfig(), k_neighbors=1)
    mapped, mapped_unc, source_fit = mapper.map(source, source_uncertainty=unc)

    assert mapped.dims == ("time", "band", "y", "x")
    assert mapped_unc.shape == mapped.shape
    assert source_fit.dims == ("time", "y", "x")
    assert tuple(mapped.coords["band"].values.tolist()) == tuple(band.name for band in _target_bands())
    assert mapped.rio.crs is not None
    assert mapped.rio.crs.to_string() == "EPSG:32615"
    assert mapped_unc.rio.crs is not None
    assert mapped_unc.rio.crs.to_string() == "EPSG:32615"
    assert mapped.rio.transform(recalc=True) == source.rio.transform(recalc=True)
    assert mapped_unc.rio.transform(recalc=True) == source.rio.transform(recalc=True)


def test_mapping_rejects_non_config_runtime_input() -> None:
    with pytest.raises(TypeError, match="SpectralMappingConfig"):
        SpectralMapper(_source_bands(), _target_bands(), spectral_library=object())  # type: ignore[arg-type]


def test_mapping_requires_supported_published_source_sensor(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_fake_runtime(monkeypatch, tmp_path, source_sensors=("terra_modis",))

    with pytest.raises(ValueError, match="does not support the requested source sensor"):
        SpectralMapper(_unsupported_source_bands(), _target_bands(), spectral_library=SpectralMappingConfig())


def test_spectral_mapper_rejects_missing_band_dimension() -> None:
    mapper = SpectralMapper(_target_bands(), _target_bands(), spectral_library=SpectralMappingConfig(), k_neighbors=1)
    data = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])

    with pytest.raises(ValueError, match="'band' dimension"):
        mapper.map(data)


def test_spectral_mapper_rejects_missing_source_bands() -> None:
    mapper = SpectralMapper(_target_bands(), _target_bands(), spectral_library=SpectralMappingConfig(), k_neighbors=1)
    data = xr.DataArray(
        np.ones((1, 1, 1), dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": ["B02"], "y": [0], "x": [0]},
    )

    with pytest.raises(KeyError, match="missing source bands"):
        mapper.map(data)


def test_uncertainty_is_at_least_floor_for_mapped_bands(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    wavelengths, spectrum = _full_spectrum()
    sample = xr.DataArray(
        spectrum.reshape(1, 1, -1),
        dims=["y", "x", "wavelength"],
        coords={"y": [0], "x": [0], "wavelength": wavelengths},
    )
    published_target = convolve_hyperspectral_reflectance(sample, wavelengths, _published_target_sensor_bands())
    calls = _install_fake_runtime(
        monkeypatch,
        tmp_path,
        target_rows_by_sensor={"sentinel-2a_msi": np.asarray(published_target.values[:, 0, 0], dtype=np.float64)},
    )
    source = xr.DataArray(
        np.array([[[0.2]], [[0.3]], [[0.4]], [[0.25]], [[0.15]]], dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": [band.name for band in _source_bands()], "y": [0], "x": [0]},
    )

    def _zero_fit(**kwargs):  # noqa: ANN003
        calls["batch"] = kwargs
        return SimpleNamespace(
            reflectance=np.repeat(
                np.asarray(published_target.values[:, 0, 0], dtype=np.float64)[np.newaxis, :],
                1,
                axis=0,
            ),
            source_fit_rmse=np.array([0.0], dtype=np.float32),
            output_columns=("ultra_blue", "blue", "green", "red", "nir", "swir1", "swir2"),
        )

    mapper = SpectralMapper(_source_bands(), _target_bands(), spectral_library=SpectralMappingConfig(), k_neighbors=1)
    mapper._package_mapper = SimpleNamespace(map_reflectance_batch_arrays_ndarray=_zero_fit)
    mapped, mapped_unc, _source_fit = mapper.map(source)

    assert mapped.shape == (len(_target_bands()), 1, 1)
    assert np.all(mapped_unc.values >= np.float32(0.005))
