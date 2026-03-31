from __future__ import annotations

import json
from datetime import datetime
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
import xarray as xr

import siac.algorithms.surface.swir_refine as swir_refine_mod
from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    MonthlyKernelWeightComposite,
)
from siac.algorithms.surface.brdf_monthly_database import build_monthly_composite_database
from siac.algorithms.surface.spectral_mapping import HyperspectralLibrary
from siac.algorithms.surface.swir_refine import (
    _forward_model_monthly_reflectance,
    _weekly_sample_dates,
    build_monthly_composites_from_brdf,
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
)
from siac.domain import SensorBand, SensorConfig
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    RTCoefficients,
)


class _IdentityRTModel:
    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
        shape = geometry.sza.shape
        dims = geometry.sza.dims
        coords = geometry.sza.coords
        one = xr.DataArray(np.ones(shape, dtype=np.float32), dims=dims, coords=coords)
        zero = xr.DataArray(np.zeros(shape, dtype=np.float32), dims=dims, coords=coords)
        return RTCoefficients(xap=one, xbp=zero, xcp=zero)

    def supports_jacobian(self) -> bool:
        return False

    @property
    def backend_name(self) -> str:
        return "identity"

    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
        return True


def _sensor_config() -> SensorConfig:
    return SensorConfig(
        sensor_id="MSI",
        satellite_id="S2A",
        bands=(
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
            SensorBand("B08", 842.0, 115.0, 10.0, 2),
            SensorBand("B11", 1610.0, 90.0, 20.0, 3),
            SensorBand("B12", 2190.0, 180.0, 20.0, 4),
        ),
    )


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


def _with_geo(data: xr.DataArray) -> xr.DataArray:
    height = int(data.sizes["y"])
    width = int(data.sizes["x"])
    resolution = 20.0
    xmin = 600000.0
    ymax = 4200000.0
    x = np.linspace(xmin + resolution / 2.0, xmin + width * resolution - resolution / 2.0, width, dtype=np.float64)
    y = np.linspace(ymax - resolution / 2.0, ymax - height * resolution + resolution / 2.0, height, dtype=np.float64)
    transform = rasterio.transform.from_origin(xmin, ymax, resolution, resolution)
    out = data.assign_coords({"x": x, "y": y}).rio.set_spatial_dims(x_dim="x", y_dim="y")
    return out.rio.write_crs("EPSG:32632").rio.write_transform(transform)


def _atmo(shape: tuple[int, int]) -> AtmosphericState:
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.2), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
    )


def _database():
    bands = ["B02", "B03", "B08", "B11", "B12"]
    composites: list[MonthlyBestPixelComposite] = []
    for month_idx in range(15):
        cube = np.zeros((len(bands), 2, 1), dtype=np.float32)
        cube[bands.index("B08"), 0, 0] = 0.40 + 0.01 * month_idx
        cube[bands.index("B11"), 0, 0] = 0.30 + 0.01 * month_idx
        cube[bands.index("B12"), 0, 0] = 0.20 + 0.01 * month_idx
        cube[bands.index("B02"), 0, 0] = 0.10 + 0.001 * month_idx
        cube[bands.index("B03"), 0, 0] = 0.15 + 0.001 * month_idx

        cube[bands.index("B08"), 1, 0] = 0.50 + 0.01 * month_idx
        cube[bands.index("B11"), 1, 0] = 0.40 + 0.01 * month_idx
        cube[bands.index("B12"), 1, 0] = 0.30 + 0.01 * month_idx
        cube[bands.index("B02"), 1, 0] = 0.40 + 0.001 * month_idx
        cube[bands.index("B03"), 1, 0] = 0.45 + 0.001 * month_idx

        composites.append(
            MonthlyBestPixelComposite(
                reflectance=xr.DataArray(
                    cube,
                    dims=["band", "y", "x"],
                    coords={"band": bands, "y": [0, 1], "x": [0]},
                ),
                quality=xr.DataArray(
                    np.zeros((2, 1), dtype=np.int16),
                    dims=["y", "x"],
                    coords={"y": [0, 1], "x": [0]},
                ),
                sample_index=xr.DataArray(
                    np.full((2, 1), month_idx, dtype=np.int16),
                    dims=["y", "x"],
                    coords={"y": [0, 1], "x": [0]},
                ),
                year=2020 + month_idx // 3,
                month=(month_idx % 3) + 6,
            )
        )
    return build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
    )


def _spectral_library() -> HyperspectralLibrary:
    wavelengths = np.arange(400.0, 2501.0, 1.0, dtype=np.float32)
    veg = np.clip(
        0.03
        + 0.04 * np.exp(-0.5 * ((wavelengths - 550.0) / 35.0) ** 2)
        - 0.03 * np.exp(-0.5 * ((wavelengths - 675.0) / 20.0) ** 2)
        + 0.45 / (1.0 + np.exp(-(wavelengths - 715.0) / 18.0)),
        0.0,
        0.9,
    )
    soil = np.clip(
        (0.09 + 1.2e-4 * (wavelengths - 400.0))
        * (1.0 - 0.03 * np.exp(-0.5 * ((wavelengths - 1900.0) / 80.0) ** 2)),
        0.0,
        0.7,
    )
    water = np.clip(
        0.02 * np.exp(-(wavelengths - 400.0) / 280.0)
        * (
            1.0
            - 0.65 * np.exp(-0.5 * ((wavelengths - 740.0) / 45.0) ** 2)
            - 0.95 * np.exp(-0.5 * ((wavelengths - 1200.0) / 70.0) ** 2)
        ),
        0.0,
        0.08,
    )
    return HyperspectralLibrary(
        wavelengths_nm=wavelengths,
        spectra=np.stack([veg, soil, water]).astype(np.float32),
        sample_ids=("veg", "soil", "water"),
    )


def _patch_kernel_to_reflectance(monkeypatch: pytest.MonkeyPatch) -> None:
    def _fake_reflectance_from_kernel_weights(weights, geometry):  # noqa: ANN001
        del geometry
        unc = weights.reflectance_unc if weights.reflectance_unc is not None else weights.f0_unc
        return weights.f0.astype(np.float32), unc.astype(np.float32)

    monkeypatch.setattr(
        "siac.algorithms.surface.swir_refine._reflectance_from_kernel_weights",
        _fake_reflectance_from_kernel_weights,
    )


def _patch_zero_kernel_geometry(monkeypatch: pytest.MonkeyPatch) -> None:
    class _ZeroKernels:
        def __init__(self, hb=2.0, br=1.0):  # noqa: ANN001
            self.hb = hb
            self.br = br

        def compute(self, vza, sza, raa):  # noqa: ANN001
            return xr.zeros_like(vza), xr.zeros_like(sza)

    monkeypatch.setattr("siac.algorithms.surface.swir_refine.BRDFKernels", _ZeroKernels)


def test_query_surface_prior_from_monthly_database_returns_visible_surface_prior() -> None:
    sensor_config = _sensor_config()
    toa = xr.Dataset(
        {
            "B02": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
            "B03": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
            "B08": xr.DataArray(np.array([[0.47], [0.57]], dtype=np.float32), dims=["y", "x"]),
            "B11": xr.DataArray(np.array([[0.37], [0.47]], dtype=np.float32), dims=["y", "x"]),
            "B12": xr.DataArray(np.array([[0.27], [0.37]], dtype=np.float32), dims=["y", "x"]),
        }
    )
    obs = ObservationBundle(
        toa=toa,
        geometry=_geometry((2, 1)),
        cloud_mask=xr.DataArray(np.zeros((2, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 2.0),
    )

    prior = query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_IdentityRTModel(),
        database=_database(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        k_neighbors=1,
    )

    assert prior.kernels is None
    assert prior.boa.dims == ("band", "y", "x")
    assert tuple(prior.boa.coords["band"].values.tolist()) == ("B02", "B03")
    np.testing.assert_allclose(
        prior.boa.sel(band="B02").values[:, 0],
        np.array([0.107, 0.407], dtype=np.float32),
        atol=1e-6,
    )
    np.testing.assert_allclose(
        prior.boa.sel(band="B03").values[:, 0],
        np.array([0.157, 0.457], dtype=np.float32),
        atol=1e-6,
    )
    assert prior.mask.shape == (2, 1)
    assert bool(prior.mask.values.all())


def test_query_surface_prior_from_monthly_database_resamples_coarse_atmo_prior() -> None:
    sensor_config = _sensor_config()
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((4, 2), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((4, 2), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((4, 2), 0.47, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((4, 2), 0.37, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((4, 2), 0.27, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((4, 2)),
        cloud_mask=xr.DataArray(np.zeros((4, 2), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 2.0, 4.0),
    )

    class _ShapeCheckingRTModel(_IdentityRTModel):
        def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
            assert atmo_state.aot.shape == geometry.sza.shape
            return super().compute_coefficients(geometry, atmo_state, band, compute_jacobian)

    prior = query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_ShapeCheckingRTModel(),
        database=_database(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        k_neighbors=1,
    )

    assert prior.boa.shape == (2, 2, 1)


def test_query_surface_prior_from_monthly_database_uses_template_backed_area_resampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sensor_config = _sensor_config()
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((4, 2), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((4, 2), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((4, 2), 0.47, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((4, 2), 0.37, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((4, 2), 0.27, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((4, 2)),
        cloud_mask=xr.DataArray(np.zeros((4, 2), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 2.0, 4.0),
    )
    calls: list[tuple[str, bool]] = []
    original = swir_refine_mod._resample_da

    def _capture(da, target_shape, method="bilinear", *, template=None):  # noqa: ANN001
        calls.append((str(method), template is not None))
        return original(da, target_shape, method, template=template)

    monkeypatch.setattr(swir_refine_mod, "_resample_da", _capture)

    prior = query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_IdentityRTModel(),
        database=_database(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        k_neighbors=1,
    )

    area_calls = [template_present for method, template_present in calls if method == "area"]
    assert area_calls
    assert all(area_calls)
    assert prior.boa.shape == (2, 2, 1)


def test_query_surface_prior_from_monthly_database_filters_high_composite_quality() -> None:
    sensor_config = _sensor_config()
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((2, 1), 0.47, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((2, 1), 0.37, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((2, 1), 0.27, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((2, 1)),
        cloud_mask=xr.DataArray(np.zeros((2, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 2.0),
    )
    database = _database()

    class _QualityFilteringDatabase:
        query_band_names = database.query_band_names
        visible_band_names = database.visible_band_names
        median_summary = database.median_summary

        def predict_visible(self, corrected_reflectance, *, k_neighbors=3):  # noqa: ANN001
            del corrected_reflectance, k_neighbors
            coords = {"band": ["B02", "B03"], "y": [0, 1], "x": [0]}
            predicted = xr.DataArray(
                np.array([[[0.1], [0.2]], [[0.15], [0.25]]], dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            )
            predicted_unc = xr.DataArray(
                np.full((2, 2, 1), 0.01, dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            )
            predicted_quality = xr.DataArray(
                np.array([[0.02], [0.08]], dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            )
            return predicted, predicted_unc, predicted_quality

    prior = query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_IdentityRTModel(),
        database=_QualityFilteringDatabase(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        k_neighbors=1,
        max_composite_quality=0.05,
    )

    np.testing.assert_array_equal(prior.mask.values[:, 0], np.array([True, False]))


def test_query_surface_prior_from_monthly_database_filters_high_prediction_uncertainty() -> None:
    sensor_config = _sensor_config()
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((2, 1), 0.47, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((2, 1), 0.37, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((2, 1), 0.27, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((2, 1)),
        cloud_mask=xr.DataArray(np.zeros((2, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 2.0),
    )
    database = _database()

    class _UncertaintyFilteringDatabase:
        query_band_names = database.query_band_names
        visible_band_names = database.visible_band_names
        median_summary = database.median_summary

        def predict_visible(self, corrected_reflectance, *, k_neighbors=3):  # noqa: ANN001
            del corrected_reflectance, k_neighbors
            coords = {"band": ["B02", "B03"], "y": [0, 1], "x": [0]}
            predicted = xr.DataArray(
                np.array([[[0.1], [0.2]], [[0.15], [0.25]]], dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            )
            predicted_unc = xr.DataArray(
                np.array([[[0.01], [0.06]], [[0.01], [0.06]]], dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            )
            predicted_quality = xr.DataArray(
                np.full((2, 1), 0.02, dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            )
            return predicted, predicted_unc, predicted_quality

    prior = query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_IdentityRTModel(),
        database=_UncertaintyFilteringDatabase(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        k_neighbors=1,
        max_prediction_uncertainty=0.05,
    )

    np.testing.assert_array_equal(prior.mask.values[:, 0], np.array([True, False]))


def test_query_surface_prior_from_monthly_database_filters_high_knn_feature_distance() -> None:
    sensor_config = _sensor_config()
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((2, 1), 0.47, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((2, 1), 0.37, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((2, 1), 0.27, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((2, 1)),
        cloud_mask=xr.DataArray(np.zeros((2, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 2.0),
    )
    database = _database()

    class _DistanceFilteringDatabase:
        query_band_names = database.query_band_names
        visible_band_names = database.visible_band_names
        median_summary = database.median_summary

        def predict_visible_with_diagnostics(self, corrected_reflectance, *, k_neighbors=3):  # noqa: ANN001
            del corrected_reflectance, k_neighbors
            coords = {"band": ["B02", "B03"], "y": [0, 1], "x": [0]}
            predicted = xr.DataArray(
                np.array([[[0.1], [0.2]], [[0.15], [0.25]]], dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            )
            predicted_unc = xr.DataArray(
                np.full((2, 2, 1), 0.01, dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            )
            predicted_quality = xr.DataArray(
                np.full((2, 1), 0.02, dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            )
            predicted_source_fit = xr.DataArray(
                np.zeros((2, 1), dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            )
            predicted_distance = xr.DataArray(
                np.array([[0.01], [0.09]], dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            )
            return SimpleNamespace(
                predicted=predicted,
                uncertainty=predicted_unc,
                quality=predicted_quality,
                source_fit_rmse=predicted_source_fit,
                knn_feature_distance=predicted_distance,
            )

    prior = query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_IdentityRTModel(),
        database=_DistanceFilteringDatabase(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        k_neighbors=1,
        max_knn_feature_distance=0.05,
    )

    np.testing.assert_array_equal(prior.mask.values[:, 0], np.array([True, False]))


def test_query_surface_prior_from_monthly_database_caches_distance_metrics(tmp_path) -> None:  # noqa: ANN001
    sensor_config = _sensor_config()
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((2, 1), 0.0, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((2, 1), 0.47, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((2, 1), 0.37, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((2, 1), 0.27, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((2, 1)),
        cloud_mask=xr.DataArray(np.zeros((2, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 2.0),
    )
    database = _database()

    class _DiagnosticDatabase:
        query_band_names = database.query_band_names
        visible_band_names = database.visible_band_names
        median_summary = database.median_summary

        def predict_visible_with_diagnostics(self, corrected_reflectance, *, k_neighbors=3):  # noqa: ANN001
            del corrected_reflectance, k_neighbors
            coords = {"band": ["B02", "B03"], "y": [0, 1], "x": [0]}
            predicted = _with_geo(xr.DataArray(
                np.array([[[0.1], [0.2]], [[0.15], [0.25]]], dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            ))
            predicted_unc = _with_geo(xr.DataArray(
                np.full((2, 2, 1), 0.01, dtype=np.float32),
                dims=["band", "y", "x"],
                coords=coords,
            ))
            predicted_quality = _with_geo(xr.DataArray(
                np.full((2, 1), 0.02, dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            ))
            predicted_source_fit = _with_geo(xr.DataArray(
                np.array([[0.01], [0.03]], dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            ))
            predicted_distance = _with_geo(xr.DataArray(
                np.array([[0.01], [0.09]], dtype=np.float32),
                dims=["y", "x"],
                coords={"y": [0, 1], "x": [0]},
            ))
            return SimpleNamespace(
                predicted=predicted,
                uncertainty=predicted_unc,
                quality=predicted_quality,
                source_fit_rmse=predicted_source_fit,
                knn_feature_distance=predicted_distance,
            )

    query_surface_prior_from_monthly_database(
        observation=obs,
        atmo_prior=_atmo((2, 1)),
        rt_model=_IdentityRTModel(),
        database=_DiagnosticDatabase(),
        query_band_names=("B08", "B11", "B12"),
        visible_band_names=("B02", "B03"),
        diagnostic_cache_dir=tmp_path,
    )

    diagnostics_dir = tmp_path / "diagnostics"
    data_paths = sorted(diagnostics_dir.glob("swir_refine_distances_*.tif"))
    metadata_paths = sorted(diagnostics_dir.glob("swir_refine_distances_*.json"))
    assert len(data_paths) == 2
    assert len(metadata_paths) == 1

    expected = {
        "source_fit_rmse": np.array([0.01, 0.03], dtype=np.float32),
        "knn_feature_distance": np.array([0.01, 0.09], dtype=np.float32),
    }
    expected_transform = rasterio.transform.from_origin(600000.0, 4200000.0, 20.0, 20.0)
    for path in data_paths:
        metric_name = path.stem.split("_", 4)[-1]
        with rasterio.open(path) as src:
            assert src.crs.to_string() == "EPSG:32632"
            assert src.transform == expected_transform
            np.testing.assert_allclose(src.read(1)[:, 0], expected[metric_name])

    metadata = json.loads(metadata_paths[0].read_text(encoding="utf-8"))
    assert metadata["query_band_names"] == ["B08", "B11", "B12"]
    assert metadata["visible_band_names"] == ["B02", "B03"]
    assert sorted(metadata["metrics"]) == sorted(expected)


def test_weekly_sample_dates_use_weekly_spacing() -> None:
    sample_dates = _weekly_sample_dates(2024, 1)
    assert [dt.day for dt in sample_dates] == [1, 8, 15, 22, 29]


def test_build_monthly_composites_from_brdf_honors_explicit_period_selection() -> None:
    class _FakeBRDFProvider:
        def __init__(self) -> None:
            self.batch_calls: list[dict[str, object]] = []
            self.source_bands = tuple(_sensor_config().bands)

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            self.batch_calls.append(kwargs)
            outputs = []
            bands = [band.name for band in kwargs["bands"]]
            for sample_dates in kwargs["sample_date_sets"]:
                sample_dates = tuple(sample_dates)
                n_time = len(sample_dates)
                data = np.full((n_time, len(bands), 2, 1), 0.2, dtype=np.float32)
                unc = np.full_like(data, 0.03, dtype=np.float32)
                coords = {
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                    "band": bands,
                    "y": [0, 1],
                    "x": [0],
                }
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    provider = _FakeBRDFProvider()

    collection = build_monthly_composites_from_brdf(
        brdf_provider=provider,
        bounds=(0.0, 0.0, 1.0, 2.0),
        crs="EPSG:32632",
        resolution=500.0,
        year_months=((2023, 8), (2022, 7)),
    )

    assert [(composite.year, composite.month) for composite in collection.composites] == [(2022, 7), (2023, 8)]
    assert len(provider.batch_calls) == 1
    requested_sample_sets = [tuple(sample_dates) for sample_dates in provider.batch_calls[0]["sample_date_sets"]]
    assert [(sample_dates[0].year, sample_dates[0].month) for sample_dates in requested_sample_sets] == [(2022, 7), (2023, 8)]


def test_build_monthly_composites_from_brdf_rejects_duplicate_periods() -> None:
    class _FakeBRDFProvider:
        source_bands = tuple(_sensor_config().bands)

    with pytest.raises(ValueError, match="Duplicate monthly composite selection"):
        build_monthly_composites_from_brdf(
            brdf_provider=_FakeBRDFProvider(),
            bounds=(0.0, 0.0, 1.0, 2.0),
            crs="EPSG:32632",
            resolution=500.0,
            year_months=((2023, 7), (2023, 7)),
        )


def test_build_monthly_surface_prior_database_requests_weekly_brdf_samples(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)

    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((2, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((2, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((2, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((2, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((2, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((2, 1)),
        cloud_mask=xr.DataArray(np.zeros((2, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 2.0),
    )

    class _FakeBRDFProvider:
        def __init__(self) -> None:
            self.batch_calls: list[dict[str, object]] = []
            self.source_bands = tuple(_sensor_config().bands)

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            self.batch_calls.append(kwargs)
            outputs = []
            bands = [band.name for band in kwargs["bands"]]
            for sample_dates in kwargs["sample_date_sets"]:
                sample_dates = tuple(sample_dates)
                n_time = len(sample_dates)
                data = np.full((n_time, len(bands), 2, 1), 0.2, dtype=np.float32)
                unc = np.full_like(data, 0.03, dtype=np.float32)
                coords = {
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                    "band": bands,
                    "y": [0, 1],
                    "x": [0],
                }
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    provider = _FakeBRDFProvider()
    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]

    database = build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=provider,
        resolution=500.0,
        geometry=_geometry((2, 1)),
        visible_bands=visible_bands,
        query_bands=query_bands,
    )

    assert database.entries_features.shape[0] == 15 * 2 * 1
    assert len(provider.batch_calls) == 5
    assert all(len(call["sample_date_sets"]) <= 3 for call in provider.batch_calls)
    all_sample_sets = [
        tuple(sample_dates)
        for call in provider.batch_calls
        for sample_dates in call["sample_date_sets"]
    ]
    assert len(all_sample_sets) == 15
    july_call = next(sample_dates for sample_dates in all_sample_sets if sample_dates[0].year == 2023 and sample_dates[0].month == 7)
    assert [dt.day for dt in july_call] == [1, 8, 15, 22, 29]


def test_forward_model_monthly_reflectance_uses_qa_based_reflectance_uncertainty(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_zero_kernel_geometry(monkeypatch)

    coords = {
        "time": np.array(["2024-07-01", "2024-07-08"], dtype="datetime64[D]"),
        "band": ["B02", "B08"],
        "y": [0],
        "x": [0],
    }
    data = np.full((2, 2, 1, 1), 0.2, dtype=np.float32)
    huge_unc = np.full((2, 2, 1, 1), 9.0, dtype=np.float32)
    reflectance_unc = np.array(
        [
            [[[0.015]], [[0.0454734]]],
            [[[0.08]], [[0.10]]],
        ],
        dtype=np.float32,
    )

    weights = BRDFKernelWeights(
        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
        f0_unc=xr.DataArray(huge_unc, dims=["time", "band", "y", "x"], coords=coords),
        f1_unc=xr.DataArray(huge_unc, dims=["time", "band", "y", "x"], coords=coords),
        f2_unc=xr.DataArray(huge_unc, dims=["time", "band", "y", "x"], coords=coords),
        reflectance_unc=xr.DataArray(reflectance_unc, dims=["time", "band", "y", "x"], coords=coords),
    )

    _reflectance, quality, _reflectance_unc = _forward_model_monthly_reflectance(
        weights,
        geometry=_geometry((1, 1)),
        year=2024,
        month=7,
    )

    np.testing.assert_allclose(
        quality.values[:, 0, 0],
        np.array([0.0302367, 0.09], dtype=np.float32),
        rtol=1e-5,
    )


def test_build_monthly_surface_prior_database_preserves_target_grid_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import rioxarray  # noqa: F401

    _patch_kernel_to_reflectance(monkeypatch)

    shape = (4, 4)
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full(shape, 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry(shape),
        cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(300000.0, 5500000.0, 302000.0, 5502000.0),
    )
    geometry = swir_refine_mod.resample_geometry_for_surface_prior(obs, resolution=500.0)

    class _CoarseBRDFProvider:
        def __init__(self) -> None:
            self.source_bands = tuple(_sensor_config().bands)

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            outputs = []
            for sample_dates in kwargs["sample_date_sets"]:
                coords = {
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                    "band": [band.name for band in self.source_bands],
                    "y": [0.0],
                    "x": [0.0],
                }
                base = np.array([0.1, 0.12, 0.42, 0.31, 0.21], dtype=np.float32).reshape(1, 5, 1, 1)
                data = np.repeat(base, len(sample_dates), axis=0)
                unc = np.full_like(data, 0.02)
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]
    database = build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=_CoarseBRDFProvider(),
        resolution=500.0,
        geometry=geometry,
        visible_bands=visible_bands,
        query_bands=query_bands,
    )

    assert database.median_summary.rio.crs is not None
    assert database.median_summary.rio.crs.to_string() == "EPSG:32632"
    assert tuple(database.median_summary.coords["x"].values.tolist()) == pytest.approx(
        tuple(geometry.sza.coords["x"].values.tolist())
    )
    assert tuple(database.median_summary.coords["y"].values.tolist()) == pytest.approx(
        tuple(geometry.sza.coords["y"].values.tolist())
    )
    assert database.median_summary.rio.transform(recalc=True) == geometry.sza.rio.transform(recalc=True)


def test_build_monthly_surface_prior_database_maps_source_basis_to_target_basis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)

    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((1, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((1, 1)),
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 1.0),
    )

    source_bands = (
        SensorBand("Band3", 469.0, 20.0, 500.0, 0),
        SensorBand("Band4", 555.0, 20.0, 500.0, 1),
        SensorBand("Band2", 858.5, 35.0, 500.0, 2),
        SensorBand("Band6", 1640.0, 24.0, 500.0, 3),
        SensorBand("Band7", 2130.0, 50.0, 500.0, 4),
    )

    class _MappedSourceBRDFProvider:
        def __init__(self) -> None:
            self.source_bands = source_bands

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            outputs = []
            assert tuple(band.name for band in kwargs["bands"]) == tuple(band.name for band in source_bands)
            for sample_dates in kwargs["sample_date_sets"]:
                sample_dates = tuple(sample_dates)
                coords = {
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                    "band": [band.name for band in source_bands],
                    "y": [0],
                    "x": [0],
                }
                base = np.array([0.08, 0.12, 0.42, 0.30, 0.22], dtype=np.float32).reshape(1, 5, 1, 1)
                data = np.repeat(base, len(sample_dates), axis=0)
                unc = np.full_like(data, 0.02)
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]
    database = build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=_MappedSourceBRDFProvider(),
        resolution=500.0,
        geometry=_geometry((1, 1)),
        visible_bands=visible_bands,
        query_bands=query_bands,
        spectral_library=_spectral_library(),
    )

    assert database.query_band_names == ("B08", "B11", "B12")
    assert database.visible_band_names == ("B02", "B03")
    assert np.isfinite(database.entries_features).all()
    assert np.isfinite(database.entries_visible).all()


def test_build_monthly_surface_prior_database_maps_kernel_composites_to_target_basis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)

    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((1, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((1, 1)),
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 1.0),
    )

    source_bands = (
        SensorBand("Band3", 469.0, 20.0, 500.0, 0),
        SensorBand("Band4", 555.0, 20.0, 500.0, 1),
        SensorBand("Band2", 858.5, 35.0, 500.0, 2),
        SensorBand("Band6", 1640.0, 24.0, 500.0, 3),
        SensorBand("Band7", 2130.0, 50.0, 500.0, 4),
    )

    class _TwoSampleProvider:
        def __init__(self) -> None:
            self.source_bands = source_bands

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            outputs = []
            for sample_dates in kwargs["sample_date_sets"]:
                times = np.array([np.datetime64(dt.date(), "D") for dt in tuple(sample_dates)[:2]])
                coords = {
                    "time": times,
                    "band": [band.name for band in source_bands],
                    "y": [0],
                    "x": [0],
                }
                data = np.array(
                    [
                        [0.08, 0.12, 0.42, 0.30, 0.22],
                        [0.09, 0.13, 0.44, 0.32, 0.24],
                    ],
                    dtype=np.float32,
                ).reshape(2, 5, 1, 1)
                unc = np.full_like(data, 0.01)
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    seen: dict[str, object] = {"calls": 0}

    class _SwitchingMapper:
        def __init__(self, _source_bands, target_bands, *, spectral_library, k_neighbors):  # noqa: ANN001
            del spectral_library, k_neighbors
            self._target_band_names = [band.name for band in target_bands]

        def map(self, reflectance, *, source_uncertainty=None):  # noqa: ANN001
            seen["calls"] += 1
            if seen["calls"] == 1:
                seen["bands"] = tuple(reflectance.coords["band"].values.tolist())
                seen["dims"] = tuple(reflectance.dims)
                seen["time"] = int(reflectance.sizes.get("time", 0))
                seen["has_unc"] = source_uncertainty is not None
            if "time" in reflectance.dims:
                mapped = xr.DataArray(
                    np.full(
                        (
                            int(reflectance.sizes["time"]),
                            len(self._target_band_names),
                            int(reflectance.sizes["y"]),
                            int(reflectance.sizes["x"]),
                        ),
                        0.2,
                        dtype=np.float32,
                    ),
                    dims=["time", "band", "y", "x"],
                    coords={
                        "time": reflectance.coords["time"],
                        "band": self._target_band_names,
                        "y": reflectance.coords["y"],
                        "x": reflectance.coords["x"],
                    },
                )
            else:
                mapped = xr.DataArray(
                    np.full(
                        (
                            len(self._target_band_names),
                            int(reflectance.sizes["y"]),
                            int(reflectance.sizes["x"]),
                        ),
                        0.2,
                        dtype=np.float32,
                    ),
                    dims=["band", "y", "x"],
                    coords={
                        "band": self._target_band_names,
                        "y": reflectance.coords["y"],
                        "x": reflectance.coords["x"],
                    },
                )
            mapped_unc = xr.full_like(mapped, 0.03)
            mapped_fit = xr.zeros_like(mapped_unc.mean(dim="band", skipna=True), dtype=np.float32)
            return mapped, mapped_unc, mapped_fit

    monkeypatch.setattr("siac.algorithms.surface.swir_refine.SpectralMapper", _SwitchingMapper)

    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]
    build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=_TwoSampleProvider(),
        resolution=500.0,
        geometry=_geometry((1, 1)),
        visible_bands=visible_bands,
        query_bands=query_bands,
        spectral_library=_spectral_library(),
    )

    assert seen["bands"] == ("Band3", "Band4", "Band2", "Band6", "Band7")
    assert seen["dims"] == ("time", "band", "y", "x")
    assert seen["time"] == 15
    assert seen["has_unc"] is True
    assert seen["calls"] == 1


def test_normalize_monthly_composites_uses_area_when_target_grid_is_coarser(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    methods: list[str] = []
    sensor_config = _sensor_config()
    band_names = [band.name for band in sensor_config.bands]
    reflectance = xr.DataArray(
        np.full((len(band_names), 2, 2), 0.2, dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": band_names, "y": [1.5, 0.5], "x": [0.5, 1.5]},
    )
    quality = xr.DataArray(
        np.full((2, 2), 0.03, dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [1.5, 0.5], "x": [0.5, 1.5]},
    )
    sample_index = xr.DataArray(
        np.arange(4, dtype=np.int16).reshape(2, 2),
        dims=["y", "x"],
        coords={"y": [1.5, 0.5], "x": [0.5, 1.5]},
    )
    composite = MonthlyBestPixelComposite(
        reflectance=reflectance,
        quality=quality,
        sample_index=sample_index,
        year=2024,
        month=7,
    )

    def _fake_resample_da(da, target_shape, method="bilinear", *, template=None):  # noqa: ANN001
        methods.append(method)
        if da.ndim == 3 and "band" in da.dims:
            band_values = da.coords["band"].values
            return xr.DataArray(
                np.full((len(band_values), target_shape[0], target_shape[1]), 0.2, dtype=np.float32),
                dims=["band", "y", "x"],
                coords={"band": band_values, "y": template.coords["y"], "x": template.coords["x"]},
            )
        return xr.DataArray(
            np.full(target_shape, float(np.asarray(da.values).mean()), dtype=np.float32),
            dims=["y", "x"],
            coords={"y": template.coords["y"], "x": template.coords["x"]},
        )

    monkeypatch.setattr(swir_refine_mod, "_resample_da", _fake_resample_da)

    normalized = swir_refine_mod._normalize_monthly_composites_to_target_basis(
        (composite,),
        geometry=_geometry((1, 1)),
        target_bands=list(sensor_config.bands),
        spectral_mapper=None,
    )

    assert len(normalized) == 1
    assert normalized[0].reflectance.shape == (len(band_names), 1, 1)
    assert methods.count("area") == len(band_names) + 1
    assert methods.count("nearest") == 1


def test_normalize_monthly_kernel_composites_uses_area_when_target_grid_is_coarser(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)
    methods: list[str] = []
    sensor_config = _sensor_config()
    band_names = [band.name for band in sensor_config.bands]
    cube = xr.DataArray(
        np.full((len(band_names), 2, 2), 0.2, dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": band_names, "y": [1.5, 0.5], "x": [0.5, 1.5]},
    )
    composite = MonthlyKernelWeightComposite(
        kernels=BRDFKernelWeights(
            f0=cube,
            f1=xr.zeros_like(cube),
            f2=xr.zeros_like(cube),
            f0_unc=xr.full_like(cube, 0.01),
            f1_unc=xr.full_like(cube, 0.01),
            f2_unc=xr.full_like(cube, 0.01),
        ),
        quality=xr.DataArray(
            np.full((2, 2), 0.03, dtype=np.float32),
            dims=["y", "x"],
            coords={"y": [1.5, 0.5], "x": [0.5, 1.5]},
        ),
        sample_index=xr.DataArray(
            np.arange(4, dtype=np.int16).reshape(2, 2),
            dims=["y", "x"],
            coords={"y": [1.5, 0.5], "x": [0.5, 1.5]},
        ),
        year=2024,
        month=7,
    )

    def _fake_resample_da(da, target_shape, method="bilinear", *, template=None):  # noqa: ANN001
        methods.append(method)
        if da.ndim == 3 and "band" in da.dims:
            band_values = da.coords["band"].values
            return xr.DataArray(
                np.full((len(band_values), target_shape[0], target_shape[1]), 0.2, dtype=np.float32),
                dims=["band", "y", "x"],
                coords={"band": band_values, "y": template.coords["y"], "x": template.coords["x"]},
            )
        return xr.DataArray(
            np.full(target_shape, float(np.asarray(da.values).mean()), dtype=np.float32),
            dims=["y", "x"],
            coords={"y": template.coords["y"], "x": template.coords["x"]},
        )

    monkeypatch.setattr(swir_refine_mod, "_resample_da", _fake_resample_da)

    normalized = swir_refine_mod._normalize_monthly_composites_to_target_basis(
        (composite,),
        geometry=_geometry((1, 1)),
        target_bands=list(sensor_config.bands),
        spectral_mapper=None,
    )

    assert len(normalized) == 1
    assert normalized[0].reflectance.shape == (len(band_names), 1, 1)
    assert methods.count("area") == (len(band_names) * 6) + 1
    assert methods.count("nearest") == 1


def test_normalize_monthly_composites_folds_source_fit_rmse_into_quality() -> None:
    source_bands = (
        SensorBand("Band3", 469.0, 20.0, 500.0, 0),
        SensorBand("Band4", 555.0, 20.0, 500.0, 1),
        SensorBand("Band2", 858.5, 35.0, 500.0, 2),
        SensorBand("Band6", 1640.0, 24.0, 500.0, 3),
        SensorBand("Band7", 2130.0, 50.0, 500.0, 4),
    )
    target_bands = list(_sensor_config().bands)
    composite = MonthlyBestPixelComposite(
        reflectance=xr.DataArray(
            np.full((len(source_bands), 1, 1), 0.2, dtype=np.float32),
            dims=["band", "y", "x"],
            coords={"band": [band.name for band in source_bands], "y": [0], "x": [0]},
        ),
        quality=xr.DataArray(
            np.full((1, 1), 0.03, dtype=np.float32),
            dims=["y", "x"],
            coords={"y": [0], "x": [0]},
        ),
        sample_index=xr.DataArray(
            np.zeros((1, 1), dtype=np.int16),
            dims=["y", "x"],
            coords={"y": [0], "x": [0]},
        ),
        year=2024,
        month=7,
    )

    class _FitAwareMapper:
        def map(self, reflectance, *, source_uncertainty=None):  # noqa: ANN001
            del source_uncertainty
            mapped = xr.DataArray(
                np.full(
                    (
                        int(reflectance.sizes.get("time", 1)),
                        len(target_bands),
                        int(reflectance.sizes["y"]),
                        int(reflectance.sizes["x"]),
                    ),
                    0.2,
                    dtype=np.float32,
                ),
                dims=["time", "band", "y", "x"] if "time" in reflectance.dims else ["band", "y", "x"],
                coords=(
                    {
                        "time": reflectance.coords["time"],
                        "band": [band.name for band in target_bands],
                        "y": reflectance.coords["y"],
                        "x": reflectance.coords["x"],
                    }
                    if "time" in reflectance.dims
                    else {
                        "band": [band.name for band in target_bands],
                        "y": reflectance.coords["y"],
                        "x": reflectance.coords["x"],
                    }
                ),
            )
            mapped_unc = xr.zeros_like(mapped, dtype=np.float32)
            fit = xr.DataArray(
                np.full(tuple(reflectance.sizes[dim] for dim in reflectance.dims if dim != "band"), 0.04, dtype=np.float32),
                dims=tuple(dim for dim in reflectance.dims if dim != "band"),
                coords={dim: reflectance.coords[dim] for dim in reflectance.dims if dim != "band"},
            )
            return mapped, mapped_unc, fit

    normalized = swir_refine_mod._normalize_monthly_composites_to_target_basis(
        (composite,),
        geometry=_geometry((1, 1)),
        target_bands=target_bands,
        spectral_mapper=_FitAwareMapper(),
    )

    assert len(normalized) == 1
    assert float(normalized[0].quality.values[0, 0]) == pytest.approx(0.05, rel=1e-6)


def test_build_monthly_surface_prior_database_reuses_cached_spectral_mapping_for_identical_months(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)

    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((1, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((1, 1)),
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 1.0),
    )

    source_bands = (
        SensorBand("Band3", 469.0, 20.0, 500.0, 0),
        SensorBand("Band4", 555.0, 20.0, 500.0, 1),
        SensorBand("Band2", 858.5, 35.0, 500.0, 2),
        SensorBand("Band6", 1640.0, 24.0, 500.0, 3),
        SensorBand("Band7", 2130.0, 50.0, 500.0, 4),
    )

    class _ConstantSourceBRDFProvider:
        def __init__(self) -> None:
            self.source_bands = source_bands

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            outputs = []
            for sample_dates in kwargs["sample_date_sets"]:
                sample_dates = tuple(sample_dates)
                coords = {
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                    "band": [band.name for band in source_bands],
                    "y": [0],
                    "x": [0],
                }
                base = np.array([0.08, 0.12, 0.42, 0.30, 0.22], dtype=np.float32).reshape(1, 5, 1, 1)
                data = np.repeat(base, len(sample_dates), axis=0)
                unc = np.full_like(data, 0.02)
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    calls = {"map": 0, "band_counts": [], "time_counts": []}

    class _CountingMapper:
        def __init__(self, _source_bands, target_bands, *, spectral_library, k_neighbors):  # noqa: ANN001
            del spectral_library, k_neighbors
            self._target_band_names = [band.name for band in target_bands]

        def map(self, reflectance, *, source_uncertainty=None):  # noqa: ANN001
            del source_uncertainty
            calls["map"] += 1
            calls["band_counts"].append(int(reflectance.sizes["band"]))
            calls["time_counts"].append(int(reflectance.sizes.get("time", 0)))
            if "time" in reflectance.dims:
                mapped = xr.DataArray(
                    np.full(
                        (
                            int(reflectance.sizes["time"]),
                            len(self._target_band_names),
                            int(reflectance.sizes["y"]),
                            int(reflectance.sizes["x"]),
                        ),
                        0.2,
                        dtype=np.float32,
                    ),
                    dims=["time", "band", "y", "x"],
                    coords={
                        "time": reflectance.coords["time"],
                        "band": self._target_band_names,
                        "y": reflectance.coords["y"],
                        "x": reflectance.coords["x"],
                    },
                )
            else:
                mapped = xr.DataArray(
                    np.full(
                        (
                            len(self._target_band_names),
                            int(reflectance.sizes["y"]),
                            int(reflectance.sizes["x"]),
                        ),
                        0.2,
                        dtype=np.float32,
                    ),
                    dims=["band", "y", "x"],
                    coords={
                        "band": self._target_band_names,
                        "y": reflectance.coords["y"],
                        "x": reflectance.coords["x"],
                    },
                )
            mapped_unc = xr.full_like(mapped, 0.01)
            mapped_fit = xr.zeros_like(mapped_unc.mean(dim="band", skipna=True), dtype=np.float32)
            return mapped, mapped_unc, mapped_fit

    monkeypatch.setattr("siac.algorithms.surface.swir_refine.SpectralMapper", _CountingMapper)

    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]
    database = build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=_ConstantSourceBRDFProvider(),
        resolution=500.0,
        geometry=_geometry((1, 1)),
        visible_bands=visible_bands,
        query_bands=query_bands,
        spectral_library=_spectral_library(),
    )

    assert calls["map"] == 1
    assert calls["band_counts"] == [5]
    assert calls["time_counts"] == [15]
    assert database.entries_features.shape[0] == 15


def test_build_monthly_surface_prior_database_accepts_custom_sequence_source_bands(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)

    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((1, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((1, 1)),
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 1.0),
    )
    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]

    class _BandSequence:
        def __init__(self, bands: tuple[SensorBand, ...]) -> None:
            self._bands = bands

        def __getitem__(self, index: int) -> SensorBand:
            return self._bands[index]

        def __len__(self) -> int:
            return len(self._bands)

    class _SequenceBRDFProvider:
        def __init__(self) -> None:
            self.source_bands = _BandSequence((*visible_bands, *query_bands))

        def get_temporal_brdf_parameters_batch(self, **kwargs):
            outputs = []
            coords = {
                "band": [band.name for band in self.source_bands],
                "y": [0],
                "x": [0],
            }
            base = np.array([0.08, 0.12, 0.42, 0.30, 0.22], dtype=np.float32).reshape(1, 5, 1, 1)
            for sample_dates in kwargs["sample_date_sets"]:
                sample_dates = tuple(sample_dates)
                time_coords = {
                    **coords,
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                }
                data = np.repeat(base, len(sample_dates), axis=0)
                unc = np.full_like(data, 0.02)
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=time_coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=time_coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=time_coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=time_coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=time_coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=time_coords),
                    )
                )
            return outputs

    database = build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=_SequenceBRDFProvider(),
        resolution=500.0,
        geometry=_geometry((1, 1)),
        visible_bands=visible_bands,
        query_bands=query_bands,
    )

    assert database.query_band_names == ("B08", "B11", "B12")
    assert database.visible_band_names == ("B02", "B03")


def test_build_monthly_surface_prior_database_falls_back_to_target_bands_when_provider_has_no_source_bands(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_kernel_to_reflectance(monkeypatch)

    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((1, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((1, 1)),
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 1.0),
    )
    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]
    batch_band_calls: list[list[str]] = []

    class _LegacyBRDFProvider:
        def get_temporal_brdf_parameters_batch(self, **kwargs):  # noqa: ANN003
            bands = [band.name for band in kwargs["bands"]]
            batch_band_calls.append(bands)
            outputs = []
            for sample_dates in kwargs["sample_date_sets"]:
                coords = {
                    "time": np.array([np.datetime64(dt.date(), "D") for dt in sample_dates]),
                    "band": bands,
                    "y": [0],
                    "x": [0],
                }
                data = np.full((len(sample_dates), len(bands), 1, 1), 0.2, dtype=np.float32)
                unc = np.full_like(data, 0.03, dtype=np.float32)
                outputs.append(
                    BRDFKernelWeights(
                        f0=xr.DataArray(data, dims=["time", "band", "y", "x"], coords=coords),
                        f1=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f2=xr.DataArray(np.zeros_like(data), dims=["time", "band", "y", "x"], coords=coords),
                        f0_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f1_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                        f2_unc=xr.DataArray(unc, dims=["time", "band", "y", "x"], coords=coords),
                    )
                )
            return outputs

    database = build_monthly_surface_prior_database(
        observation=obs,
        brdf_provider=_LegacyBRDFProvider(),
        resolution=500.0,
        geometry=_geometry((1, 1)),
        visible_bands=visible_bands,
        query_bands=query_bands,
    )

    assert batch_band_calls
    assert batch_band_calls[0] == ["B02", "B03", "B08", "B11", "B12"]
    assert database.visible_band_names == ("B02", "B03")
    assert database.query_band_names == ("B08", "B11", "B12")


def test_build_monthly_surface_prior_database_requires_batch_method() -> None:
    obs = ObservationBundle(
        toa=xr.Dataset(
            {
                "B02": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B03": xr.DataArray(np.full((1, 1), 0.1, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
                "B12": xr.DataArray(np.full((1, 1), 0.2, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        geometry=_geometry((1, 1)),
        cloud_mask=xr.DataArray(np.zeros((1, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 1.0),
    )
    sensor_config = _sensor_config()
    visible_bands = [sensor_config.get_band("B02"), sensor_config.get_band("B03")]
    query_bands = [sensor_config.get_band("B08"), sensor_config.get_band("B11"), sensor_config.get_band("B12")]

    with pytest.raises(TypeError, match="get_temporal_brdf_parameters_batch"):
        build_monthly_surface_prior_database(
            observation=obs,
            brdf_provider=object(),
            resolution=500.0,
            geometry=_geometry((1, 1)),
            visible_bands=visible_bands,
            query_bands=query_bands,
        )
