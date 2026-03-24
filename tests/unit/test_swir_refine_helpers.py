from __future__ import annotations

from datetime import datetime
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.surface import swir_refine
from siac.algorithms.surface.swir_refine import (
    _deduplicate_bands,
    _get_temporal_brdf_parameters_batch,
    _iter_history_months,
    _month_center_datetime,
    _native_observation_resolution,
    _observation_shape,
    _resample_atmo_to_observation_grid,
    _resample_dataset,
    _resolve_provider_source_bands,
    _select_month_mask,
    _weekly_sample_dates,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)
from siac.domain import SensorBand, SensorConfig
from siac.runtime import AtmosphericState, GeometryAngles, ObservationBundle, RTCoefficients


def _sensor_config() -> SensorConfig:
    return SensorConfig(
        sensor_id="MSI",
        satellite_id="S2A",
        bands=(
            SensorBand("B02", 490.0, 65.0, 10.0, 0),
            SensorBand("B03", 560.0, 35.0, 10.0, 1),
            SensorBand("B08", 842.0, 115.0, 10.0, 2),
            SensorBand("B11", 1610.0, 90.0, 20.0, 3),
        ),
    )


def _geometry(shape: tuple[int, int]) -> GeometryAngles:
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"]),
    )


def _observation(shape: tuple[int, int] = (2, 2)) -> ObservationBundle:
    toa = xr.Dataset(
        {
            "B02": xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            "B03": xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            "B08": xr.DataArray(np.full(shape, 0.4, dtype=np.float32), dims=["y", "x"]),
            "B11": xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        }
    )
    return ObservationBundle(
        toa=toa,
        geometry=_geometry(shape),
        cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        sensor_config=_sensor_config(),
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, float(shape[1]), float(shape[0])),
    )


def _atmo(shape: tuple[int, int]) -> AtmosphericState:
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
    )


class _IdentityRTModel:
    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # type: ignore[no-untyped-def]
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


class _StubDatabase:
    query_band_names = ("B08", "B11")
    visible_band_names = ("B02",)

    def __init__(self) -> None:
        self.median_summary = xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=["y", "x"],
            coords={"y": [0, 1], "x": [0, 1]},
        )

    def predict_visible(self, corrected_query: xr.Dataset, *, k_neighbors: int) -> tuple[xr.DataArray, xr.DataArray]:
        assert k_neighbors == 4
        assert set(corrected_query.data_vars) == {"B08", "B11"}
        coords = {"band": ["B02"], "y": [0, 1], "x": [0, 1]}
        predicted = xr.DataArray(
            np.array([[[0.2, 0.3], [np.nan, 0.4]]], dtype=np.float32),
            dims=["band", "y", "x"],
            coords=coords,
        )
        uncertainty = xr.DataArray(
            np.array([[[0.05, 0.06], [0.07, np.nan]]], dtype=np.float32),
            dims=["band", "y", "x"],
            coords=coords,
        )
        return predicted, uncertainty


def test_query_surface_prior_validates_band_ordering() -> None:
    database = _StubDatabase()
    obs = _observation()
    atmo = _atmo((2, 2))

    with pytest.raises(ValueError, match="query-band ordering"):
        query_surface_prior_from_monthly_database(
            observation=obs,
            atmo_prior=atmo,
            rt_model=object(),
            database=database,
            query_band_names=("B11", "B08"),
        )

    with pytest.raises(ValueError, match="visible-band ordering"):
        query_surface_prior_from_monthly_database(
            observation=obs,
            atmo_prior=atmo,
            rt_model=object(),
            database=database,
            visible_band_names=("B03",),
        )


def test_query_surface_prior_combines_cloud_mask_resampling_and_invalid_pixels(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    correction = SimpleNamespace(
        boa=xr.Dataset(
            {
                "B08": xr.DataArray(np.full((1, 1), 0.4, dtype=np.float32), dims=["y", "x"]),
                "B11": xr.DataArray(np.full((1, 1), 0.3, dtype=np.float32), dims=["y", "x"]),
            }
        ),
        cloud_mask=xr.DataArray(
            np.array([[[False]], [[True]]], dtype=bool),
            dims=["band", "y", "x"],
            coords={"band": ["B08", "B11"], "y": [0], "x": [0]},
        ),
    )
    def fake_correct(self, toa, geometry, aligned_atmo, cloud_mask=None):  # type: ignore[no-untyped-def]
        return correction

    def fake_resample_cloud_mask(_cloud_mask, _target_shape):  # type: ignore[no-untyped-def]
        return xr.DataArray(
            np.array([[False, True], [False, False]], dtype=bool),
            dims=["y", "x"],
            coords={"y": [0, 1], "x": [0, 1]},
        )

    monkeypatch.setattr(swir_refine.AtmosphericCorrector, "correct", fake_correct)
    monkeypatch.setattr(swir_refine, "_resample_cloud_mask", fake_resample_cloud_mask)

    prior = query_surface_prior_from_monthly_database(
        observation=_observation(),
        atmo_prior=_atmo((2, 2)),
        rt_model=_IdentityRTModel(),
        database=_StubDatabase(),
        query_band_names=("B08", "B11"),
        visible_band_names=("B02",),
        k_neighbors=4,
    )

    assert prior.boa.shape == (1, 2, 2)
    assert prior.mask.shape == (2, 2)
    assert prior.mask.values.tolist() == [[True, False], [False, False]]


def test_query_surface_prior_corrects_only_resampled_query_bands(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, object] = {}

    def fake_correct(self, toa, geometry, aligned_atmo, cloud_mask=None):  # type: ignore[no-untyped-def]
        seen["bands"] = tuple(toa.data_vars)
        seen["toa_shape"] = toa["B08"].shape
        seen["geometry_shape"] = geometry.sza.shape
        seen["atmo_shape"] = aligned_atmo.aot.shape
        seen["cloud_shape"] = None if cloud_mask is None else cloud_mask.shape
        return SimpleNamespace(boa=toa, cloud_mask=cloud_mask)

    monkeypatch.setattr(swir_refine.AtmosphericCorrector, "correct", fake_correct)

    prior = query_surface_prior_from_monthly_database(
        observation=_observation((4, 4)),
        atmo_prior=_atmo((2, 2)),
        rt_model=_IdentityRTModel(),
        database=_StubDatabase(),
        query_band_names=("B08", "B11"),
        visible_band_names=("B02",),
        k_neighbors=4,
    )

    assert seen["bands"] == ("B08", "B11")
    assert seen["toa_shape"] == (2, 2)
    assert seen["geometry_shape"] == (2, 2)
    assert seen["atmo_shape"] == (2, 2)
    assert seen["cloud_shape"] == (2, 2)
    assert prior.boa.shape == (1, 2, 2)


def test_query_surface_prior_uses_database_grid_for_knn_query(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, object] = {}

    class _MismatchedGridDatabase:
        query_band_names = ("B08", "B11")
        visible_band_names = ("B02",)

        def __init__(self) -> None:
            self.median_summary = xr.DataArray(
                np.zeros((2, 4, 4), dtype=np.float32),
                dims=["feature", "y", "x"],
                coords={"feature": ["median_query_0", "median_query_1"], "y": np.arange(4), "x": np.arange(4)},
            )

        def predict_visible(self, corrected_query: xr.Dataset, *, k_neighbors: int) -> tuple[xr.DataArray, xr.DataArray]:
            seen["knn_shape"] = corrected_query["B08"].shape
            assert k_neighbors == 4
            coords = {"band": ["B02"], "y": np.arange(4), "x": np.arange(4)}
            values = np.full((1, 4, 4), 0.2, dtype=np.float32)
            return (
                xr.DataArray(values, dims=["band", "y", "x"], coords=coords),
                xr.DataArray(np.full_like(values, 0.05), dims=["band", "y", "x"], coords=coords),
            )

    def fake_correct(self, toa, geometry, aligned_atmo, cloud_mask=None):  # type: ignore[no-untyped-def]
        seen["toa_shape"] = toa["B08"].shape
        seen["geometry_shape"] = geometry.sza.shape
        seen["atmo_shape"] = aligned_atmo.aot.shape
        return SimpleNamespace(boa=toa, cloud_mask=cloud_mask)

    monkeypatch.setattr(swir_refine.AtmosphericCorrector, "correct", fake_correct)

    prior = query_surface_prior_from_monthly_database(
        observation=_observation((8, 8)),
        atmo_prior=_atmo((2, 2)),
        rt_model=_IdentityRTModel(),
        database=_MismatchedGridDatabase(),
        query_band_names=("B08", "B11"),
        visible_band_names=("B02",),
        k_neighbors=4,
    )

    assert seen["toa_shape"] == (4, 4)
    assert seen["geometry_shape"] == (4, 4)
    assert seen["atmo_shape"] == (4, 4)
    assert seen["knn_shape"] == (4, 4)
    assert prior.boa.shape == (1, 4, 4)


def test_resample_geometry_for_surface_prior_uses_native_resolution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: list[tuple[tuple[int, int], float, float]] = []

    def fake_target_shape(shape: tuple[int, int], native_resolution: float, resolution: float) -> tuple[int, int]:
        seen.append((shape, native_resolution, resolution))
        return (1, 2)

    def fake_resample(data: xr.DataArray, target_shape: tuple[int, int], method: str) -> xr.DataArray:
        return xr.DataArray(np.full(target_shape, 1.0, dtype=np.float32), dims=["y", "x"])

    monkeypatch.setattr(swir_refine, "_compute_target_shape", fake_target_shape)
    monkeypatch.setattr(swir_refine, "_resample_da", fake_resample)

    out = resample_geometry_for_surface_prior(_observation(), resolution=20.0)

    assert seen == [((2, 2), 10.0, 20.0)]
    assert out.sza.shape == (1, 2)
    assert out.vaa.shape == (1, 2)


def test_history_month_helpers_cover_wraparound_and_month_end_logic() -> None:
    january_months = _iter_history_months(datetime(2024, 1, 15, 10, 30))
    december_months = _iter_history_months(datetime(2024, 12, 15, 10, 30))

    assert january_months[:3] == [(2022, 12), (2023, 1), (2023, 2)]
    assert december_months[:3] == [(2023, 11), (2023, 12), (2024, 1)]
    assert _month_center_datetime(2024, 2, datetime(2024, 7, 1, 9, 45)) == datetime(2024, 2, 15, 9, 45)
    assert [dt.day for dt in _weekly_sample_dates(2023, 2)] == [1, 8, 15, 22, 28]


def test_swir_refine_helper_functions_cover_error_and_fallback_paths() -> None:
    sensor_config = _sensor_config()
    duplicate_bands = [
        sensor_config.get_band("B02"),
        sensor_config.get_band("B02"),
        sensor_config.get_band("B08"),
    ]
    assert [band.name for band in _deduplicate_bands(duplicate_bands)] == ["B02", "B08"]

    month_mask = _select_month_mask(
        np.array(["2024-06-30", "2024-07-02", "2024-07-30"], dtype="datetime64[D]"),
        year=2024,
        month=7,
    )
    assert month_mask.tolist() == [False, True, True]

    obs = _observation()
    assert _native_observation_resolution(obs) == 10.0
    atmo = _atmo((2, 2))
    assert _resample_atmo_to_observation_grid(obs, atmo) is atmo

    no_matching_obs = ObservationBundle(
        toa=xr.Dataset({"OTHER": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])}),
        geometry=_geometry((2, 2)),
        cloud_mask=xr.DataArray(np.zeros((2, 2), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 2.0, 2.0),
    )
    assert _native_observation_resolution(no_matching_obs) == 10.0

    bad_shape_obs = ObservationBundle(
        toa=xr.Dataset({"B02": xr.DataArray(np.ones(3, dtype=np.float32), dims=["y"])}),
        geometry=_geometry((3, 1)),
        cloud_mask=xr.DataArray(np.zeros((3, 1), dtype=bool), dims=["y", "x"]),
        sensor_config=sensor_config,
        metadata={"observation_time": datetime(2024, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(0.0, 0.0, 1.0, 3.0),
    )
    with pytest.raises(ValueError, match="2-D"):
        _observation_shape(bad_shape_obs)

    with pytest.raises(ValueError, match="No query bands"):
        _resample_dataset(
            xr.Dataset({"B02": xr.DataArray(np.ones((2, 2)), dims=["y", "x"])}),
            band_names=("B08",),
            target_shape=(1, 1),
        )

    fallback_bands = _resolve_provider_source_bands(object(), duplicate_bands)
    assert [band.name for band in fallback_bands] == ["B02", "B08"]

    with pytest.raises(TypeError, match="sequence"):
        _resolve_provider_source_bands(SimpleNamespace(source_bands="B02"), duplicate_bands)
    with pytest.raises(TypeError, match="sequence"):
        _resolve_provider_source_bands(SimpleNamespace(source_bands=object()), duplicate_bands)
    with pytest.raises(TypeError, match="SensorBand"):
        _resolve_provider_source_bands(SimpleNamespace(source_bands=[sensor_config.get_band("B02"), "bad"]), duplicate_bands)

    class _Provider:
        def get_temporal_brdf_parameters_batch(self, **kwargs):  # type: ignore[no-untyped-def]
            return ["ok", kwargs["target_resolution"]]

    assert _get_temporal_brdf_parameters_batch(
        _Provider(),
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        obs_times=[datetime(2024, 7, 15)],
        target_resolution=500.0,
        bands=duplicate_bands,
        temporal_windows=[15],
        sample_date_sets=[_weekly_sample_dates(2024, 7)],
    ) == ["ok", 500.0]

    with pytest.raises(TypeError, match="get_temporal_brdf_parameters_batch"):
        _get_temporal_brdf_parameters_batch(
            SimpleNamespace(get_temporal_brdf_parameters_batch=None),
            bounds=(0.0, 0.0, 1.0, 1.0),
            crs="EPSG:4326",
            obs_times=[datetime(2024, 7, 15)],
            target_resolution=500.0,
            bands=duplicate_bands,
            temporal_windows=[15],
            sample_date_sets=[_weekly_sample_dates(2024, 7)],
        )
