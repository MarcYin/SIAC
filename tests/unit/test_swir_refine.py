from __future__ import annotations

from datetime import datetime

import numpy as np
import xarray as xr

from siac.core.types import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    RTCoefficients,
    SensorBand,
    SensorConfig,
)
from siac.priors.surface.brdf_monthly_composite import MonthlyBestPixelComposite
from siac.priors.surface.brdf_monthly_database import build_monthly_composite_database
from siac.priors.surface.swir_refine import (
    _forward_model_monthly_reflectance,
    _weekly_sample_dates,
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
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


def test_weekly_sample_dates_use_weekly_spacing() -> None:
    sample_dates = _weekly_sample_dates(2024, 1)
    assert [dt.day for dt in sample_dates] == [1, 8, 15, 22, 29]


def test_build_monthly_surface_prior_database_requests_weekly_brdf_samples() -> None:
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
    assert len(provider.batch_calls) == 1
    july_call = tuple(provider.batch_calls[0]["sample_date_sets"][1])
    assert [dt.day for dt in july_call] == [1, 8, 15, 22, 29]


def test_forward_model_monthly_reflectance_uses_qa_based_reflectance_uncertainty() -> None:
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

    _reflectance, quality = _forward_model_monthly_reflectance(
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
