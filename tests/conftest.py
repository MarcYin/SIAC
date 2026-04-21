"""
Pytest configuration and fixtures for SIAC tests.
"""

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

# =============================================================================
# Path fixtures
# =============================================================================

@pytest.fixture
def fixtures_dir() -> Path:
    """Path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def sample_data_dir(fixtures_dir: Path) -> Path:
    """Path to sample data directory."""
    return fixtures_dir / "sample_data"


# =============================================================================
# Geometry fixtures
# =============================================================================

@pytest.fixture
def sample_geometry() -> dict:
    """Sample viewing/solar geometry arrays."""
    shape = (100, 100)

    # Create realistic angle ranges
    sza = np.full(shape, 30.0) * np.pi / 180  # 30 degrees solar zenith
    saa = np.full(shape, 150.0) * np.pi / 180  # 150 degrees solar azimuth
    vza = np.full(shape, 5.0) * np.pi / 180  # 5 degrees view zenith
    vaa = np.full(shape, 100.0) * np.pi / 180  # 100 degrees view azimuth

    # Add some spatial variation
    x = np.linspace(0, 1, shape[1])
    y = np.linspace(0, 1, shape[0])
    xx, yy = np.meshgrid(x, y)

    vza = vza + (xx * 0.05)  # Small VZA variation across swath

    return {
        "sza": sza,
        "saa": saa,
        "vza": vza,
        "vaa": vaa,
    }


@pytest.fixture
def sample_geometry_xr(sample_geometry: dict) -> dict:
    """Sample geometry as xarray DataArrays."""
    coords = {
        "y": np.arange(sample_geometry["sza"].shape[0]),
        "x": np.arange(sample_geometry["sza"].shape[1]),
    }

    return {
        key: xr.DataArray(val, dims=["y", "x"], coords=coords)
        for key, val in sample_geometry.items()
    }


# =============================================================================
# Atmospheric state fixtures
# =============================================================================

@pytest.fixture
def sample_atmo_state() -> dict:
    """Sample atmospheric state arrays."""
    shape = (100, 100)

    return {
        "aot": np.full(shape, 0.15),  # AOT at 550nm
        "tcwv": np.full(shape, 2.5),  # Total column water vapor (g/cm^2)
        "tco3": np.full(shape, 0.3),  # Total column ozone (atm-cm)
        "aot_unc": np.full(shape, 0.05),
        "tcwv_unc": np.full(shape, 0.3),
        "tco3_unc": np.full(shape, 0.01),
        "elevation": np.full(shape, 0.1),  # Elevation (km)
    }


# =============================================================================
# BRDF fixtures
# =============================================================================

@pytest.fixture
def sample_brdf_weights() -> dict:
    """Sample BRDF kernel weights."""
    shape = (7, 100, 100)  # 7 MODIS bands

    # Typical vegetation BRDF parameters
    return {
        "f0": np.random.default_rng().uniform(0.02, 0.3, shape),  # Isotropic
        "f1": np.random.default_rng().uniform(0.01, 0.1, shape),  # Volumetric
        "f2": np.random.default_rng().uniform(0.005, 0.05, shape),  # Geometric
        "f0_unc": np.full(shape, 0.01),
        "f1_unc": np.full(shape, 0.005),
        "f2_unc": np.full(shape, 0.002),
    }


# =============================================================================
# TOA reflectance fixtures
# =============================================================================

@pytest.fixture
def sample_toa() -> xr.Dataset:
    """Sample TOA reflectance dataset."""
    bands = ["B02", "B03", "B04", "B08", "B11", "B12"]
    shape = (100, 100)

    data_vars = {}
    for band in bands:
        # Realistic TOA reflectance values
        toa = np.random.default_rng().uniform(0.05, 0.4, shape).astype(np.float32)
        data_vars[band] = xr.DataArray(
            toa,
            dims=["y", "x"],
            attrs={"units": "reflectance", "wavelength_nm": 490 if band == "B02" else 560},
        )

    ds = xr.Dataset(data_vars)
    ds.attrs["observation_time"] = datetime(2023, 7, 15, 10, 30, 0).isoformat()
    ds.attrs["sensor"] = "MSI"
    ds.attrs["satellite"] = "S2A"

    return ds


# =============================================================================
# Emulator fixtures
# =============================================================================

@pytest.fixture
def sample_emulator_weights(tmp_path: Path) -> Path:
    """Create sample emulator weights file."""
    rng = np.random.default_rng(42)
    hidden = 64
    input_dim = 7
    output_dim = 1

    # Random weights for testing
    w1 = rng.standard_normal((input_dim, hidden)).astype(np.float32) * 0.1
    b1 = np.zeros(hidden, dtype=np.float32)
    w2 = rng.standard_normal((hidden, hidden)).astype(np.float32) * 0.1
    b2 = np.zeros(hidden, dtype=np.float32)
    w3 = rng.standard_normal((hidden, output_dim)).astype(np.float32) * 0.1
    b3 = np.zeros(output_dim, dtype=np.float32)

    path = tmp_path / "test_emulator.npz"
    np.savez(
        path,
        w1=w1, b1=b1,
        w2=w2, b2=b2,
        w3=w3, b3=b3,
        input_scale=np.array([1.0, 0.0]),
        output_scale=np.array([1.0, 0.0]),
    )

    return path


# =============================================================================
# Configuration fixtures
# =============================================================================

@pytest.fixture
def sample_config_dict() -> dict:
    """Sample configuration dictionary."""
    return {
        "sensor": "auto",
        "atmo_prior": {
            "provider": "cams",
        },
        "brdf": {
            "provider": "mcd43",
            "temporal_window": 16,
        },
        "surface_prior": {
            "method": "kernel_model",
            "psf_sigma_x": 29.75,
            "psf_sigma_y": 39.0,
        },
        "rt_model": {
            "backend": "emulator",
        },
        "solver": {
            "aot_gamma": 10.0,
            "tcwv_gamma": 5.0,
            "aerosol_resolution": 120.0,
        },
        "output": {
            "format": "cog",
            "include_uncertainty": True,
        },
    }


# =============================================================================
# Mock RT model (shared across tests)
# =============================================================================

class MockRTModel:
    """Mock RT model that satisfies RTModelBackend protocol."""

    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):
        """Return mock coefficients."""
        from siac.runtime import RTCoefficients
        shape = geometry.sza.shape

        xap = xr.DataArray(np.full(shape, 0.95), dims=["y", "x"])
        xbp = xr.DataArray(np.full(shape, 0.02), dims=["y", "x"])
        xcp = xr.DataArray(np.full(shape, 0.1), dims=["y", "x"])

        d_xap = d_xbp = d_xcp = None
        if compute_jacobian:
            d_xap = xr.concat(
                [xr.DataArray(np.full(shape, -0.5), dims=["y", "x"]),
                 xr.DataArray(np.full(shape, -0.01), dims=["y", "x"])],
                dim="param",
            ).assign_coords(param=["aot", "tcwv"])
            d_xbp = xr.concat(
                [xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
                 xr.DataArray(np.full(shape, 0.005), dims=["y", "x"])],
                dim="param",
            ).assign_coords(param=["aot", "tcwv"])
            d_xcp = xr.concat(
                [xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
                 xr.DataArray(np.full(shape, 0.002), dims=["y", "x"])],
                dim="param",
            ).assign_coords(param=["aot", "tcwv"])

        return RTCoefficients(xap=xap, xbp=xbp, xcp=xcp,
                              d_xap=d_xap, d_xbp=d_xbp, d_xcp=d_xcp)

    def supports_jacobian(self):
        return True

    @property
    def backend_name(self):
        return "mock"

    def is_available_for_sensor(self, sensor_id, satellite_id):
        return True


@pytest.fixture
def mock_rt_model():
    """Shared mock RT model fixture."""
    return MockRTModel()


# =============================================================================
# Pipeline contract fixtures (new)
# =============================================================================

@pytest.fixture
def mock_sensor_config():
    """3-band sensor config (Blue, Green, Red) for pipeline contract tests."""
    from siac.domain import SensorBand, SensorConfig
    return SensorConfig(
        sensor_id="MOCK",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
        ),
        default_ref_scale=1.0 / 10000.0,
        default_ref_offset=0.0,
    )


PIPELINE_SHAPE = (32, 32)


@pytest.fixture
def mock_geometry():
    """GeometryAngles at pipeline test resolution."""
    from siac.runtime import GeometryAngles
    shape = PIPELINE_SHAPE
    return GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )


@pytest.fixture
def mock_observation_bundle(mock_sensor_config, mock_geometry):
    """Complete valid ObservationBundle with 32x32 synthetic TOA."""
    from siac.runtime import ObservationBundle
    shape = PIPELINE_SHAPE
    toa_ds = xr.Dataset({
        "B02": xr.DataArray(np.random.RandomState(42).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
        "B03": xr.DataArray(np.random.RandomState(43).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
        "B04": xr.DataArray(np.random.RandomState(44).uniform(0.05, 0.3, shape).astype(np.float32), dims=["y", "x"]),
    })
    cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
    return ObservationBundle(
        toa=toa_ds,
        geometry=mock_geometry,
        cloud_mask=cloud_mask,
        sensor_config=mock_sensor_config,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30, 0)},
        crs="EPSG:32632",
        bounds=(300000.0, 5500000.0, 309780.0, 5600040.0),
    )


@pytest.fixture
def mock_atmospheric_state():
    """Spatially uniform AtmosphericState at 32x32."""
    from siac.runtime import AtmosphericState
    shape = PIPELINE_SHAPE
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
    )


@pytest.fixture
def mock_surface_prior():
    """Uniform SurfacePrior at 32x32."""
    from siac.runtime import BRDFKernelWeights, SurfacePrior
    shape = PIPELINE_SHAPE
    brdf = BRDFKernelWeights(
        f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
        f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
    )
    return SurfacePrior(
        boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
        boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        kernels=brdf,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
    )


@pytest.fixture
def mock_solved_atmosphere(mock_atmospheric_state):
    """Dummy solved output (converged, 5 iterations)."""
    from siac.runtime import SolvedAtmosphere
    return SolvedAtmosphere(
        atmo_state=mock_atmospheric_state,
        aot=mock_atmospheric_state.aot,
        tcwv=mock_atmospheric_state.tcwv,
        aot_unc=mock_atmospheric_state.aot_unc,
        tcwv_unc=mock_atmospheric_state.tcwv_unc,
        cost_final=0.001,
        n_iterations=5,
        converged=True,
    )


@pytest.fixture
def mock_solver_input_bundle(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_rt_model,
):
    """Pre-assembled SolverInputBundle (skips M4)."""
    from siac.runtime import SolverInputBundle
    obs = mock_observation_bundle
    bands = obs.sensor_config.default_aerosol_solver_bands()
    band_names = [b.name for b in bands]
    toa_arrays = [obs.toa[bn] for bn in band_names if bn in obs.toa.data_vars]
    toa_da = xr.concat(toa_arrays, dim="band") if toa_arrays else obs.toa[list(obs.toa.data_vars)[0]].expand_dims("band")
    toa_da = toa_da.assign_coords(band=[data.name for data in toa_arrays] or [list(obs.toa.data_vars)[0]])
    return SolverInputBundle(
        toa=toa_da,
        geometry=obs.geometry,
        cloud_mask=obs.cloud_mask,
        sensor_config=obs.sensor_config,
        bands=bands,
        atmo_prior=mock_atmospheric_state,
        surface_prior=mock_surface_prior,
        rt_model=mock_rt_model,
        aux_resolution_m=500.0,
        aerosol_resolution_m=1000.0,
    )


# =============================================================================
# Mock module callables for pipeline tests
# =============================================================================

@pytest.fixture
def mock_preprocessor(mock_observation_bundle):
    """(Path, AOI) -> ObservationBundle. Returns fixed bundle."""
    def _preprocess(input_path, aoi=None):
        return mock_observation_bundle
    return _preprocess


@pytest.fixture
def mock_atmo_provider(mock_atmospheric_state):
    """(bounds, crs, time, res) -> AtmosphericState."""
    def _get_prior(bounds, crs, obs_time, resolution):
        return mock_atmospheric_state
    return _get_prior


@pytest.fixture
def mock_surface_prior_provider(mock_surface_prior):
    """(observation, atmo_prior|None, rt_model, res) -> SurfacePrior."""
    def _get_surface_prior(observation, atmo_prior, rt_model, resolution):
        _ = (observation, atmo_prior, rt_model, resolution)
        return mock_surface_prior
    _get_surface_prior.requires_atmo_prior = False
    return _get_surface_prior


@pytest.fixture
def mock_grid_assembler(mock_solver_input_bundle):
    """Grid assembler that returns a pre-built SolverInputBundle."""
    def _assemble(
        obs,
        atmo,
        surface,
        rt_model,
        aux_resolution_m=500.0,
        aerosol_resolution_m=1000.0,
        **_kwargs,
    ):
        _ = (obs, atmo, surface, rt_model, aux_resolution_m, aerosol_resolution_m)
        return mock_solver_input_bundle
    return _assemble


@pytest.fixture
def mock_solver_fn(mock_solved_atmosphere):
    """Solver that returns a pre-built SolvedAtmosphere."""
    def _solve(inputs, config):
        return mock_solved_atmosphere
    return _solve


@pytest.fixture
def mock_corrector_fn(mock_observation_bundle, mock_solved_atmosphere):
    """Corrector that returns a synthetic CorrectionResult."""
    from siac.runtime import CorrectionDiagnostics, CorrectionResult
    def _correct(obs, solved, rt_model):
        return CorrectionResult(
            boa=obs.toa,
            boa_unc=None,
            aot=solved.aot,
            tcwv=solved.tcwv,
            cloud_mask=obs.cloud_mask,
            diagnostics=CorrectionDiagnostics(processing_time_s=0.01),
        )
    return _correct


# =============================================================================
# Cost function input fixtures
# =============================================================================

@pytest.fixture
def cost_function_inputs():
    """Standard inputs for cost function tests."""
    from siac.domain import SensorBand
    from siac.runtime import AtmosphericState, BRDFKernelWeights, GeometryAngles, SurfacePrior

    shape = (16, 16)

    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )

    atmo_prior = AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
    )

    toa = xr.DataArray(
        np.random.RandomState(42).uniform(0.05, 0.3, (3, *shape)).astype(np.float32),
        dims=["band", "y", "x"],
    )

    brdf_weights = BRDFKernelWeights(
        f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
        f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
    )

    surface_prior = SurfacePrior(
        boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
        boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        kernels=brdf_weights,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
    )

    mask = xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"])

    bands = [
        SensorBand("B02", 490.0, 65.0, 10.0, 0),
        SensorBand("B03", 560.0, 35.0, 10.0, 1),
        SensorBand("B04", 665.0, 30.0, 10.0, 2),
    ]

    return {
        "toa": toa,
        "surface_prior": surface_prior,
        "geometry": geometry,
        "atmo_prior": atmo_prior,
        "bands": bands,
        "mask": mask,
    }


# =============================================================================
# Markers
# =============================================================================

def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line("markers", "slow: marks tests as slow")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line("markers", "regression: marks tests as regression tests")
    config.addinivalue_line("markers", "gee: marks tests requiring Google Earth Engine")
