"""
Pytest configuration and fixtures for SIAC tests.
"""

import numpy as np
import pytest
import xarray as xr
from pathlib import Path
from datetime import datetime


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
        "f0": np.random.uniform(0.02, 0.3, shape),  # Isotropic
        "f1": np.random.uniform(0.01, 0.1, shape),  # Volumetric
        "f2": np.random.uniform(0.005, 0.05, shape),  # Geometric
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
        toa = np.random.uniform(0.05, 0.4, shape).astype(np.float32)
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
    hidden = 64
    input_dim = 7
    output_dim = 1

    # Random weights for testing
    w1 = np.random.randn(input_dim, hidden).astype(np.float32) * 0.1
    b1 = np.zeros(hidden, dtype=np.float32)
    w2 = np.random.randn(hidden, hidden).astype(np.float32) * 0.1
    b2 = np.zeros(hidden, dtype=np.float32)
    w3 = np.random.randn(hidden, output_dim).astype(np.float32) * 0.1
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
            "aerosol_resolution": 1000.0,
        },
        "output": {
            "format": "cog",
            "include_uncertainty": True,
        },
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
    config.addinivalue_line("markers", "py6s: marks tests requiring Py6S")
