# SIAC v2 - Sensor-Invariant Atmospheric Correction

A modular, extensible framework for atmospheric correction of satellite imagery.

## Features

- **Modular Architecture**: Pluggable components for atmospheric priors, BRDF products, RT models, and sensors
- **Modern Stack**: Built with xarray, rioxarray, Pydantic, and Rust (via PyO3)
- **Extensible**: Easy to add new sensors, data sources, and processing backends
- **Fast**: Rust-accelerated BRDF kernels, PSF convolution, and neural network inference
- **Well-tested**: Comprehensive test suite with unit, integration, and regression tests

## Supported Satellites

- Sentinel-2 A/B (MSI)
- Landsat 8/9 (OLI/OLI-2)
- Extensible to other sensors via LUT or Py6S backends

## Installation

```bash
# From PyPI (when released)
pip install siac

# Development installation
git clone https://github.com/MarcYin/SIAC.git
cd SIAC/siac_v2
pip install -e ".[dev]"

# Build Rust extensions
maturin develop --release
```

## Quick Start

```python
from siac import SIAC
from siac.core.config import SIACConfig

# Load configuration
config = SIACConfig.from_yaml("siac_config.yaml")

# Or use defaults
config = SIACConfig()

# Process Sentinel-2 data
siac = SIAC(config)
result = siac.process("/path/to/S2_SAFE/")

# Access results
print(result.boa)  # BOA reflectance dataset
print(result.aot.mean())  # Mean AOT
```

## Configuration

SIAC uses a hierarchical configuration system. Create a `siac_config.yaml`:

```yaml
sensor: auto  # s2, l8, or auto

atmo_prior:
  provider: cams  # cams, merra2, era5, or user

brdf:
  provider: mcd43  # mcd43, vnp43, mcd19, gee, zarr
  temporal_window: 16

surface_prior:
  method: kernel_model
  psf_sigma_x: 29.75
  psf_sigma_y: 39.0

rt_model:
  backend: emulator  # emulator, lut, or py6s

solver:
  aot_gamma: 10.0
  tcwv_gamma: 5.0
  aerosol_resolution: 1000.0

output:
  format: cog
  include_uncertainty: true
```

## Architecture

```
siac/
├── core/           # Data types, protocols, configuration
├── satellite/      # Sensor-specific preprocessors (S2, L8)
├── priors/
│   ├── atmospheric/  # CAMS, MERRA-2, ERA5 providers
│   ├── brdf/         # MCD43, VNP43, MCD19 providers
│   └── surface/      # Surface prior derivation
├── rt/
│   ├── emulator/     # Neural network emulators (fast)
│   ├── lut/          # Look-up table backend (medium)
│   └── direct/       # Py6S simulation (slow)
├── solver/         # Multi-grid L-BFGS-B optimization
└── correction/     # TOA to BOA conversion
```

```mermaid

					graph LR
    subgraph SOLVER_INPUTS["Aerosol & WV Retrieval Inputs"]
        TOA["TOA reflectance<br/>(6 bands at solver res)"]
        GEOM["GeometryAngles<br/>cos(SZA), cos(VZA), cos(RAA)"]
        ATMO_PRIOR["AtmosphericState prior<br/>AOT₀±σ, TCWV₀±σ, TCO3"]
        SURF_PRIOR["SurfacePrior<br/>BOA_prior±σ per band"]
        CLOUD["Cloud mask<br/>(valid pixel mask)"]
        ELEV["Elevation<br/>(km, from DEM)"]
        RT["RT Model<br/>(emulator or LUT)"]
    end

    subgraph SOLVER_PROCESS["Solver"]
        COST["CostFunction<br/>J_obs + J_prior + J_smooth"]
        LBFGSB["L-BFGS-B<br/>multi-grid"]
    end

    TOA --> COST
    GEOM --> COST
    ATMO_PRIOR --> COST
    SURF_PRIOR --> COST
    CLOUD --> COST
    ELEV --> COST
    RT --> COST
    COST --> LBFGSB

    LBFGSB --> AOT_OUT["AOT_solved ± σ"]
    LBFGSB --> TCWV_OUT["TCWV_solved ± σ"]

    subgraph CORRECTION_INPUTS["Atmospheric Correction Inputs"]
        TOA2["TOA reflectance<br/>(all 13 bands, full res)"]
        GEOM2["GeometryAngles<br/>(full res)"]
        SOLVED["Solved AtmosphericState<br/>AOT_solved, TCWV_solved, TCO3"]
        RT2["RT Model"]
    end

    AOT_OUT --> SOLVED
    TCWV_OUT --> SOLVED

    subgraph CORRECTION_PROCESS["Correction"]
        RTCOEFF["RT coefficients<br/>xap, xbp, xcp per band"]
        BOA_EQ["boa = (xap·toa−xbp)/(1+xcp·y)"]
    end

    TOA2 --> RTCOEFF
    GEOM2 --> RTCOEFF
    SOLVED --> RTCOEFF
    RT2 --> RTCOEFF
    RTCOEFF --> BOA_EQ
    BOA_EQ --> BOA_FINAL["BOA reflectance"]
```

```mermaid
					graph TD
    subgraph INPUT["User Input"]
        Q["S2 Query / Path + Config"]
    end

    Q --> RESOLVE["Resolve AOI & obs_time"]

    subgraph PARALLEL["Parallel Data Fetch & Preprocessing"]
        direction LR
        subgraph G1["Group A: Satellite Data"]
            S2_FETCH["Fetch S2 L1C<br/>(CDSE / GCS / local)"]
            S2_LOAD["Load TOA reflectance<br/>(13 bands, JP2→xarray)"]
            S2_GEOM["Extract sun/view angles<br/>(SZA, SAA, VZA, VAA)"]
            S2_CLOUD["Extract cloud mask<br/>(cirrus + bright pixel)"]
            S2_FETCH --> S2_LOAD --> S2_GEOM
            S2_LOAD --> S2_CLOUD
        end

        subgraph G2["Group B: Atmospheric Priors"]
            ATMO_FETCH["Fetch atmo data<br/>(CAMS / MERRA-2)"]
            ATMO_PROC["Interpolate to AOI grid<br/>→ AOT, TCWV, TCO3 + σ"]
            ATMO_FETCH --> ATMO_PROC
        end

        subgraph G3["Group C: Surface BRDF Priors"]
            BRDF_FETCH["Fetch MCD43A1<br/>(earthaccess / GEE)"]
            BRDF_PROC["Temporal composite<br/>+ reproject to AOI"]
            KERN["Kernel model<br/>→ surface reflectance prior<br/>(boa_prior, boa_unc)"]
            BRDF_FETCH --> BRDF_PROC --> KERN
        end

        subgraph G4["Group D: Ancillary"]
            DEM_FETCH["Fetch DEM<br/>(Copernicus GLO-30)"]
            DEM_CROP["Crop/resample<br/>to AOI → elevation"]
            DEM_FETCH --> DEM_CROP
        end

        subgraph G5["Group E: RT Model"]
            RT_LOAD["Load emulator weights<br/>(.npz) or LUT (Zarr)"]
            RT_LOAD
        end
    end

    RESOLVE --> G1
    RESOLVE --> G2
    RESOLVE --> G3
    RESOLVE --> G4
    RESOLVE --> G5

    subgraph COLLECT["Collect & Validate"]
        BARRIER["Synchronization Barrier<br/>Wait for all groups"]
    end

    G1 --> BARRIER
    G2 --> BARRIER
    G3 --> BARRIER
    G4 --> BARRIER
    G5 --> BARRIER

    subgraph SOLVE["Aerosol / WV Retrieval (Sequential)"]
        COST["Build CostFunction<br/>J = J_obs + J_prior + J_smooth"]
        MG["Multi-grid L-BFGS-B Solve<br/>→ AOT_solved, TCWV_solved"]
        COST --> MG
    end

    BARRIER --> SOLVE

    subgraph CORRECT["Atmospheric Correction"]
        RT_APPLY["RT Model: compute<br/>xap, xbp, xcp per band"]
        BOA["Apply correction<br/>boa = (xap·toa − xbp)/(1 + xcp·y)"]
        RT_APPLY --> BOA
    end

    SOLVE --> CORRECT

    subgraph OUT["Output"]
        WRITE["Write BOA + AOT + TCWV<br/>(COG / Zarr / NetCDF)"]
    end

    CORRECT --> OUT

```
				


## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=siac --cov-report=html

# Build Rust extensions
maturin develop --release

# Type checking
mypy python/siac

# Linting
ruff check python/siac
```

## Citation

```bibtex
@article{yin2019sensor,
  title={A sensor-invariant atmospheric correction method: application to Sentinel-2/MSI and Landsat 8/OLI},
  author={Yin, Feng and Lewis, Philip E and Gomez-Dans, Jose and Wu, Qiang},
  journal={EarthArXiv},
  year={2019},
  doi={10.31223/osf.io/ps957}
}
```

## License

AGPL-3.0 - See [LICENSE](LICENSE) for details.
