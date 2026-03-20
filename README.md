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

Remote-auth can be provided in `config.toml` or via the centralized env overlay:

- Earthdata / earthaccess: `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`
- CDS API: `CDSAPI_KEY`
- AWS/S3: `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`
- GCS: `GOOGLE_APPLICATION_CREDENTIALS`
- CDSE: `SIAC_CDSE_USERNAME`, `SIAC_CDSE_PASSWORD`

For CDSE object-store datasets, SIAC can mint temporary S3 credentials from the
same CDSE username/password and use them against `s3://eodata/...`.

## Quick Start

```python
from siac import SIAC
from siac.config import SIACConfig

# Load configuration
config = SIACConfig.from_file("siac.toml")

# Or use defaults
config = SIACConfig()

# Process Sentinel-2 data
siac = SIAC(config)
result = siac.process("/path/to/S2_SAFE/")

# Access results
print(result.boa)  # BOA reflectance dataset
print(result.aot.mean())  # Mean AOT
```

## Full Sentinel-2 Run (Real Data)

Use the helper script to run the full Sentinel-2 workflow (query/download + M1-M6):

```bash
PYTHONPATH=python pixi run python tools/run_full_s2.py \
  --query "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433" \
  --backend gcs \
  --cache-dir ~/.cache/siac/s2 \
  --aoi-bbox 119.35 -25.10 119.42 -25.03 \
  --aoi-crs EPSG:4326 \
  --output-dir ./outputs/s2_full_run
```

Notes:
- `--query` accepts a local SAFE path, product ID, or tile-date shorthand (`T31UDQ_20210801`).
- Use `--aoi-bbox` to keep runs tractable on local machines.
- Outputs are written under `output-dir/boa/`, `output-dir/auxiliary/`, plus summary/STAC files from the helper script.

## Configuration

SIAC uses a hierarchical TOML configuration system. Create a `siac_config.toml`:

```toml
[paths]
lut_path = "https://gws-access.jasmin.ac.uk/public/nceo_isp/libradtran_continental_average_lut_1nm.zarr.zip"
spectral_library_root = "/data/siac/spectral-library"
rsrf_root = "/data/siac/RSRF"

[providers.atmo]
kind = "cams"

[providers.brdf]
kind = "mcd43"
temporal_window = 16

[providers.s2]
backend = "cdse"
max_cloud_cover = 40.0

[algorithms.surface_prior]
method = "kernel_model"

[algorithms.rt]
backend = "emulator"

[algorithms.solver]
aot_gamma = 10.0
tcwv_gamma = 5.0
aerosol_resolution = 1000.0

[output.defaults]
format = "cog"
include_uncertainty = true
```

## Architecture

```
siac/
├── api/            # Public entrypoints
├── config/         # System config, TOML loading, resolution
├── catalog/        # Built-in sensor catalog data
├── domain/         # Pure domain types and spectral helpers
├── adapters/       # External systems and backend adapters
├── algorithms/     # Retrieval, correction, RT, cloud, surface logic
├── app/            # Runtime planning and component assembly
├── workflows/      # End-to-end orchestration
├── runtime/        # Xarray-backed execution payloads
├── geo/            # Geometry/reprojection utilities
└── storage/        # Raster and product I/O helpers
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

# Run with strict coverage gates (global >=95%, each file >90%)
pytest tests/ --cov=siac --cov-report=json:coverage.json --cov-report=xml:coverage.xml --cov-report=html
python tools/check_coverage_thresholds.py --min-total 95 --min-file 90

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
