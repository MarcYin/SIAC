# Priors And Mapping

SIAC depends on more than raw scene observations. It uses atmospheric priors, surface priors, and scale/spectral mapping steps to make the retrieval tractable and physically constrained.

## Atmospheric priors

The repository supports multiple atmospheric-prior kinds through configuration, including:

- `cams`
- `merra2`
- `mcd19`
- `vnp19`
- `era5`
- `user`

These are represented in `AtmoProviderConfig` and resolved through the atmospheric provider assembly path. In the workflow, the atmospheric prior produces an `AtmosphericState` containing:

- `aot`
- `tcwv`
- `tco3`
- uncertainty fields
- elevation

## Surface priors

Surface priors estimate plausible BOA reflectance before the final correction stage. They help constrain the retrieval by providing a surface-side expectation rather than forcing the atmospheric state to explain everything.

Configured surface-prior methods currently include:

- `kernel_model`
- `whittaker`
- `monthly_database`
- `neural`
- `direct`

Main implementation areas:

- `python/siac/algorithms/surface/kernel_model.py`
- `python/siac/algorithms/surface/brdf_whittaker.py`
- `python/siac/algorithms/surface/brdf_monthly_database.py`
- `python/siac/algorithms/surface/swir_refine.py`

## BRDF prior products

The paper emphasizes MODIS BRDF products as a key surface prior. The current repository generalizes the provider layer and supports multiple BRDF provider kinds, including `mcd43` and `vnp43`, with provider resolution handled through config plus assembly.

## Spectral mapping

Spectral mapping matters when the source prior product and the target sensor do not share the same effective spectral response. The repository treats this as an explicit algorithmic concern rather than a hidden conversion.

Main implementation area:

- `python/siac/algorithms/surface/spectral_mapping.py`

Relevant config area:

- `algorithms.surface_prior.spectral_mapping`

## Spatial mapping and ePSF

The paper's spatial mapping idea appears in the repository as PSF-aware or scale-aware handling when matching coarser prior information to finer scene observations. In practical terms, this is why SIAC separates solver resolution from scene resolution and why PSF and related utilities appear both in Python and Rust-backed code.

Relevant areas:

- `python/siac/algorithms/surface/swir_refine.py`
- `src/siac_rs/src/psf.rs`

For the Route-B monthly-database query specifically, the initial NIR/SWIR
correction is intentionally done on the coarse query grid: SIAC area-resamples
the native TOA query bands to the target solving/database resolution first,
then applies the first-pass correction on that coarse grid. This is separate
from M6, where the final atmospheric correction is still applied back at the
native band resolutions.

## Cloud masking

Cloud handling is part of the upstream observation preparation and affects which pixels contribute to later stages. The current repository exposes cloud-mask configuration through:

- `algorithms.cloud_mask.mode`
- `algorithms.cloud_mask.provider`
- `algorithms.cloud_mask.target_resolution_m`

Main implementation area:

- `python/siac/algorithms/cloud/mask.py`

## Grid resampling between stages

Moving data between spatial grids (native sensor resolution → solver grid →
correction grid) is a recurring concern across the pipeline. All resampling is
handled by the canonical functions in `python/siac/geo/resample.py`:

- `resample_field_to_template` — bilinear interpolation for continuous fields
  (atmospheric state, geometry angles), with coordinate-based interpolation
  when possible and `scipy.ndimage.zoom` as fallback. NaN gap-fill is applied
  by default.
- `resample_mask_to_template` — conservative resampling for boolean masks
  (cloud, shadow), applying a maximum-filter dilation before interpolation to
  avoid losing masked pixels when downsampling.
- `resample_coefficients_to_template` — resamples all fields of an
  `RTCoefficients` bundle, handling both 2-D and 3-D (with `param` dimension)
  arrays.
- `resample_field_for_correction` — combines resampling with a finiteness
  guarantee (nearest-neighbour fill then source-mean fallback) for the
  correction stage.

These functions are used by the grid assembler (M4), the solver (M5), and the
atmospheric corrector (M6).

## M1-M6 stage mapping

| Stage | Prior or mapping role |
| --- | --- |
| M1 | prepares TOA, geometry, and cloud mask |
| M2 | materializes atmospheric priors |
| M3 | materializes surface priors, possibly using BRDF routes and mapping |
| M4 | moves all pieces onto the solver grid |
| M5 | solves the atmospheric state using those constraints |
| M6 | applies the solved state back into correction space |

## What changed from the paper-era implementation

The paper should be treated as the scientific baseline, but the repository is more modular than a paper-oriented narrative:

- provider choice is surfaced through config rather than embedded in one path
- RT backend selection is explicit
- spectral mapping is exposed as its own configurable component
- the execution flow is split into typed runtime payloads and assembly/orchestration layers
