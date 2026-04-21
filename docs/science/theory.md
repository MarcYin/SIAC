# Theory

SIAC frames atmospheric correction as a constrained estimation problem. The goal is to recover a physically meaningful atmospheric state and then use it to correct top-of-atmosphere observations into bottom-of-atmosphere reflectance.

The scientific baseline for this repository is the 2022 SIAC paper:

- Yin, F., Lewis, P. E., and Gomez-Dans, J. L. (2022), [Bayesian atmospheric correction over land: Sentinel-2/MSI and Landsat 8/OLI](https://gmd.copernicus.org/articles/15/7933/2022/)

## High-level idea

SIAC combines:

- observed TOA reflectance from the scene
- solar and viewing geometry
- atmospheric priors such as AOT and TCWV
- surface priors derived from BRDF products
- a radiative transfer representation through an emulator, LUT, or native 6S backend
- spatial smoothness constraints in the atmospheric retrieval stage

Instead of treating atmospheric correction as a purely forward transform, SIAC treats the atmospheric state as something to be estimated under observation and prior constraints.

## Retrieval then correction

The codebase follows the same broad pattern as the paper:

1. build a normalized observation bundle from the scene
2. estimate or fetch atmospheric priors
3. estimate or fetch a surface prior
4. assemble solver inputs at aerosol retrieval resolution
5. solve for atmospheric state
6. apply correction at the output resolution

```mermaid
flowchart LR
    TOA["TOA reflectance + geometry"] --> PriorA["Atmospheric prior"]
    TOA --> PriorS["Surface prior"]
    PriorA --> Solve["Bayesian-style retrieval"]
    PriorS --> Solve
    Solve --> State["Solved atmosphere"]
    State --> Correct["Atmospheric correction"]
    TOA --> Correct
    Correct --> BOA["BOA reflectance + uncertainty"]
```

## Observation, prior, and smoothness terms

At a conceptual level, the retrieval balances three things:

| Term | Role |
| --- | --- |
| Observation term | Penalizes mismatch between observations and the modeled BOA-to-TOA relationship |
| Prior term | Penalizes departures from prior atmospheric and surface expectations |
| Smoothness term | Encourages spatially coherent atmospheric retrievals at the solver scale |

In the current repository this logic appears across:

- `python/siac/algorithms/grid/assembler.py`
- `python/siac/algorithms/solver/cost.py`
- `python/siac/algorithms/solver/multigrid.py`
- `python/siac/algorithms/correction/atmospheric.py`

### Cost function detail

The cost function in `algorithms/solver/cost.py` combines the three terms into a single scalar minimized by L-BFGS-B:

$$J = J_\text{obs} + J_\text{prior} + J_\text{smooth}$$

**Observation term** — For each solver band, the forward model predicts TOA reflectance from the current atmospheric state and surface prior via `RTModelBackend.simulate_toa`. The residual is weighted by band-dependent weights proportional to $\lambda^{\alpha}$ where $\alpha$ is the `band_weight_power` parameter (default −1.6, configured as `algorithms.solver.alpha`). This gives shorter wavelengths — more sensitive to aerosol — higher influence.

**Prior term** — Penalizes departure of AOT and TCWV from their prior values, weighted by prior uncertainty. Regularization strengths `aot_gamma` and `tcwv_gamma` scale the prior cost.

**Smoothness term** — A DCT-based (Discrete Cosine Transform) spatial regularization. Rather than constructing and inverting a full Laplacian matrix, the implementation builds a sparse Laplacian with boundary-aware diagonal, multiplies in DCT space, and penalizes high-frequency structure in the AOT and TCWV fields. This encourages smooth atmospheric fields without requiring expensive matrix operations.

### Multi-grid solver strategy

The default solver (`algorithms/solver/multigrid.py`) uses a coarse-to-fine multi-grid approach:

1. Start at a coarse grid (e.g. 4× the target aerosol resolution).
2. Solve for AOT and TCWV using L-BFGS-B at this resolution.
3. Interpolate the solution to the next finer grid level.
4. Repeat until the target `aerosol_resolution` is reached.

Each level uses the previous level's solution as the initial guess, improving convergence speed and avoiding local minima. The `MultiGridConfig` controls bounds, regularization strengths, and the `band_weight_power`.

## Why coarse retrieval comes before full-resolution correction

The paper emphasizes retrieving atmospheric parameters at a coarser resolution than the native scene bands. The same architectural idea is visible in the code:

- observations are read at scene resolution
- atmospheric and surface priors are assembled onto the solver grid
- the solver retrieves atmospheric state at the configured aerosol resolution
- the solved state is then used to correct BOA at the output grid

This is important because the atmospheric state varies more smoothly than raw TOA observations, and the prior products used for SIAC are often coarser than the scene itself.

## Uncertainty in SIAC

SIAC is designed to carry uncertainty explicitly rather than treat the correction as a black box. The paper treats uncertainty as part of the scientific case for SIAC, and the repository reflects that by carrying uncertainty-bearing fields through runtime payloads such as:

- `AtmosphericState`
- `SurfacePrior`
- `SolvedAtmosphere`
- `CorrectionResult`

The current docs should be read as describing the uncertainty structure present in the codebase, while the paper remains the authoritative scientific explanation for the derivation and validation logic.

## Code-level mapping

| Theory concept | Main code area |
| --- | --- |
| scene normalization | `adapters.satellite`, `runtime.models.ObservationBundle` |
| atmospheric priors | `adapters.atmo` |
| surface priors | `algorithms.surface`, `adapters.brdf` |
| solver grid construction | `algorithms.grid` |
| state retrieval | `algorithms.solver` |
| atmospheric correction | `algorithms.correction`, `adapters.rt`, `algorithms.rt` |

## Reading path

- Continue to [Priors and Mapping](priors-and-mapping.md) for the main scientific building blocks.
- Continue to [Paper Companion](paper-companion.md) for a paper-to-code crosswalk.
