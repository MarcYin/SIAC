# Paper Companion

This page helps scientific users move between the 2022 SIAC paper and the current SIAC v2 repository layout.

Reference paper:

- [Bayesian atmospheric correction over land: Sentinel-2/MSI and Landsat 8/OLI](https://gmd.copernicus.org/articles/15/7933/2022/)

## How to use this page

- Read the paper for the scientific rationale, derivation, validation, and terminology.
- Use this page to find where those ideas live in the current codebase.
- Use the architecture docs to understand how the modular v2 implementation packages those ideas into layers.

## Paper concept to code map

| Paper topic | Current repository area |
| --- | --- |
| atmospheric correction scheme in SIAC | `python/siac/workflows/pipeline.py`, `python/siac/algorithms/correction/atmospheric.py` |
| Bayesian retrieval framing | `python/siac/algorithms/solver/cost.py`, `python/siac/algorithms/solver/multigrid.py` |
| atmospheric priors (AOT, TCWV, ozone) | `python/siac/adapters/atmo/` |
| BRDF-based surface prior | `python/siac/adapters/brdf/`, `python/siac/algorithms/surface/` |
| spectral mapping | `python/siac/algorithms/surface/spectral_mapping.py` |
| spatial mapping / ePSF | `python/siac/algorithms/surface/swir_refine.py`, `src/siac_rs/src/psf.rs` |
| radiative transfer approximation | `python/siac/adapters/rt.py`, `python/siac/algorithms/rt/`, `src/siac_rs/src/emulator.rs` |
| grid resampling between stages | `python/siac/geo/resample.py` |
| implementation details | `python/siac/app/`, `python/siac/workflows/`, `python/siac/runtime/` |
| uncertainty handling | `python/siac/runtime/models.py`, solver/correction paths |

## Section-style crosswalk

| Paper area | Read in the repo |
| --- | --- |
| Introduction and motivation | [Overview](../concepts/overview.md), [Theory](theory.md) |
| atmospheric correction scheme | [Execution Flow](../architecture/execution-flow.md), `algorithms/solver`, `algorithms/correction` |
| materials and methods | `adapters/`, `algorithms/`, `runtime/`, [Priors And Mapping](priors-and-mapping.md) |
| appendices on spectral and spatial mapping | `algorithms/surface/spectral_mapping.py`, PSF-aware surface prior logic |
| implementation details appendix | `app/`, `workflows/`, `runtime/`, [Package Map](../architecture/package-map.md) |
| validation and conclusions | [Validation And Limitations](validation-and-limitations.md) |

## Important differences to keep in mind

The repository and the paper are related but not identical artifacts.

| Paper-first view | Repo-first view |
| --- | --- |
| organized around method explanation | organized around modular package boundaries |
| emphasizes scientific derivation and validation | emphasizes configuration, runtime flow, and extension points |
| focused on the published SIAC method | exposes a broader pluggable framework for providers and backends |

## Reading sequence for scientific users

1. Read the paper sections on the SIAC scheme and appendices of interest.
2. Read [Theory](theory.md) and [Priors And Mapping](priors-and-mapping.md).
3. Read [Execution Flow](../architecture/execution-flow.md) and [Data Flow](../architecture/data-flow.md).
4. Use the table above to jump into the relevant modules.
