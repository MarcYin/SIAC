# Validation And Limitations

SIAC should be understood through two lenses:

- the 2022 paper as the scientific validation baseline
- the current repository as the living implementation and architecture

## What the paper supports

The paper presents SIAC as a Bayesian atmospheric correction method for Sentinel-2/MSI and Landsat 8/OLI and validates it using several lines of evidence, including:

- comparisons against AERONET atmospheric measurements
- comparisons against RadCalNet surface reflectance
- inter-sensor consistency checks between Sentinel-2 and Landsat 8
- analysis of uncertainty behavior

Those claims belong to the publication and should be cited as such.

## What the current repository clearly exposes

From the current codebase, the following are directly visible:

- a modular config-driven implementation
- a staged M1-M6 runtime pipeline
- pluggable atmospheric and BRDF provider choices
- explicit output control for BOA, uncertainty, and auxiliary products
- a Rust acceleration boundary for selected numerical primitives
- a primary user-facing workflow for Sentinel-2

## Limits and caution points

| Area | Caution |
| --- | --- |
| Scientific claims | Use the paper for validated scientific claims, not only the code structure. |
| Sensor coverage | Sentinel-2 is the clearest user-facing path in the current repo; Landsat is present in the API but should be documented carefully. |
| Backend maturity | The schema supports multiple providers and backends, but not every option should automatically be treated as equally mature operationally. |
| Uncertainty interpretation | The repository carries uncertainty fields, but the paper remains the main source for the derivation and interpretation of uncertainty behavior. |
| Operational assumptions | Credentials, caches, remote access, and optional dependencies introduce real environment-dependent behavior. |

## Practical guidance

- Cite the paper when discussing performance, validation, and method quality.
- Cite the repository when discussing architecture, package layout, and current interfaces.
- Keep tutorial claims conservative unless they are directly supported by tested paths in the current repository.

## Recommended wording in downstream docs

- Say "the SIAC method is described and validated in the 2022 paper".
- Say "the current repository implements a modular SIAC v2 architecture".
- Avoid implying that every configurable backend has the same validation status as the paper's primary published workflow.
