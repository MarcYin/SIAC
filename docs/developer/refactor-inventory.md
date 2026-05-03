# Refactor Inventory

This inventory records behavior that must stay stable while SIAC internals are
modernized. Treat it as the checklist for small, reviewable refactor passes.

## Guardrails

- Keep the documented public API stable: `siac.SIAC`, `siac.SIACConfig`,
  `process_sentinel2`, `resolve_s2_input`, `search_sentinel2`,
  `siac_process_s2`, and the `siac process-s2` CLI.
- Do not reintroduce removed compatibility import paths:
  `siac.config.schema`, `siac.app.assembly`, `siac.app._assembly_runtime`,
  `siac.algorithms.rt.lut.srf_kernel`, `siac.adapters.brdf.vnp43_earthaccess`,
  or `siac.process_landsat8`.
- Move implementation into owning modules and update callers to canonical paths
  in the same pass.
- Every structural pass should have a targeted unit test plus the fast suite as
  the behavior check.

## Proposed Passes

| Pass | Current behavior | Structural improvement | Validation check |
| --- | --- | --- | --- |
| Removed-path guardrails | Compatibility modules and aliases are absent. | Keep negative import/export tests so deleted paths stay deleted. | `pytest tests/unit/test_removed_compatibility_paths.py` |
| Pipeline helper extraction | `workflows.pipeline` owns orchestration plus diagnostics, output-stream handling, and result shaping. | Move diagnostics and result-building helpers into focused workflow modules without changing `run_pipeline`. | Pipeline, orchestration, STAC, and output adapter unit tests. |
| Assembly extraction | Focused `_assembly_*` modules own preprocessor, provider, surface, solver, correction, RT, I/O, and S2 backend assembly. | Keep callers on the focused modules and avoid umbrella assembly re-exports. | Assembly helper tests and planning tests. |
| BRDF provider split | `adapters.brdf.mcd43_earthaccess` contains product specs, stack loading, and several providers. | Extract product metadata and stack-loading helpers before moving concrete providers. | MCD43/VNP43 Earthaccess path tests. |
| Solver split | `MultiGridSolver` and grid-search helpers are oversized. | Extract input coercion, QA/result building, and grid-search refinement helpers. | Solver unit tests, Rust compatibility tests, and synthetic pipeline tests. |
| Native RT split | 6S build/patch/native execution code is concentrated in large modules. | Separate source patching, build discovery, runner setup, and execution adapters. | `tests/unit/test_sixs_native.py`, RT backend tests, optional native smoke. |
| Surface-prior split | `swir_refine` exposes monthly database, Route-B query, and workflow helpers. | Move monthly database/query helpers behind explicit modules with direct canonical imports. | SWIR refine, monthly composite, and spectral mapping tests. |

## Migration-Only Work

These are not behavior-preserving refactors and should be split into separate
tasks:

- Adding an end-to-end Landsat workflow.
- Changing Dask/thread execution semantics.
- Upgrading dependency bounds for runtime behavior.
- Reshaping the Rust extension API or removing `_rust_compat`.
