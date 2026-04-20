# Native 6SV2.1 Backend Architecture

This page documents the implemented native 6SV2.1 backend used by
`algorithms.rt.backend = "sixs"`.

For user-facing setup and configuration, start with
[Native 6SV2.1 Backend](../user-guide/native-sixs-backend.md). This page
focuses on the backend architecture and runtime contract.

## Current Scope

The native backend is implemented and production-facing within the current SIAC
RT adapter surface. It provides:

- native batched 6SV2.1 execution from Python
- `mode = "direct" | "scene_lut" | "auto"`
- OpenMP-backed native execution
- `worker_libraries` fallback execution
- configurable atmospheric, aerosol, surface, and atmospheric-correction inputs
- the SIAC correction coefficients `xap`, `xbp`, `xcp`
- additional print-equivalent 6S outputs as named arrays

The backend entry points are:

- `python/siac/algorithms/rt/direct/sixs.py`
- `python/siac/algorithms/rt/direct/sixs_native.py`
- `python/siac/algorithms/rt/direct/sixs_build.py`

## High-Level Flow

```mermaid
flowchart LR
    Config["algorithms.rt.backend = sixs"] --> Runner["SixSBackend / SixSNativeRunner"]
    Runner --> Build["ensure_native_sixs_module(...)"]
    Build --> Source["fetch or reuse 6SV2.1 source"]
    Source --> Patch["patch upstream Fortran"]
    Patch --> Compile["build F2PY native module"]
    Compile --> Session["native session(s)"]
    Session --> Direct["direct batch run"]
    Session --> SceneLUT["scene-LUT run + interpolate"]
    Direct --> Outputs["xap / xbp / xcp + extras"]
    SceneLUT --> Outputs
```

## Build Model

At runtime, the backend resolves a native Python extension module. The resolved
module may come from:

- `algorithms.rt.sixs.module_path`
- the deprecated compatibility alias `algorithms.rt.sixs.library_path`
- an automatic build in the configured or default build root

The builder:

1. resolves the build root, defaulting to `~/.cache/siac/rt6s/<profile>`
2. downloads or reuses the upstream 6SV2.1 source tree
3. patches selected upstream files
4. injects a callable bridge around the 6S case core
5. generates an explicit F2PY signature that exposes only
   `sixs_f2py_run_batch`
6. compiles an F2PY extension with OpenMP-enabled Fortran flags

`build_profile = "release"` uses optimized flags. `build_profile = "parity"`
switches to lower optimization for parity work against the original executable.

## Source Patches

The native builder patches copied upstream sources before compilation.

### `DISCOM.f`

The original boundary checks can read beyond valid indices because Fortran
logical operators do not short-circuit. SIAC rewrites those checks into nested
`if` blocks.

### `KERNELPOL.f`

The upstream `is = 0` branch leaves zero-order slots uninitialized. SIAC
initializes those entries explicitly to avoid invalid values leaking into
polarized kernels.

### Additional Native-Safety Patches

The native build also carries implementation-level patches to support the Python
bridge and safe native execution, including widened file-path handling for
user-supplied aerosol models and OpenMP-related adjustments in the copied
Fortran tree. The explicit signature step is important because it prevents
F2PY from auto-generating wrapper glue for legacy COMMON blocks that the native
OpenMP build marks `THREADPRIVATE`.

## Execution Paths

### OpenMP Path

`parallel_backend = "openmp"` is the default path.

The backend:

1. loads one isolated native extension session
2. prepares contiguous scene arrays for the valid-pixel subset
3. calls one batched native entry point
4. lets the compiled extension parallelize case evaluation internally
5. writes outputs back into scene-shaped arrays, masking invalid cases with
   `NaN`

### Worker-Libraries Fallback

`parallel_backend = "worker_libraries"` is implemented as a fallback execution
model, not as the primary runtime architecture.

The backend:

1. creates multiple isolated native extension sessions
2. splits valid cases into chunks
3. runs those chunks through Python worker threads with one native session per
   worker
4. merges outputs back into scene-shaped arrays

This is the supported fallback when OpenMP is not the desired execution model.

## Route Selection

The backend implements all three route modes:

- `direct`: evaluate every valid case natively
- `scene_lut`: run a reduced RT grid and interpolate back to the scene
- `auto`: select scene-LUT only when the reduced case count is small enough to
  justify interpolation

`auto` currently checks:

- minimum valid-pixel count
- whether the LUT case count is strictly smaller than the direct case count
- whether the estimated compression ratio exceeds
  `scene_lut_required_speedup`

This is why `scene_lut` performs very well on smooth scenes and very poorly on
geometry-rich scenes: the architecture optimizes for scene compressibility, not
for every possible scene shape.

## Output Contract

The backend preserves SIAC's core correction contract:

- `xap`
- `xbp`
- `xcp`

Additional outputs are selected through `algorithms.rt.sixs.output_variables`.
Internally the native runner normalizes requests so the mandatory SIAC
coefficients are always present, even when the caller asks for a reduced output
set.

Non-core native outputs are surfaced as named extras on the runtime
coefficients object rather than changing the downstream correction contract.

## Validation Hooks

The implemented validation paths are:

- upstream parity harness:
  `python/siac/sixs_upstream_parity.py`
- integration tests:
  `tests/integration/test_sixs_backend_integration.py`

These cover:

- parity against the original upstream 6SV2.1 executable
- multi-band consistency
- repeated native-session reuse
- `scene_lut` agreement with `direct` on aligned scenes
- `worker_libraries` agreement with `openmp`

## Current Caveats

- the native backend does not expose Jacobians yet
- route quality still depends strongly on scene compressibility; `scene_lut` is
  not universally faster
- the public user workflow should treat `module_path` as the current module
  override field and `library_path` as compatibility input
