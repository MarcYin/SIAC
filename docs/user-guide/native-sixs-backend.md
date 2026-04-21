# Native 6SV2.1 Backend

SIAC includes a native 6SV2.1 backend exposed as `algorithms.rt.backend = "sixs"`.
This backend is separate from the LUT and emulator backends:

- `sixs`: native compiled 6SV2.1 with batched array execution, OpenMP support,
  scene-LUT routing, and selectable 6S report outputs.
- `lut` and `emulator`: alternative SIAC RT backends that do not require a
  native 6S build.

Use the native backend when you want direct access to 6SV2.1 behavior, parity
validation against upstream 6S, or 6S print-equivalent diagnostics beyond the
core correction coefficients.

For a measured comparison between native 6S direct, native 6S `scene_lut`, and
the remote libRadtran ZIP/Zarr LUT, see
[RT Model Differences](rt-model-differences.md).
That remote LUT route is now modeled explicitly as a fixed libRadtran preset,
not as a configurable semantic peer of native 6S.

## What The Backend Returns

The native backend always returns the three SIAC correction coefficients:

- `xap`
- `xbp`
- `xcp`

Those arrays are always present internally even if you request a smaller
`output_variables` list. Additional named outputs are exposed through the same
result object and can be requested selectively.

By default, `algorithms.rt.sixs.output_variables` includes the full currently
supported native output surface.

SIAC uses those three coefficients in the BOA correction equation:

```text
y = xap * toa - xbp
boa = y / (1 + xcp * y)
```

That is why cross-backend comparisons in the user guide focus on `xap`, `xbp`,
and `xcp`: they are the direct interface between the RT backend and the final
surface-reflectance calculation.

## Build And Setup

### Prerequisites

The native backend requires:

- a Fortran compiler such as `gfortran`
- `meson` and `ninja` for the F2PY Meson build path when that backend is available
- the normal SIAC Python package environment

The Pixi workspace includes an `rt6s` feature/environment for this toolchain.
That environment pins Python 3.11 and `setuptools < 60`, which keeps the
legacy F2PY `distutils` backend available for the native 6S extension. Outside
that environment, SIAC prefers the Meson backend and only offers `distutils`
when NumPy and `setuptools` still support it.

### Pixi Path

Use the `rt6s` environment when you want the repo-managed toolchain:

```bash
pixi install
pixi run -e rt6s build-6s-native
```

That task runs `tools/build_6s_native.py`, which builds the native Python
extension. In the repo-managed environment, the hosted Linux smoke workflow now
uses `SIAC_SIXS_F2PY_BACKEND=meson,distutils`, so GitHub Actions validates the
Meson path first and only falls back to the legacy `distutils` backend if Meson
fails unexpectedly. Outside that environment, the builder still tries the
Meson backend first and only exposes `distutils` when the current NumPy and
`setuptools` combination supports it.

When that dedicated Linux smoke workflow fails, SIAC now writes the native-build
diagnostics under `tmp/rt6s_ci_smoke/diagnostics/`, prints the key log tails in
the job output, emits check annotations for the key Meson or distutils stderr
tails plus the aggregated build failure summary, and uploads the full
`tmp/rt6s_ci_smoke` directory as a GitHub Actions artifact. Use that artifact
first when a failure reproduces only on the hosted runner, and use the
annotations when you need a quick public summary.

To remove the native build cache:

```bash
pixi run -e rt6s clean-6s-native
```

### Manual Path

If you are not using Pixi, install SIAC in an environment that already has the
native build prerequisites, then run:

```bash
python tools/build_6s_native.py
```

Useful overrides:

```bash
python tools/build_6s_native.py \
  --source-dir /path/to/6sV2.1 \
  --build-dir /tmp/siac-rt6s \
  --compiler gfortran \
  --build-profile release
```

### Build Profiles

- `release`: optimized production build. This is the default.
- `parity`: lower-optimization build for one-to-one parity work against the
  upstream executable.

The default build cache root is:

```text
~/.cache/siac/rt6s/<build_profile>
```

### Source Handling

The native builder can either:

- download the upstream 6SV2.1 tarball from `algorithms.rt.sixs.source_url`, or
- use an existing unpacked source tree from `algorithms.rt.sixs.source_dir`

During build, SIAC patches the copied 6S sources before compiling the Python
extension. It also generates an explicit F2PY signature for
`sixs_f2py_run_batch`, so the compiled Python module only exposes the batched
array entrypoint rather than auto-wrapping legacy 6S COMMON blocks. The build
path is part of the supported user workflow; you do not need to patch 6S
manually.

## Minimal Configuration

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.sixs]
mode = "auto"
parallel_backend = "openmp"
build_profile = "release"
native_threads = 8
```

If you want SIAC to use a prebuilt module instead of auto-building it, set
`module_path` and disable auto-build:

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.sixs]
auto_build = false
module_path = "/opt/siac/_siac_rt6s_native.so"
```

## Shared RT Setup Layer

The native 6S backend now sits under a generic RT setup layer:

- `algorithms.rt.setup.*` describes the requested RT semantics
- `algorithms.rt.sixs.*` carries native execution/build controls

This is now the only supported semantic split. Atmosphere, aerosol, surface,
and correction setup belong under `algorithms.rt.setup.*`; runtime choices
such as build path, route mode, and thread count belong under
`algorithms.rt.sixs.*`.

### Example: Generic RT Semantics Plus Native 6S Runtime Controls

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.setup.atmosphere]
profile = "us_standard_62"
columns_mode = "input_columns"

[algorithms.rt.setup.aerosol]
profile = "continental"

[algorithms.rt.setup.surface]
mode = "homogeneous_lambertian"

[algorithms.rt.setup.surface.target]
kind = "constant"
constant = 0.1

[algorithms.rt.setup.atmospheric_correction]
mode = "lambertian_reflectance"
value = 0.1

[algorithms.rt.sixs]
mode = "auto"
parallel_backend = "openmp"
native_threads = 8
build_profile = "release"
```

This split is deliberate:

- the RT setup can be shared conceptually across backends
- native 6S still needs its own build/runtime controls
- advanced 6S modes such as radiosonde atmospheres, layered aerosols, and
  custom Mie-model files are expressed through the shared `rt.setup.*` layer

## Starting From User Atmospheric State

The full SIAC pipeline is retrieval-first: it derives an atmospheric state and
then applies correction. The native 6S backend does **not** require you to use
that full path.

If you already have:

- TOA reflectance
- viewing and illumination geometry
- atmospheric state fields such as AOT, TCWV, TCO3, and elevation
- a native 6S surface and atmospheric-correction configuration

you can start directly from the RT/correction stage and bypass M1 to M5.

That is the right entrypoint when:

- you trust an external atmospheric-state source
- you want to compare SIAC's correction contract against another retrieval
  system
- you want to run native 6S with custom BRDF or heterogeneous-Lambertian
  surface settings without running the full SIAC retrieval workflow

### Coefficients-Only Entry Point

Use this path when you want `xap`, `xbp`, `xcp`, or other selected 6S outputs
but you do not want SIAC to apply BOA correction for you yet.

For direct programmatic use, construct `SixSBackend` with both:

- a `SixSAlgorithmConfig` for native runtime controls
- an `RTSetupConfig` for the RT semantics

```python
from datetime import datetime

import numpy as np
import xarray as xr

from siac.algorithms.rt.direct.sixs import SixSBackend
from siac.config import RTSetupConfig, SixSAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.runtime import AtmosphericState, GeometryAngles

shape = (256, 256)

geometry = GeometryAngles.from_degrees(
    xr.DataArray(np.full(shape, 30.0, dtype=np.float32), dims=("y", "x")),
    xr.DataArray(np.full(shape, 150.0, dtype=np.float32), dims=("y", "x")),
    xr.DataArray(np.full(shape, 5.0, dtype=np.float32), dims=("y", "x")),
    xr.DataArray(np.full(shape, 110.0, dtype=np.float32), dims=("y", "x")),
)

atmo = AtmosphericState(
    aot=xr.DataArray(np.full(shape, 0.15, dtype=np.float32), dims=("y", "x")),
    tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=("y", "x")),
    tco3=xr.DataArray(np.full(shape, 0.30, dtype=np.float32), dims=("y", "x")),
    aot_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
    tcwv_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=("y", "x")),
    tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=("y", "x")),
    elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=("y", "x")),
)

band = SensorBand(
    name="B04",
    center_wavelength=665.0,
    bandwidth=30.0,
    resolution=10.0,
    band_index=3,
)

backend = SixSBackend(
    sixs_config=SixSAlgorithmConfig(
        source_dir="tmp/6s_upstream",
        build_profile="release",
        mode="direct",
        output_variables=("xap", "xbp", "xcp", "tgasm", "sutott", "sast"),
    ),
    rt_setup=RTSetupConfig(
        atmosphere={"profile": "us_standard_62", "columns_mode": "input_columns"},
        aerosol={"profile": "continental"},
        surface={
            "mode": "homogeneous_lambertian",
            "target": {"kind": "constant", "constant": 0.1},
        },
        atmospheric_correction={"mode": "lambertian_reflectance", "value": 0.1},
        reference_reflectance=0.1,
    ),
)
backend.set_observation_time(datetime(2025, 7, 12, 10, 30))

coeffs = backend.compute_coefficients(geometry, atmo, band)

xap = coeffs.xap
xbp = coeffs.xbp
xcp = coeffs.xcp
tgasm = coeffs.get_output("tgasm")
```

This path uses your provided `AtmosphericState` directly. No atmospheric
retrieval is run.

### Direct TOA-To-BOA Entry Point

If you already have TOA reflectance on the same grid, you can apply correction
immediately with the returned coefficients:

```python
toa_b04 = xr.DataArray(
    np.full(shape, 0.12, dtype=np.float32),
    dims=("y", "x"),
)

boa_b04 = coeffs.apply_correction(toa_b04)
```

For multiple bands, compute a coefficient object per band and apply the same
pattern:

```python
coeff_list = backend.compute_coefficients_multi(geometry, atmo, [band1, band2, band3])
boa_band1 = coeff_list[0].apply_correction(toa_band1)
```

This is the simplest way to start “early” from user-provided atmospheric state:
you use native 6S as the forward RT engine and SIAC's `RTCoefficients` object
as the correction contract.

### Using The M6 Corrector With Your Own Atmospheric State

If you already have a multi-band TOA dataset and a matching `SensorConfig`, you
can also call the correction stage directly:

```python
from siac.algorithms.correction.atmospheric import AtmosphericCorrector

corrector = AtmosphericCorrector(
    rt_model=backend,
    sensor_config=sensor_config,
    correction_workers=4,
)

result = corrector.correct(
    toa=toa_dataset,
    geometry=geometry,
    atmo_state=atmo,
    cloud_mask=cloud_mask,
)

boa_dataset = result.boa
```

This still bypasses the full retrieval workflow. The difference is only that
the M6 corrector handles the band loop, coefficient resampling, BOA masking,
and `CorrectionResult` assembly for you.

### What You Must Provide

When you bypass M1 to M5, SIAC no longer estimates the atmosphere for you. You
must provide:

- geometry arrays on the grid you want to correct
- `AtmosphericState` arrays with finite `aot`, `tcwv`, `tco3`, and `elevation`
- TOA reflectance on a compatible grid if you want BOA output
- any native 6S surface/BRDF configuration required by your chosen correction
  mode

The native 6S backend will still validate combinations such as:

- BRDF correction modes requiring `surface.mode = "homogeneous_brdf"`
- heterogeneous Lambertian surfaces requiring `setup.surface.environment`
- aerosol/user-profile modes requiring their companion config blocks

## Surface Configuration For Early-Start Use

The native 6S surface configuration is available whether you use the full SIAC
pipeline or start directly from a user-provided atmospheric state.

### Heterogeneous Lambertian Example

```toml
[algorithms.rt.setup.surface]
mode = "heterogeneous_lambertian"
radius_km = 0.5

[algorithms.rt.setup.surface.target]
kind = "constant"
constant = 0.12

[algorithms.rt.setup.surface.environment]
kind = "built_in"
built_in = "green_vegetation"

[algorithms.rt.setup.atmospheric_correction]
mode = "lambertian_reflectance"
value = 0.12
```

### BRDF Example

```toml
[algorithms.rt.setup.surface]
mode = "homogeneous_brdf"

[algorithms.rt.setup.surface.target]
kind = "constant"
constant = 0.15

[algorithms.rt.setup.surface.brdf]
model = "rahman"

[algorithms.rt.setup.surface.brdf.parameters]
intensity = 0.12
asymmetry_factor = 0.03
structural_parameter = 0.45

[algorithms.rt.setup.atmospheric_correction]
mode = "brdf_reflectance"
value = 0.1
```

In both cases, the same configuration can be used for:

- end-to-end SIAC runs, or
- manual direct-correction flows where you supply `AtmosphericState` yourself

## Configuration Surface

### Execution And Build Controls

Key `algorithms.rt.sixs` fields for native execution:

| Field | Purpose |
| --- | --- |
| `source_url` | Upstream 6S tarball URL used when `source_dir` is not supplied |
| `source_dir` | Existing unpacked 6SV2.1 source tree |
| `build_dir` | Override the native build root |
| `module_path` | Use a specific compiled Python extension path |
| `auto_build` | Build automatically on first use when the module is missing |
| `compiler` | Fortran compiler executable name or path |
| `build_profile` | `release` or `parity` |
| `mode` | `direct`, `scene_lut`, or `auto` |
| `parallel_backend` | `openmp` or `worker_libraries` |
| `native_threads` | Native OpenMP worker count, or the worker count used to size fallback sessions |
| `worker_libraries` | Number of isolated native library copies for `parallel_backend = "worker_libraries"` |
| `chunk_size` | Cases per chunk in `worker_libraries` mode |

### Route Selection Controls

| Field | Purpose |
| --- | --- |
| `scene_lut_min_pixels` | Minimum valid-pixel count before `auto` will consider scene-LUT routing |
| `scene_lut_max_nodes_per_axis` | Initial grid resolution per interpolated axis |
| `scene_lut_max_cases` | Hard cap on scene-LUT cases after axis reduction |
| `scene_lut_required_speedup` | Minimum `direct_case_count / lut_case_count` ratio for `auto` to switch to scene-LUT |

### Generic RT Semantics

The following semantic inputs are configured under `algorithms.rt.setup.*`.

Supported `setup.atmosphere.profile` values:

- `no_gas`
- `tropical`
- `midlatitude_summer`
- `midlatitude_winter`
- `subarctic_summer`
- `subarctic_winter`
- `us_standard_62`
- `auto_latitude_date`
- `user_water_ozone`
- `user_profile`

Operational notes:

- `auto_latitude_date` requires `setup.atmosphere.profile_latitude`. The backend
  selects a built-in 6S atmosphere from latitude and month.
- `user_profile` requires `setup.atmosphere.radiosonde_profile`.
- `user_water_ozone` uses SIAC-provided atmospheric state values rather than a
  built-in profile deck.
- `setup.atmosphere.columns_mode = "input_columns"` keeps the selected atmospheric
  profile shape but scales total water vapour and ozone to the scene `tcwv` and
  `tco3` inputs. This matches the original 6S capability where preset
  atmospheric profiles can still be combined with user-specified column totals.
  It is the default behavior for the native 6S backend.
- `setup.atmosphere.columns_mode = "profile_default"` keeps the selected profile's
  built-in water-vapour and ozone totals.

In practice, this means you can use:

```toml
[algorithms.rt.setup.atmosphere]
profile = "us_standard_62"
columns_mode = "input_columns"
```

and still have SIAC drive the native 6S run with scene-specific `tcwv` and
`tco3`, while preserving the US62 vertical profile shape.

### Aerosol Profiles

Supported `setup.aerosol.profile` values:

- `none`
- `continental`
- `maritime`
- `urban`
- `desert`
- `biomass_burning`
- `stratospheric`
- `user_mixture`
- `multimodal_log_normal`
- `modified_gamma`
- `junge_power_law`
- `sun_photometer`
- `layered_profile`
- `user_model`

Required companion fields:

- `user_mixture` requires `setup.aerosol.mixture`, which must sum to `1.0`
  in dust, water, oceanic, soot order.
- `multimodal_log_normal`, `modified_gamma`, and `junge_power_law` require
  `setup.aerosol.distribution`.
- `sun_photometer` requires `setup.aerosol.sun_photometer_aerosol`.
- `layered_profile` requires `setup.aerosol.layer_profile`. All layers must use the
  same `aerosol_type`.
- `user_model` requires `setup.aerosol.model_path`.

### Surface And Atmospheric-Correction Modes

Supported surface modes:

- `homogeneous_lambertian`
- `heterogeneous_lambertian`
- `homogeneous_brdf`

Surface reflectance targets use:

- `kind = "constant"` with `constant`
- `kind = "built_in"` with one of `green_vegetation`, `clear_water`, `sand`,
  or `lake_water`
- `kind = "spectrum"` with `values` and optional `wavelengths_um`

Supported BRDF models under `setup.surface.brdf.model`:

- `user_defined`
- `hapke`
- `verstraete`
- `roujean`
- `walthall`
- `minnaert`
- `ocean`
- `iaquinta_pinty`
- `rahman`
- `kuusk`
- `modis`
- `ross_li_maignan`

Operational rules enforced by the schema:

- `heterogeneous_lambertian` requires `setup.surface.environment`.
- `homogeneous_brdf` requires `setup.surface.brdf`.
- BRDF atmospheric-correction modes require
  `setup.surface.mode = "homogeneous_brdf"`.

Supported atmospheric-correction modes:

- `none`
- `lambertian_reflectance`
- `lambertian_radiance`
- `brdf_reflectance`
- `brdf_radiance`

The default SIAC-facing setting is Lambertian reflectance correction with a
reference reflectance of `0.1`.

Example BRDF configuration:

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.sixs]
mode = "auto"

[algorithms.rt.setup.surface]
mode = "homogeneous_brdf"

[algorithms.rt.setup.surface.target]
kind = "constant"
constant = 0.15

[algorithms.rt.setup.surface.brdf]
model = "rahman"

[algorithms.rt.setup.surface.brdf.parameters]
intensity = 0.12
asymmetry_factor = 0.03
structural_parameter = 0.45

[algorithms.rt.setup.atmospheric_correction]
mode = "brdf_reflectance"
value = 0.1
```

## Route Selection: `direct`, `scene_lut`, `auto`

### `direct`

`direct` evaluates native 6S for every valid scene case.

Use it when:

- geometry varies strongly from pixel to pixel
- the scene is small
- you need exact per-case native evaluation
- a scene-LUT would explode into too many distinct interpolation nodes

### `scene_lut`

`scene_lut` builds a reduced RT grid across the scene axes, runs native 6S on
that reduced case set, then interpolates the requested outputs back to the full
scene.

Use it when:

- the scene is spatially smooth in geometry and atmospheric state
- many pixels collapse onto a small number of representative RT cases
- you want large speedups for broad, compressible scenes

### `auto`

`auto` uses a simple heuristic:

1. never use scene-LUT if `mode = "direct"`
2. always use scene-LUT if `mode = "scene_lut"`
3. in `auto`, only consider scene-LUT when:
   - valid pixel count is at least `scene_lut_min_pixels`
   - the reduced LUT case count is positive and smaller than the direct case count
   - `direct_case_count / lut_case_count >= scene_lut_required_speedup`

Default thresholds:

- `scene_lut_min_pixels = 512`
- `scene_lut_max_nodes_per_axis = 4`
- `scene_lut_max_cases = 4096`
- `scene_lut_required_speedup = 1.5`

## Benchmark Guidance

Measured one-band examples show that route selection depends on how compressible
the scene geometry and atmosphere are:

| Scene | Direct | Scene-LUT | LUT cases vs direct cases | Operational reading |
| --- | ---: | ---: | ---: | --- |
| Smooth Lambertian `16x16` | `5.26 s` | `0.38 s` | `16` vs `256` | scene-LUT is strongly favorable |
| Smooth Lambertian `32x32` | `22.62 s` | `0.44 s` | `16` vs `1024` | scene-LUT scales much better |
| Smooth Rahman BRDF `16x16` | `5.78 s` | `0.48 s` | `16` vs `256` | scene-LUT still favorable on smooth BRDF scenes |
| Mixed-geometry Lambertian `8x8` | `1.27 s` | `91.24 s` | `3888` vs `64` | direct is decisively better |

Accuracy guidance from the same comparisons:

- smooth scenes: scene-LUT stayed very close to direct, with surface-reflectance
  deltas around `1e-5` to `1e-4`
- geometry-rich scenes: scene-LUT was much less favorable operationally, both
  in runtime and approximation quality

Practical guidance:

- prefer `auto` as the default operational mode
- pin `mode = "scene_lut"` only for scenes you already know are smooth and
  compressible
- pin `mode = "direct"` for small scenes, strongly varying geometry, or QA and
  parity work
- benchmark with the same scene, band set, and thread count while changing only
  `algorithms.rt.sixs.mode`

## Parallel Execution

### `parallel_backend = "openmp"`

This is the default and preferred execution path. The native extension loads one
isolated session and executes batched cases with OpenMP threads.

Use this for normal production runs unless you have an environment-specific
reason to avoid OpenMP.

### `parallel_backend = "worker_libraries"`

This fallback mode creates multiple isolated copies of the native extension and
fans out case chunks across Python worker threads.

Use it when:

- OpenMP is unavailable or unstable in your runtime environment
- you need a conservative fallback path for native execution

Relevant knobs:

- `worker_libraries`
- `chunk_size`
- `native_threads`

If the runtime cannot create isolated module copies because of filesystem
permissions, SIAC logs a warning and falls back to the OpenMP path instead of
failing the whole run.

## Output Variables

`output_variables` accepts any of the currently supported native names below.
`xap`, `xbp`, and `xcp` are always retained by SIAC even if you request a
smaller list.

### Common Operational Outputs

These are the most useful outputs for day-to-day SIAC runs:

- `xap`, `xbp`, `xcp`
- `xa`, `xb`, `xc`
- `rapp`, `xrad`, `rog`
- `tgasm`
- `sdtott`, `sutott`, `sttott`
- `sast`

### Full Supported Output Surface

```text
xap, xbp, xcp, xa, xb, xc, rapp, xrad, rog, rogbrdf, refet, alumet, refet1,
refet2, refet3, rpfet, plumet, xpol, rpfet_over_refet, aini_1_1, aini_1_2,
aini_1_3, aini_2_1, aini_2_2, aini_2_3, ainr_1_1, ainr_1_2, ainr_1_3,
ainr_2_1, ainr_2_2, ainr_2_3, sb, seb, dgasm, ugasm, tgasm, sdwava, suwava,
stwava, sdozon, suozon, stozon, sddica, sudica, stdica, sdoxyg, suoxyg, stoxyg,
sdniox, suniox, stniox, sdmeth, sumeth, stmeth, sdmoca, sumoca, stmoca,
sdtotr, sutotr, sttotr, sdtota, sutota, sttota, sdtott, sutott, sttott, sasr,
sasa, sast, sodray, sodaer, sodtot, sodrayp, sodaerp, sodtotp, sroray, sroaer,
srotot, srqray, srqaer, srqtot, sruray, sruaer, srutot, srpray, srpaer, srptot,
sdpray, sdpaer, sdptot, sdppray, sdppaer, sdpptot, fophsr, fophsa, fophst,
foqhsr, foqhsa, foqhst, fouhsr, fouhsa, fouhst, spdpray, spdpaer, spdptot,
pizerr, pizera, pizert, rocave, robar1_over_xnorm1, robar2_over_xnorm2, rbard,
albbrdf, rfoamave, rwatave, rglitave, rooceaw
```

### Example: Restrict Outputs

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.sixs]
output_variables = ["xap", "xbp", "xcp", "tgasm", "sutott", "sast"]
```

That keeps the output payload smaller while preserving the mandatory correction
coefficients.

## Parity Validation

The repo contains two concrete validation paths for the native backend:

- upstream parity harness: `python/siac/sixs_upstream_parity.py`
- integration coverage: `tests/integration/test_sixs_backend_integration.py`

The parity harness builds:

- the original upstream 6SV2.1 executable
- the SIAC native extension in `build_profile = "parity"`

It then compares native outputs against parsed upstream stdout for a fixed case
suite. This is the right validation path when you are checking that native 6S
behavior still matches upstream 6S after backend changes.

The integration tests additionally cover:

- multi-band output consistency
- repeated native-session reuse
- `scene_lut` versus `direct` agreement on compressible scenes
- `worker_libraries` versus `openmp` agreement

## Troubleshooting

### `Unknown 6S output variable`

The requested name is not in the native output surface. Use one of the exact
names listed on this page.

### Native build fails because compiler tools are missing

Install a Fortran compiler, `meson`, and `ninja`, or use the Pixi `rt6s`
environment.

### `scene_lut` is slower than `direct`

That usually means the scene is not compressible enough for interpolation.
Switch to `mode = "direct"` or leave `mode = "auto"` and tune:

- `scene_lut_min_pixels`
- `scene_lut_max_nodes_per_axis`
- `scene_lut_max_cases`
- `scene_lut_required_speedup`

### OpenMP is problematic in the target environment

Switch to:

```toml
[algorithms.rt.sixs]
parallel_backend = "worker_libraries"
worker_libraries = 4
chunk_size = 4096
```

If the worker-library path still cannot create isolated native module copies in
the configured filesystem location, SIAC warns and falls back to the OpenMP
path. In that situation, either use a writable build/module directory or stay
on the default OpenMP backend explicitly.

### Native backend returns some all-`NaN` pixels

The backend masks invalid or failed native cases with `NaN` and logs a warning
when the native status array contains non-zero values. Check the scene inputs
and requested configuration first, then retry with:

- `mode = "direct"`
- `parallel_backend = "openmp"`
- a smaller requested `output_variables` list if you are debugging output scope

## Related Pages

- [Installation](../getting-started/installation.md)
- [Configuration Reference](../reference/configuration.md)
- [Troubleshooting](troubleshooting.md)
- [Native 6SV2.1 Backend Architecture](../architecture/native-sixs-backend.md)
