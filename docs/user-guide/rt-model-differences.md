# RT Model Differences

This page is a measured comparison of the three RT routes currently exposed in
SIAC:

- native direct 6S
- native 6S `scene_lut`
- the remote ZIP/Zarr LUT at `paths.lut_path`

These routes are not interchangeable in the same sense:

- native `direct` and native `scene_lut` use the same native 6SV2.1 backend
  and differ mainly in how often 6S is evaluated
- the remote ZIP/Zarr LUT is a different RT family in the current public SIAC
  setup: it is a libRadtran-generated spectral LUT with SIAC's coefficient
  derivation layer on top

Because of that, users should separate two questions:

1. Is `scene_lut` an accurate acceleration of native 6S for this scene?
2. How different is the remote libRadtran LUT from native 6S for the same
   atmospheric inputs?

## Executive Summary

- Native `scene_lut` versus native `direct` is usually a routing decision, not
  a model-family decision. On smooth scenes it is much faster and remains very
  close to direct 6S.
- The remote ZIP/Zarr LUT is not “6S, but cached”. It is a different RT model
  family, so non-zero coefficient differences are expected even when geometry,
  AOT, water vapour, ozone, and elevation are matched.
- In the current measured suite, native-versus-remote differences are modest
  for `xap`, larger for `xcp`, and consistently largest for `xbp`.
- The tested NIR band (`B08`) shows the clearest and most stable divergence
  between remote libRadtran LUT and native 6S.
- If you need the closest practical agreement to upstream 6S, use native 6S.
  If you need operational convenience and accept a different RT family, the
  remote LUT remains a valid option.

## What Is Actually Being Compared

### Native 6S Direct

This route evaluates the native 6SV2.1 backend directly at each required scene
point.

### Native 6S Scene-LUT

This route still uses native 6SV2.1, but SIAC first evaluates 6S on a compact
scene-local support grid and interpolates coefficients across the scene.

### Remote ZIP/Zarr LUT

This route uses the remote spectral LUT configured by `paths.lut_path`. In the
current public setup, that LUT is libRadtran-derived and packaged as ZIP/Zarr.
It is therefore best understood as a different RT backend with the same SIAC
coefficient interface, not as a remote cache of 6S outputs.

## Why `xap`, `xbp`, and `xcp` Matter

SIAC applies the atmospheric-correction coefficients through the BOA formula:

```text
y = xap * toa - xbp
boa = y / (1 + xcp * y)
```

That is the reason this page focuses on `xap`, `xbp`, and `xcp` rather than on
internal RT quantities. Differences in these three arrays directly change the
surface-reflectance result used downstream.

In the current experiments, `xbp` is the most sensitive coefficient across RT
families. That does not automatically mean the BOA reflectance difference is of
the same size, but it is a useful early indicator that two RT routes are not
equivalent.

## Experiment Inputs And Methodology

The numbers below come from these reproducible scripts:

- `tools/benchmark_6s_routes.py`
- `tools/compare_6s_route_coefficients.py`
- `tools/compare_native_6s_to_remote_lut.py`

The current checked results were produced from:

- `tmp/6s_route_benchmark_land_smooth.json`
- `tmp/6s_route_benchmark_land_mixed8.json`
- `tmp/6s_route_benchmark_brdf_smooth16.json`
- `tmp/6s_route_coefficient_comparison.json`
- `tmp/rt_model_family_report/report.json`

### Local Native 6S Configuration Used For Cross-Model Comparison

For the native-versus-remote comparison, the local 6S side was configured as:

- `atmospheric_profile = "us_standard_62"`
- `atmospheric_columns_mode = "input_columns"`
- `aerosol_profile = "continental"`
- homogeneous Lambertian surface reflectance `0.1`

That setup is deliberate:

- it keeps the local 6S atmospheric profile shape aligned with the packaged
  US62 assumption used by the remote LUT
- it still lets SIAC inject the scene `tcwv` and `tco3` values into the native
  6S run

This is important for interpretation. The local 6S route is not forced into the
fixed built-in water-vapour and ozone totals of the US62 profile; it preserves
the US62 profile shape while using scene-specific column totals.

### Remote LUT Identity

The current remote ZIP/Zarr LUT report records:

- generator: `libRadtran`
- aerosol family: `continental_average`
- atmospheric profile shape: `us_standard_62`

This already implies a small but real mismatch against local 6S:

- the RT solver family is different
- the aerosol family is close but not identical (`continental_average` versus
  native 6S `continental`)

The experiments below do not try to decompose the total difference into solver,
spectral integration, aerosol family, or interpolation sub-terms. They only
measure the net coefficient difference seen by SIAC.

## Native 6S Direct Versus Native 6S Scene-LUT

This is the safer comparison because both routes share the same 6S solver. The
main question is whether the scene is smooth enough for interpolation to be both
cheap and accurate.

### Timing

| Case | Direct | Scene-LUT | Direct cases | Scene-LUT cases | Practical reading |
| --- | ---: | ---: | ---: | ---: | --- |
| Smooth Lambertian `16x16` | `5.26 s` | `0.38 s` | `256` | `16` | scene-LUT clearly favorable |
| Smooth Lambertian `32x32` | `22.62 s` | `0.44 s` | `1024` | `16` | scene-LUT strongly favorable |
| Smooth Rahman BRDF `16x16` | `5.78 s` | `0.48 s` | `256` | `16` | scene-LUT still favorable |
| Mixed-geometry Lambertian `8x8` | `1.27 s` | `91.24 s` | `64` | `3888` | direct is the only sensible route |

### Coefficient Drift

| Case | Mean relative `xap` error | Mean relative `xbp` error | Mean relative `xcp` error | Max BOA reflectance delta |
| --- | ---: | ---: | ---: | ---: |
| Smooth Lambertian `16x16` | `0.018%` | `0.055%` | `0.043%` | `2.8e-05` |
| Smooth Rahman BRDF `16x16` | `0.017%` | `0.084%` | `0.053%` | `6.7e-05` |
| Mixed-geometry Lambertian `8x8` | `0.075%` | `0.237%` | `0.138%` | `9.9e-05` |

### Interpretation

- If the scene compresses to a small support grid, `scene_lut` is a
  high-confidence acceleration of native 6S.
- If the support grid grows toward the direct case count, `scene_lut` loses its
  advantage and can become slower than `direct`.
- The current `mode = "auto"` heuristic is appropriate as the default route
  selector because it evaluates exactly that compression ratio.

## Native 6S Versus Remote ZIP/Zarr LUT

This is a model-family comparison. The remote LUT is not “6S but precomputed”.
The right question is not whether the coefficients are identical, but how large
the differences are and whether they are acceptable for the intended analysis.

## Case Suite

### Point Cases

| Case | Description | Geometry and atmosphere pattern |
| --- | --- | --- |
| `reference_lowland` | low aerosol, moderate water vapour, low elevation | near-nadir, mild geometry |
| `humid_highland` | higher AOT, higher water vapour, elevated terrain | same basic geometry as reference case |
| `oblique_dry_geometry` | drier atmosphere with more oblique geometry | stronger view/azimuth offset |

### Small Scene Cases

| Case | Description | Pattern |
| --- | --- | --- |
| `reference_lowland_smooth_2x2` | small scene around the lowland reference case | smooth spatial gradients in geometry, AOT, water vapour, ozone, and elevation |
| `oblique_dry_mixed_2x2` | small scene around the oblique dry case | stronger mixed spatial gradients to stress route stability |

### Bands Covered

The current cross-model report uses:

- `B02`
- `B03`
- `B04`
- `B08`

That gives a useful spread from visible to NIR without making the remote lookup
cost unreasonably large during documentation updates.

## Point-Case Results

### Point-Case Relative Difference By Band

The table below summarizes the remote libRadtran LUT versus native 6S
difference across the three point cases. Each cell is shown as
`min / mean / max`.

| Band | `xap` | `xbp` | `xcp` |
| --- | ---: | ---: | ---: |
| `B02` | `0.27 / 1.02 / 1.60 %` | `4.18 / 9.83 / 17.91 %` | `0.84 / 1.12 / 1.40 %` |
| `B03` | `0.65 / 1.32 / 2.62 %` | `3.20 / 9.15 / 16.52 %` | `0.37 / 1.40 / 2.10 %` |
| `B04` | `0.63 / 1.35 / 2.55 %` | `5.58 / 11.23 / 16.83 %` | `1.22 / 2.79 / 4.16 %` |
| `B08` | `2.42 / 2.96 / 3.26 %` | `11.65 / 15.83 / 18.54 %` | `7.45 / 8.65 / 10.47 %` |

Across all twelve point-and-band combinations, the overall mean relative
differences were:

- `xap`: about `1.66%`
- `xbp`: about `11.51%`
- `xcp`: about `3.49%`

### What The Point Cases Show

- `xap` stays relatively close in the visible bands and diverges more in `B08`.
- `xbp` is consistently the most sensitive coefficient, even in the visible.
- `xcp` is fairly close in `B02` to `B04`, but rises sharply in the tested NIR
  case.

## Small-Scene Results

### Scene-Mean Relative Difference By Band

The next table summarizes the remote-versus-native difference over the two small
`2x2` scenes. Each cell is the range of scene-mean relative differences across
the two tested scenes, shown as `min / mean / max`.

| Band | `xap` | `xbp` | `xcp` |
| --- | ---: | ---: | ---: |
| `B02` | `1.06 / 1.39 / 1.73 %` | `3.78 / 12.83 / 21.89 %` | `0.62 / 0.86 / 1.10 %` |
| `B03` | `0.82 / 1.67 / 2.53 %` | `2.99 / 11.73 / 20.46 %` | `0.57 / 0.70 / 0.83 %` |
| `B04` | `1.01 / 1.58 / 2.15 %` | `5.60 / 13.15 / 20.70 %` | `1.16 / 1.77 / 2.38 %` |
| `B08` | `3.35 / 3.53 / 3.72 %` | `11.98 / 17.48 / 22.98 %` | `7.70 / 7.90 / 8.10 %` |

Across all eight scene-and-band combinations, the mean of the scene-mean
relative differences was:

- `xap`: about `2.04%`
- `xbp`: about `13.80%`
- `xcp`: about `2.81%`

### Local Maxima Matter

Scene means do not tell the whole story. In the mixed `2x2` scene, local
cell-level differences were larger than the scene mean:

- `xbp` reached about `29%` max relative difference in the mixed `B08` case
- `xap` reached about `5.36%` max relative difference in the same case
- `xcp` reached about `8.18%` max relative difference there

That does not automatically translate one-to-one into BOA reflectance error, but
it shows that the cross-model difference is not washed out by simple spatial
averaging.

## Worked Example

For the `reference_lowland` point case in `B04`:

- `aot550 = 0.15`
- `tcwv = 2.0`
- `tco3 = 0.30`
- `elevation = 0.1 km`
- `sza / saa / vza / vaa = 30 / 150 / 5 / 110 deg`

the measured coefficients were:

| Route | `xap` | `xbp` | `xcp` |
| --- | ---: | ---: | ---: |
| Native direct 6S | `1.16139` | `0.02812` | `0.07191` |
| Remote libRadtran LUT | `1.15122` | `0.02655` | `0.07104` |

That corresponds to:

- `xap`: `0.88%` relative difference
- `xbp`: `5.58%` relative difference
- `xcp`: `1.22%` relative difference

This is a representative visible-band case: the routes are directionally close,
but not numerically identical.

## Why The Remote LUT And Native 6S Differ

The current experiments support these conclusions:

- The difference is not just interpolation error, because it remains visible in
  direct point comparisons.
- The difference is not removed by matching scene `tcwv` and `tco3`, because it
  persists even when the local 6S route uses the same profile shape plus
  scene-specific column totals.
- The dominant causes are the RT-family change itself and the packaged aerosol
  family mismatch between remote `continental_average` and local native 6S
  `continental`.

What this page does **not** claim:

- it does not isolate exactly how much of the difference comes from solver
  physics versus aerosol-family naming differences
- it does not claim that one route is universally “better”
- it does not claim that coefficient-relative difference equals BOA-relative
  difference

## Runtime Behavior

The remote LUT timing is dominated by remote I/O and cache locality, not by RT
physics. In the measured run behind `tmp/rt_model_family_report/report.json`:

- the first `B02` remote lookup in a process took roughly `32` to `41` seconds
- later `B03`, `B04`, and `B08` lookups in the same process were a few
  milliseconds

Treat those times as infrastructure-dependent and cache-dependent. They should
not be interpreted as the intrinsic cost of libRadtran-style RT evaluation.

## Practical Route Guidance

### Prefer native direct 6S when

- you want the closest practical agreement with upstream 6SV2.1
- you are validating or benchmarking against original 6S
- the scene is too heterogeneous for a compact scene-LUT
- cross-model differences of a few percent in `xap` or much larger differences
  in `xbp` are unacceptable

### Prefer native 6S scene-LUT when

- you want native 6S semantics
- the atmosphere and geometry are spatially smooth
- the scene-LUT planner collapses the workload to a small support grid
- you want a substantial runtime reduction without changing RT family

### Prefer the remote ZIP/Zarr LUT when

- you want an immediately usable remote spectral LUT route
- you accept that the backend is libRadtran-based rather than 6S-based
- the convenience of the hosted LUT matters more than exact agreement with
  native 6S

The remote ZIP/Zarr LUT should not be described as “the same as local 6S but
faster”. The accurate description is:

> a libRadtran-based remote RT model with SIAC’s coefficient derivation layer

## Reproducing The Current Reports

Route timing and native direct-versus-scene-LUT coefficient comparison:

```bash
PYTHONPATH=python pixi run -e rt6s python tools/benchmark_6s_routes.py
PYTHONPATH=python pixi run -e rt6s python tools/compare_6s_route_coefficients.py
```

Native 6S versus remote ZIP/Zarr LUT:

```bash
PYTHONPATH=python pixi run -e rt6s python tools/compare_native_6s_to_remote_lut.py \
  --source-dir tmp/6s_upstream \
  --build-dir tmp/rt_model_family_report \
  --module-path tmp/rt6s_profile_columns_smoke/_siac_rt6s_native.cpython-311-darwin.so \
  --output tmp/rt_model_family_report/report.json
```

If the remote report is being regenerated on a cold machine or in a new
process, expect the first remote lookup to dominate wall-clock time.

## Bottom Line

- Use native 6S when you need 6S behavior.
- Use native `scene_lut` when you want to keep 6S behavior and the scene is
  smooth enough to compress well.
- Use the remote ZIP/Zarr LUT when you want a different, convenient RT route
  and the observed coefficient differences are acceptable for your use case.
