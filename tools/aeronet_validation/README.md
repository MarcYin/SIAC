# AERONET Validation Experiment

Validates SIAC aerosol retrieval (AOT at 550 nm) against global AERONET
direct-sun measurements, comparing the surface-prior approaches
(`kernel_model`, `whittaker`, `monthly_database`) while every other config
knob (solver, RT backend, priors, cloud masking) stays identical.

Data root: `/gws/ssde/j25a/nceo_isp/public/siac_refactor`

```text
<data-root>/
├── aeronet/          # site list, raw web-service responses, parsed AOD per site
├── matchups/         # per-site S2 search cache + matchups.csv
├── cache/            # SIAC caches shared by all runs (s2, cams, brdf)
├── runs/             # <approach>/<matchup_id>/{config.toml,result.json,fields.nc}
├── results/          # joined CSVs, summary.csv, scatter_<approach>.png
└── slurm/            # generated sbatch script + array logs
```

## Prerequisites

- Pixi environment with the Rust extension built:
  `pixi install && pixi run bootstrap && pixi run build-rust`
- NASA Earthdata credentials for MCD43 (BRDF prior): `EARTHDATA_USERNAME` /
  `EARTHDATA_PASSWORD`. Legacy `Earthdata_user` / `Earthdata_pass` variables
  are mapped automatically.
- CAMS: read directly from the local JASMIN mirror
  `/gws/ssde/j25b/nceo_ard/public/cams/` (`YYYY-MM-DD.nc`, 2015->present) when
  present. `~/.cdsapirc` with an Atmosphere Data Store key is only the
  fallback for dates missing from the mirror — note ADS requests queue and
  can stall a run for tens of minutes.
- No credentials needed for Sentinel-2: catalog search uses anonymous CDSE
  OData, scene download uses the public GCS bucket (`--s2-backend gcs`).
  Switch to `--s2-backend cdse` only if you have `SIAC_CDSE_USERNAME` /
  `SIAC_CDSE_PASSWORD` exported.
- The radiative-transfer LUT is read from the local GWS copy
  (`/gws/ssde/j25a/nceo_isp/public/libradtran_continental_average_lut_1nm.zarr.zip`).

## Running

All commands from the repository root:

```bash
export PYTHONPATH=.:python

# 1. AERONET measurements for the experiment period (resumable)
pixi run python -m tools.aeronet_validation.cli fetch-aeronet \
    --start-date 2024-01-01 --end-date 2024-12-31 --level AOD20

# 2. Sentinel-2 matchups (anonymous CDSE search, +/-30 min window, resumable)
pixi run python -m tools.aeronet_validation.cli matchup \
    --start-date 2024-01-01 --end-date 2024-12-31 --max-cloud-cover 80

# 3. Retrievals run on LOTUS via SLURM (do NOT run retrievals on the
#    interactive sci nodes): manifest -> sbatch array. A seasonally balanced
#    subsample (--per-site-per-month) keeps download/compute volume sane.
pixi run python -m tools.aeronet_validation.cli build-manifest --per-site-per-month 1
pixi run python -m tools.aeronet_validation.cli make-slurm
sbatch /gws/ssde/j25a/nceo_isp/public/siac_refactor/slurm/submit_runs.sbatch
# JASMIN LOTUS note: qos standard/short/long enforce cpu=1 per node; the
# template defaults to --qos=high (per-user cap cpu=576, mem=4500G).

# 4. Score against AERONET
pixi run python -m tools.aeronet_validation.cli compare

# progress at any time
pixi run python -m tools.aeronet_validation.cli status
```

## Protocol

- **Matchup**: S2 acquisition within the search period whose overpass has at
  least one AERONET measurement inside +/-30 min (`--window-minutes`). One
  product per acquisition (tile duplicates and older baselines dropped).
- **Retrieval**: AOI of +/-0.05 deg (~11 km box) centered on the site
  (`--aoi-half-width-deg`), aerosol grid 120 m (`--aerosol-resolution`),
  LUT RT backend, CAMS atmospheric prior, MCD43-based surface priors.
- **Extraction**: nearest-pixel and window statistics (mean/median/std/count)
  of retrieved AOT and TCWV within +/-1.5 km of the site
  (`--extract-radius-m`), plus AOI cloud fraction.
- **AERONET AOD550**: from AOD_500nm scaled by the 440-870 Angstrom exponent,
  falling back to a log-log fit across 440/500/675/870 nm channels.
- **Scoring** (after filtering: AOT valid fraction >= 0.5, AOI cloud fraction
  <= 0.4): N, bias, MAE, RMSE, median |error|, Pearson R, regression
  slope/intercept, and fraction within EE = +/-(0.05 + 0.15 x AOD).

## Notes

- Every stage is resumable: completed sites/searches/runs are skipped unless
  `--refetch`/`--overwrite` is passed.
- `runs/<approach>/<matchup_id>/config.toml` is the exact resolved config
  payload for that task; `traceback.txt` appears next to failed runs.
- Add a new approach by extending `APPROACHES` in
  `tools/aeronet_validation/common.py` with a config overlay (for example a
  `bestpixel` monthly-composite provider once that optional dependency is
  installed).
