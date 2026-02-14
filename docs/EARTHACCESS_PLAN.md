# Earthaccess Integration Plan (Pre-Implementation)

## 1. Purpose

This plan defines how SIAC will integrate NASA Earthdata access through
`earthaccess` before implementation starts.

Primary goals:

- Add robust Earthdata access for Landsat, atmospheric priors, and BRDF products.
- Keep Earthaccess logic isolated in I/O/provider modules.
- Preserve existing module contracts (`AtmosphericState`, `SurfacePrior`, etc.).

This is a planning document only; no runtime behavior is changed by this file.

---

## 2. Scope

### In scope

- Landsat data access path (search + access/download) for M1 ingestion.
- Atmospheric prior path using NASA products, including MCD19 AOD.
- BRDF product paths for MODIS MCD43 and VIIRS VNP43 equivalents.
- Auth/cache/query architecture for Earthaccess.
- Test strategy, rollout stages, and acceptance criteria.

### Out of scope (this phase)

- Rewriting solver/corrector math.
- Replacing existing CAMS path.
- Adding non-NASA providers (except existing providers already in code).

---

## 3. Earthaccess API Baseline

Planned Earthaccess primitives to standardize on:

- `earthaccess.login(...)`
- `earthaccess.search_datasets(...)`
- `earthaccess.search_data(...)`
- `earthaccess.open(...)`
- `earthaccess.download(...)`

Implementation rule: SIAC modules should not call these functions directly.
All direct interactions go through a single wrapper module.

---

## 4. Data Product Mapping

The exact CMR `short_name` values must be discovered and validated in code during
M0 using `search_datasets(...)` before defaults are finalized.

| SIAC Need | Candidate Product Family | Planned Use |
|---|---|---|
| Landsat scene ingestion (M1) | Landsat C2 / Earthdata-available Landsat families | Search by AOI/time, access/download granules for preprocessor input |
| Atmospheric prior AOD (M2) | MCD19 AOD family | Populate AOT prior in `AtmosphericState`; combine with TCWV/TCO3 strategy |
| Surface prior BRDF (M3) | MCD43 BRDF family | BRDF kernels/parameters for MODIS-based path |
| Surface prior BRDF (M3, VIIRS equivalent) | VNP43 BRDF family | VIIRS-based BRDF alternative or fallback |

Open technical decision to close in M0:

- For M2, whether TCWV/TCO3 come from the same Earthaccess product family,
  another Earthaccess family, or a hybrid fallback (e.g., CAMS).

---

## 5. Target Module Architecture

### 5.1 Earthaccess Wrapper (single integration boundary)

Planned module: `python/siac/io/earthaccess_source.py`

Responsibilities:

- Lazy auth (`login` only when first network call is needed).
- Dataset discovery helpers for startup validation.
- AOI/time query normalization (`bounds + crs + datetime -> CMR query args`).
- Access policy abstraction:
  - stream via `open(...)`
  - localize via `download(...)`
- Retry/error normalization (auth failures, empty search results, transient errors).

### 5.2 Product Registry

Planned module: `python/siac/io/earthaccess_catalog.py`

Responsibilities:

- Central mapping from logical product keys to discovered/approved CMR metadata.
- Version pin policy for product family defaults.
- Optional provider override by config.

### 5.3 Provider Modules

Planned modules:

- `python/siac/priors/atmospheric/mcd19_earthaccess.py`
- `python/siac/priors/atmospheric/merra2.py` (or equivalent Earthaccess atmo source)
- `python/siac/priors/brdf/mcd43_earthaccess.py`
- `python/siac/priors/brdf/vnp43_earthaccess.py`

Responsibilities:

- Convert Earthdata granules to SIAC contract outputs.
- Handle reprojection/crop to AOI.
- Keep provider-specific parsing localized.

### 5.4 Landsat Data Source

Planned modules:

- `python/siac/io/landsat_data_source.py`
- `python/siac/satellite/landsat.py`

Responsibilities:

- Discover Landsat scenes by AOI + time window.
- Resolve acquisition/local access strategy.
- Convert to `ObservationBundle` using the same M1 contract as S2.

---

## 6. Configuration Plan

Planned config additions (shape, exact names may vary):

```yaml
nasa:
  earthaccess:
    enabled: true
    strategy: download   # download | open
    cache_dir: /path/to/cache
    provider: null       # optional CMR provider constraint

atmo_prior:
  provider: mcd19        # cams | merra2 | mcd19 | user
  mcd19:
    short_name: null     # resolved/validated default

brdf:
  provider: mcd43        # mcd43 | vnp43 | gee | zarr | user
  mcd43:
    short_name: null
  vnp43:
    short_name: null

satellite:
  landsat:
    enabled: true
    access_strategy: earthaccess
```

Design rule: product identifiers should be overrideable in config, with defaults
coming from the validated product registry.

---

## 7. Rollout Phases

### Phase M0: Discovery + foundation

Deliverables:

- `earthaccess_source.py` wrapper skeleton with lazy auth.
- `earthaccess_catalog.py` and product validation helper.
- A one-time script/test that verifies required product families are discoverable.

Acceptance:

- Can authenticate and perform `search_datasets` + `search_data` for target AOI/time.
- Product keys resolve to at least one valid dataset.

### Phase M1: Landsat ingestion path

Deliverables:

- `landsat_data_source.py` query/access path.
- `satellite/landsat.py` preprocessor returning valid `ObservationBundle`.

Acceptance:

- AOI/time query returns deterministic scene selection.
- Pipeline can run M1 for Landsat inputs without bypass hacks.

### Phase M2: Atmospheric priors via Earthaccess

Deliverables:

- MCD19-based atmospheric provider.
- Explicit TCWV/TCO3 strategy (Earthaccess source or fallback policy).

Acceptance:

- Provider returns valid `AtmosphericState` on AOI.
- Uncertainties and units documented and tested.

### Phase M3: BRDF providers (MCD43 + VNP43)

Deliverables:

- MCD43 provider module.
- VNP43 provider module.
- Shared BRDF parser/normalizer utilities if needed.

Acceptance:

- Both providers return BRDF-compatible inputs for M3.
- AOI/time window behavior is reproducible and tested.

### Phase M4: Pipeline + docs hardening

Deliverables:

- Config wiring in resolver paths.
- End-to-end examples for Landsat + Earthaccess providers.
- Updated architecture docs and test plan coverage.

Acceptance:

- All targeted workflows are selectable by config only.
- Existing S2/CAMS baseline remains green.

---

## 8. Testing Strategy

### 8.1 Unit tests

- Wrapper tests with mocked Earthaccess calls.
- Query normalization tests (AOI, temporal windows).
- Provider parsing tests with fixture granule metadata/files.
- Contract validation tests for `AtmosphericState` and `SurfacePrior`.

### 8.2 Integration tests

- Network-marked tests for one product from each family:
  - Landsat candidate
  - MCD19 candidate
  - MCD43 candidate
  - VNP43 candidate
- Skip cleanly when auth/network unavailable.

### 8.3 Regression tests

- Existing Sentinel-2 + CAMS paths must remain unchanged.
- Existing LUT/emulator backends unaffected.

---

## 9. Operational Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Auth/session failures | Pipeline startup failure | Lazy auth + explicit error messages + retries |
| Product naming/version drift | Provider breakage | Discovery-based validation in M0; registry-driven defaults |
| Granule heterogeneity | Parser instability | Per-product parser adapters + fixture-backed tests |
| Large AOI/time queries | Excess data and latency | Mandatory AOI/time constraints; pagination and limits |
| Remote instability | Flaky tests | Mark network tests and keep deterministic local unit fixtures |

---

## 10. Definition of Done

The Earthaccess initiative is complete when:

- Landsat, MCD19, MCD43, and VNP43 paths are all selectable by config.
- All providers return contract-valid outputs on AOI-scoped queries.
- Product discovery validation is automated (not manual-only).
- Unit/integration/regression test suites are green.
- Docs include runnable configuration examples and troubleshooting notes.
