"""Tests for the libRadtran (uvspec) RT backend.

The real ``uvspec`` engine is an optional, compiled-from-source dependency that
reaches out to disk/data at call time, so these tests exercise the pure pieces
(deck builder, stdout parser, unit conversions, config/rt-setup wiring, build
helpers) and the ``ZarrLUTBackend`` reuse path via an in-memory synthetic LUT —
never the engine itself. An ``@pytest.mark.integration`` test that actually runs
uvspec lives separately and is skip-guarded.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.rt.direct.libradtran import LibRadtranBackend
from siac.algorithms.rt.direct.libradtran_build import (
    LibRadtranPaths,
    _archive_name_from_url,
    _configure_args,
    _reptran_resolution,
    _required_reptran_resolutions,
    ensure_libradtran,
)
from siac.algorithms.rt.direct.libradtran_runner import (
    _DEFAULT_MAX_WORKERS,
    DEEP_WATER_H2O_BANDS_NM,
    LibRadtranRunner,
    _axis_from_scene,
    _InMemorySpectralLUTBackend,
    build_uvspec_deck,
    parse_uvspec_table,
)
from siac.config.algorithms import LibRadtranAlgorithmConfig, RTAlgorithmConfig
from siac.domain.sensors import SensorBand
from siac.rt_setup import resolve_backend_rt_setup, resolve_effective_rt_setup
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients


def _da(vals: list[list[float]]) -> xr.DataArray:
    return xr.DataArray(np.asarray(vals, dtype=np.float64), dims=("y", "x"))


def _scene() -> tuple[GeometryAngles, AtmosphericState]:
    geom = GeometryAngles.from_degrees(
        _da([[30.0, 31.0], [30.5, 30.2]]),
        _da([[150.0, 150.0], [150.0, 150.0]]),
        _da([[5.0, 5.0], [5.0, 5.0]]),
        _da([[110.0, 110.0], [110.0, 110.0]]),
    )
    zeros = _da([[0.0, 0.0], [0.0, 0.0]])
    atmo = AtmosphericState(
        aot=_da([[0.15, 0.25], [0.20, 0.30]]),
        tcwv=_da([[1.2, 1.6], [1.4, 1.5]]),
        tco3=_da([[0.30, 0.30], [0.30, 0.30]]),
        aot_unc=zeros,
        tcwv_unc=zeros,
        tco3_unc=zeros,
        elevation=_da([[0.05, 0.05], [0.05, 0.05]]),
    )
    return geom, atmo


def _band() -> SensorBand:
    return SensorBand(
        name="B02", center_wavelength=492.0, bandwidth=66.0, resolution=10.0, band_index=1
    )


# --------------------------------------------------------------------------- #
# Backend with an injected fake runner (no engine)
# --------------------------------------------------------------------------- #


class _FakeRunner:
    """A runner stand-in that returns a known coefficient dict."""

    def __init__(self) -> None:
        self.calls = 0
        self.observation_time = None

    def set_observation_time(self, observation_time: object) -> None:
        self.observation_time = observation_time

    def compute_coefficients(self, *, geometry, atmo_state, band):  # noqa: ANN001, ARG002
        self.calls += 1
        ones = xr.DataArray(np.full((2, 2), 1.3, dtype=np.float32), dims=("y", "x"))
        return {
            "xap": ones,
            "xbp": ones * 0.05,
            "xcp": ones * 0.1,
            "spherical_albedo": ones * 0.1,  # extra -> goes into extras
        }


def test_backend_coerces_runner_dict_to_rt_coefficients() -> None:
    runner = _FakeRunner()
    backend = LibRadtranBackend(libradtran_config=LibRadtranAlgorithmConfig(), runner=runner)
    assert backend.backend_name == "libradtran"
    assert backend.supports_jacobian() is False
    assert backend.is_available_for_sensor("MSI", "S2A") is True

    geom, atmo = _scene()
    coeffs = backend.compute_coefficients(geom, atmo, _band())
    assert isinstance(coeffs, RTCoefficients)
    assert runner.calls == 1
    np.testing.assert_allclose(coeffs.xap.values, 1.3, rtol=1e-6)
    assert "spherical_albedo" in coeffs.extras
    assert "xap" not in coeffs.extras


def test_backend_set_observation_time_and_jacobian_guard() -> None:
    runner = _FakeRunner()
    backend = LibRadtranBackend(libradtran_config=LibRadtranAlgorithmConfig(), runner=runner)
    backend.set_observation_time("2026-03-29")
    assert runner.observation_time == "2026-03-29"
    geom, atmo = _scene()
    with pytest.raises(NotImplementedError, match="Jacobians"):
        backend.compute_coefficients(geom, atmo, _band(), compute_jacobian=True)


# --------------------------------------------------------------------------- #
# Deck builder (pure)
# --------------------------------------------------------------------------- #


def test_build_uvspec_deck_unit_conversions_and_keywords() -> None:
    deck = build_uvspec_deck(
        data_dir="/data",
        aot550=0.2,
        tcwv_cm=1.4,
        tco3_atmcm=0.3,
        elevation_km=0.05,
        sza_deg=30.0,
        vza_deg=0.0,
        raa_deg=60.0,
        albedo=0.15,
        wavelength_min_nm=400.0,
        wavelength_max_nm=2500.0,
        mol_abs_param="reptran",
        number_of_streams=16,
    )
    # cm -> mm precipitable water; atm-cm -> Dobson.
    assert "mol_modify H2O 14.000000 MM" in deck
    assert "mol_modify O3 300.000000 DU" in deck
    # umu = cos(0) = 1 at nadir.
    assert "umu 1.000000" in deck
    assert "output_quantity reflectivity" in deck
    assert "output_user lambda eglo uu" in deck
    assert "aerosol_species_file continental_average" in deck
    assert "aerosol_set_tau_at_wvl 550 0.200000" in deck
    assert "rte_solver disort" in deck
    assert "number_of_streams 16" in deck


def test_build_uvspec_deck_umu_clamped_off_grazing() -> None:
    deck = build_uvspec_deck(
        data_dir="/d",
        aot550=0.1,
        tcwv_cm=1.0,
        tco3_atmcm=0.3,
        elevation_km=0.0,
        sza_deg=30.0,
        vza_deg=89.99,
        raa_deg=0.0,
        albedo=0.5,
        wavelength_min_nm=400.0,
        wavelength_max_nm=900.0,
        mol_abs_param="reptran",
        number_of_streams=8,
    )
    umu_line = next(line for line in deck.splitlines() if line.startswith("umu "))
    assert float(umu_line.split()[1]) > 0.0


# --------------------------------------------------------------------------- #
# Output parser (pure)
# --------------------------------------------------------------------------- #


def test_parse_uvspec_table_splits_surface_and_toa() -> None:
    # Two rows per wavelength: surface (zout sur) then TOA (zout toa).
    stdout = (
        "  489.000  0.876906  0.131536\n"
        "  489.000  1.000000  0.184093\n"
        "  490.000  0.878184  0.131728\n"
        "  490.000  1.000000  0.183962\n"
    )
    wl, eg_sur, toa_uu = parse_uvspec_table(stdout)
    np.testing.assert_allclose(wl, [489.0, 490.0])
    np.testing.assert_allclose(eg_sur, [0.876906, 0.878184])  # surface eglo
    np.testing.assert_allclose(toa_uu, [0.184093, 0.183962])  # TOA reflectance


def test_parse_uvspec_table_rejects_odd_rows() -> None:
    with pytest.raises(ValueError, match="even count"):
        parse_uvspec_table("  489.000  0.87  0.13\n")


# --------------------------------------------------------------------------- #
# ZarrLUTBackend reuse via in-memory synthetic LUT (no engine)
# --------------------------------------------------------------------------- #


def _synthetic_spectral_lut() -> xr.Dataset:
    dims = ("sza", "vza", "raa", "ozone", "altitude", "aot", "tcwv", "wavelength")
    aot = np.linspace(0.05, 0.5, 4, dtype=np.float32)
    tcwv = np.linspace(0.5, 3.0, 4, dtype=np.float32)
    wl = np.linspace(400.0, 900.0, 51, dtype=np.float32)
    aot_b = aot[:, None, None]
    tcwv_b = tcwv[None, :, None]
    # Physically-ordered terms: brighter surface -> higher TOA; Eg ~ transmittance.
    toa1 = (0.08 + 0.25 * aot_b + 0.0 * tcwv_b + 0.0 * wl).astype(np.float32)
    toa2 = (0.34 + 0.22 * aot_b + 0.0 * tcwv_b).astype(np.float32)
    eg1 = (0.86 - 0.15 * aot_b - 0.01 * tcwv_b + 0.0 * wl).astype(np.float32)
    eg2 = (0.90 - 0.12 * aot_b - 0.01 * tcwv_b).astype(np.float32)

    def _cube(block: np.ndarray) -> np.ndarray:
        block = np.broadcast_to(block, (4, 4, 51)).astype(np.float32)
        return block[None, None, None, None, None, :, :, :]

    return xr.Dataset(
        {
            "TOA_rho1": (dims, _cube(toa1)),
            "TOA_rho2": (dims, _cube(toa2)),
            "Eg_rho1": (dims, _cube(eg1)),
            "Eg_rho2": (dims, _cube(eg2)),
        },
        coords={
            "sza": np.array([30.0], dtype=np.float32),
            "vza": np.array([5.0], dtype=np.float32),
            "raa": np.array([40.0], dtype=np.float32),
            "ozone": np.array([300.0], dtype=np.float32),
            "altitude": np.array([0.05], dtype=np.float32),
            "aot": aot,
            "tcwv": tcwv,
            "wavelength": wl,
        },
        attrs={"rho1": 0.15, "rho2": 0.5},
    )


def test_in_memory_spectral_backend_returns_finite_coefficients() -> None:
    backend = _InMemorySpectralLUTBackend(_synthetic_spectral_lut())
    geom, atmo = _scene()
    coeffs = backend.compute_coefficients(geom, atmo, _band())
    assert isinstance(coeffs, RTCoefficients)
    assert coeffs.xap.shape == (2, 2)
    assert np.all(np.isfinite(coeffs.xap.values))
    assert np.all(coeffs.xap.values > 0.0)  # xap = 1 / T_total
    assert np.all(np.isfinite(coeffs.xbp.values))
    assert np.all(np.isfinite(coeffs.xcp.values))


def test_in_memory_backend_does_not_use_disk_scene_subset_cache(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Regression: distinct in-memory scene LUTs must not collide on the on-disk
    scene-subset cache.

    That disk cache keys on ``lut_path`` + grid-snapped scene coords, but for an
    in-memory LUT ``lut_path`` is a constant ``"<in-memory-libradtran>"``
    placeholder - NOT content-identifying. With the cache enabled, two backends
    fed different datasets but queried on the same scene resolve to the same disk
    path, so the second loads the first's STALE subset and returns its
    coefficients - silently freezing/aliasing the RT (observed as a libRadtran
    AOT retrieval whose aot interpolation collapsed). The in-memory backend now
    disables that disk cache (the within-process subset cache still serves every
    band of the scene, and the costly uvspec runs are cached at the run-cache
    layer). It exists only to skip re-fetching the REMOTE LUT's HTTP chunks.
    """
    monkeypatch.setenv("SIAC_RT_RUN_CACHE_ROOT", str(tmp_path))

    ds_a = _synthetic_spectral_lut()
    ds_b = _synthetic_spectral_lut()
    # Same coords (=> identical scene key) but markedly different RT: brighter
    # path radiance (drives xbp/xcp) and lower transmittance (drives xap).
    ds_b["TOA_rho1"] = ds_b["TOA_rho1"] + 0.05
    ds_b["TOA_rho2"] = ds_b["TOA_rho2"] + 0.05
    ds_b["Eg_rho1"] = ds_b["Eg_rho1"] - 0.08
    ds_b["Eg_rho2"] = ds_b["Eg_rho2"] - 0.08

    backend_a = _InMemorySpectralLUTBackend(ds_a)
    backend_b = _InMemorySpectralLUTBackend(ds_b)
    assert backend_a._scene_cache_dir is None  # disk subset cache disabled
    assert backend_b._scene_cache_dir is None

    geom, atmo = _scene()
    band = _band()
    coeffs_a = backend_a.compute_coefficients(geom, atmo, band)
    coeffs_b = backend_b.compute_coefficients(geom, atmo, band)

    # Distinct RT content => distinct coefficients (no cross-backend cache bleed).
    # Before the fix, backend_b loaded backend_a's stale disk subset and returned
    # backend_a's coefficients, so every term matched.
    assert not np.allclose(coeffs_a.xap.values, coeffs_b.xap.values)
    assert not np.allclose(coeffs_a.xbp.values, coeffs_b.xbp.values)
    # And the in-memory backends wrote no subset into the shared cache root.
    assert not list(tmp_path.rglob("*.subset.nc"))


# --------------------------------------------------------------------------- #
# RT setup + config wiring
# --------------------------------------------------------------------------- #


def test_rt_setup_permits_continental_average_for_libradtran() -> None:
    setup = resolve_backend_rt_setup("libradtran", None)
    assert setup.aerosol is not None
    assert setup.aerosol.profile == "continental_average"
    assert setup.atmosphere is not None
    assert setup.atmosphere.profile == "us_standard_62"
    # And via the RTAlgorithmConfig dispatch.
    eff = resolve_effective_rt_setup(RTAlgorithmConfig(backend="libradtran"), "libradtran")
    assert eff.aerosol is not None and eff.aerosol.profile == "continental_average"


def test_libradtran_config_validators() -> None:
    with pytest.raises(ValueError, match="even"):
        LibRadtranAlgorithmConfig(number_of_streams=15)
    with pytest.raises(ValueError, match="rho1"):
        LibRadtranAlgorithmConfig(rho1=0.6, rho2=0.5)
    with pytest.raises(ValueError, match="wavelength_max"):
        LibRadtranAlgorithmConfig(wavelength_min_nm=900.0, wavelength_max_nm=800.0)


def test_libradtran_requires_https_urls() -> None:
    # The build harness downloads + compiles these; non-https is rejected at load.
    with pytest.raises(ValueError, match="https"):
        LibRadtranAlgorithmConfig(source_url="http://example.org/libRadtran-2.0.6.tar.gz")
    with pytest.raises(ValueError, match="https"):
        LibRadtranAlgorithmConfig(reptran_url="ftp://example.org/reptran.tar.gz")


def test_axis_from_scene_never_negative_for_constant_zero_field() -> None:
    # A uniform aot=0 (or flat prior) must NOT pad into a negative grid node:
    # uvspec rejects negative aerosol optical depth / column water.
    zeros = xr.DataArray(np.zeros((2, 2), dtype=np.float64), dims=("y", "x"))
    axis = _axis_from_scene(zeros, 3)
    assert axis.min() >= 0.0
    assert axis[-1] > axis[0]  # still strictly increasing


# --------------------------------------------------------------------------- #
# Build harness helpers (pure)
# --------------------------------------------------------------------------- #


def test_archive_name_from_url_handles_plain_and_fetchphp() -> None:
    assert (
        _archive_name_from_url("https://www.libradtran.org/download/libRadtran-2.0.6.tar.gz")
        == "libRadtran-2.0.6.tar.gz"
    )
    assert (
        _archive_name_from_url(
            "https://www.libradtran.org/lib/exe/fetch.php?media=download:optprop_v2.1.tar.gz"
        )
        == "optprop_v2.1.tar.gz"
    )
    # A crafted URL cannot escape the cache dir: the result is a bare basename.
    assert (
        _archive_name_from_url("https://x/lib/fetch.php?media=d:../../../../etc/passwd") == "passwd"
    )
    assert _archive_name_from_url("https://x/a/b/../../evil.tar.gz") == "evil.tar.gz"
    with pytest.raises(ValueError, match="safe archive filename"):
        _archive_name_from_url("https://www.libradtran.org/download/")


def test_configure_args_contains_netcdf_prefix() -> None:
    args = _configure_args()
    assert any(a.startswith("--with-netcdf4=") for a in args)


def test_reptran_resolution_gates_only_non_coarse() -> None:
    # fine/medium need the separate reptran archive...
    assert _reptran_resolution("reptran fine") == "fine"
    assert _reptran_resolution("reptran medium") == "medium"
    # ...coarse + plain reptran are bundled, and non-reptran models need nothing.
    assert _reptran_resolution("reptran") is None
    assert _reptran_resolution("reptran coarse") is None
    assert _reptran_resolution("lowtran") is None
    assert _reptran_resolution("crs") is None


def test_ensure_libradtran_raises_when_unavailable_and_no_autobuild(tmp_path) -> None:  # noqa: ANN001
    cfg = LibRadtranAlgorithmConfig(build_dir=tmp_path / "empty", auto_build=False)
    with pytest.raises(RuntimeError, match="auto_build is disabled"):
        ensure_libradtran(cfg)


# --------------------------------------------------------------------------- #
# Concurrency / memory cap (the >100 GB OOM fix)
# --------------------------------------------------------------------------- #


def _patch_available_ram(monkeypatch: pytest.MonkeyPatch, gib: float) -> None:
    import types as _types

    import psutil

    monkeypatch.setattr(
        psutil, "virtual_memory", lambda: _types.SimpleNamespace(available=gib * 1e9)
    )


def test_estimate_uvspec_gb_fine_is_multi_gb() -> None:
    fine = LibRadtranRunner(
        libradtran_config=LibRadtranAlgorithmConfig(
            mol_abs_param="reptran fine", wavelength_min_nm=340.0, wavelength_max_nm=2500.0
        )
    )
    coarse = LibRadtranRunner(
        libradtran_config=LibRadtranAlgorithmConfig(
            mol_abs_param="reptran", wavelength_min_nm=440.0, wavelength_max_nm=520.0
        )
    )
    assert fine._estimate_uvspec_gb() > 5.0  # fine over the full window is heavy
    assert coarse._estimate_uvspec_gb() < fine._estimate_uvspec_gb()


def test_worker_count_never_scales_with_cores_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    # Huge RAM + (implicitly) many cores must still NOT spawn many workers by
    # default — that core-count scaling was the original >100 GB blow-up.
    _patch_available_ram(monkeypatch, 1024.0)
    runner = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig())
    assert runner._resolve_worker_count(100) <= _DEFAULT_MAX_WORKERS


def test_worker_count_caps_explicit_threads_by_memory(monkeypatch: pytest.MonkeyPatch) -> None:
    runner = LibRadtranRunner(
        libradtran_config=LibRadtranAlgorithmConfig(
            mol_abs_param="reptran fine",
            wavelength_min_nm=340.0,
            wavelength_max_nm=2500.0,
            native_threads=12,  # what a 12-core box would have used -> OOM
        )
    )
    # 64 GB available with ~13 GB/process must cap well under the requested 12.
    _patch_available_ram(monkeypatch, 64.0)
    capped = runner._resolve_worker_count(100)
    assert 1 <= capped < 12
    # A tiny machine collapses to a single process.
    _patch_available_ram(monkeypatch, 4.0)
    assert runner._resolve_worker_count(100) == 1


def test_memory_budget_gb_is_hard_ceiling(monkeypatch: pytest.MonkeyPatch) -> None:
    # Even on a huge-RAM box, total concurrent memory stays within memory_budget_gb.
    # This is the >120 GB fix: 3 x ~45 GB fine-full-window processes used to slip
    # through because 0.7 x (huge host RAM) clamped only at _DEFAULT_MAX_WORKERS.
    _patch_available_ram(monkeypatch, 1024.0)
    runner = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig(memory_budget_gb=30.0))
    assert runner._resolve_worker_count(100, per_proc_gb=9.0) == 3  # floor(30/9)
    assert runner._resolve_worker_count(100, per_proc_gb=16.0) == 1  # floor(30/16) -> budget binds
    assert runner._resolve_worker_count(100, per_proc_gb=50.0) == 1  # one proc > budget -> 1


def test_slurm_mem_limit_caps_budget_below_host_ram(monkeypatch: pytest.MonkeyPatch) -> None:
    # psutil reports a huge host, but a SLURM --mem allocation must bound the pool
    # (psutil ignores cgroup/SLURM limits -> the cluster-node over-subscription bug).
    _patch_available_ram(monkeypatch, 1024.0)
    monkeypatch.setenv("SLURM_MEM_PER_NODE", str(20 * 1024))  # 20 GB, reported in MB
    runner = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig(memory_budget_gb=None))
    # effective available = min(1024, 20) = 20 -> budget 0.7*20=14 -> floor(14/9)=1.
    assert runner._resolve_worker_count(100, per_proc_gb=9.0) == 1


# --------------------------------------------------------------------------- #
# Per-region resolution (mol_abs_regions): partition, stitching, memory, guards
# --------------------------------------------------------------------------- #


def test_segments_adaptive_default_adds_fine_water_bands() -> None:
    # Default behaviour (adaptive_deep_water_fine=True, no explicit regions):
    # a cheap base in the windows, reptran fine only in the deep H2O bands.
    runner = LibRadtranRunner(
        libradtran_config=LibRadtranAlgorithmConfig(
            mol_abs_param="reptran medium",
            wavelength_min_nm=400.0,
            wavelength_max_nm=2300.0,
        )
    )
    segs = runner._segments()
    fine_segs = [(lo, hi) for lo, hi, m in segs if "fine" in m]
    expected = [(lo, hi) for lo, hi in DEEP_WATER_H2O_BANDS_NM if lo > 400.0 and hi < 2300.0]
    assert fine_segs == expected  # every in-window deep-water band runs fine
    assert all("medium" in m for _lo, _hi, m in segs if "fine" not in m)  # base elsewhere
    assert segs[0][0] == 400.0 and segs[-1][1] == 2300.0  # spans the window
    for a, b in zip(segs, segs[1:]):
        assert a[1] == b[0]  # contiguous, no gap/overlap


def test_segments_adaptive_disabled_is_single_window() -> None:
    runner = LibRadtranRunner(
        libradtran_config=LibRadtranAlgorithmConfig(
            mol_abs_param="reptran",
            wavelength_min_nm=400.0,
            wavelength_max_nm=2300.0,
            adaptive_deep_water_fine=False,
        )
    )
    assert runner._segments() == [(400.0, 2300.0, "reptran")]


def test_required_reptran_resolutions_union_of_base_and_regions() -> None:
    # Adaptive default = medium base + fine water bands -> the build must fetch BOTH.
    assert _required_reptran_resolutions(
        LibRadtranAlgorithmConfig(mol_abs_param="reptran medium")
    ) == ["fine", "medium"]
    # Coarse base + explicit fine region -> only fine is non-bundled.
    assert _required_reptran_resolutions(
        LibRadtranAlgorithmConfig(
            mol_abs_param="reptran", mol_abs_regions=[(915.0, 985.0, "reptran fine")]
        )
    ) == ["fine"]
    # Coarse base, adaptive off, no regions -> nothing extra (coarse is bundled).
    assert (
        _required_reptran_resolutions(
            LibRadtranAlgorithmConfig(mol_abs_param="reptran", adaptive_deep_water_fine=False)
        )
        == []
    )


def test_segments_partition_is_contiguous_and_non_overlapping() -> None:
    cfg = LibRadtranAlgorithmConfig(
        mol_abs_param="reptran",
        wavelength_min_nm=400.0,
        wavelength_max_nm=2300.0,
        mol_abs_regions=[(890.0, 985.0, "reptran fine"), (1085.0, 1180.0, "reptran fine")],
    )
    segs = LibRadtranRunner(libradtran_config=cfg)._segments()
    # base/fine/base/fine/base, edges meet exactly, models alternate.
    assert [s[2].split()[-1] for s in segs] == ["reptran", "fine", "reptran", "fine", "reptran"]
    assert segs[0][0] == 400.0 and segs[-1][1] == 2300.0
    for a, b in zip(segs, segs[1:]):
        assert a[1] == b[0]  # contiguous, no gap/overlap


def test_estimate_narrow_fine_segment_is_far_cheaper_than_full_fine() -> None:
    # The whole point: confining fine to a ~95 nm deep-water band is ~7-9 GB,
    # not the ~45-55 GB of a full-range fine run.
    narrow = LibRadtranRunner._estimate_uvspec_gb_for("reptran fine", 95.0)
    full = LibRadtranRunner._estimate_uvspec_gb_for("reptran fine", 1900.0)
    assert narrow < 10.0 < 40.0 < full


def test_build_scene_lut_stitches_per_segment_values(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = LibRadtranAlgorithmConfig(
        mol_abs_param="reptran",
        wavelength_min_nm=400.0,
        wavelength_max_nm=700.0,
        wavelength_step_nm=1.0,
        scene_lut_max_nodes_per_axis=2,
        scene_lut_max_cases=8,
        mol_abs_regions=[(500.0, 600.0, "reptran fine")],
    )
    runner = LibRadtranRunner(libradtran_config=cfg)
    monkeypatch.setattr(
        runner,
        "_ensure_paths",
        lambda: LibRadtranPaths(uvspec=Path("/x/uvspec"), data_dir=Path("/x/data")),
    )

    def _fake_uvspec(_paths, deck: str):
        wl_line = next(line for line in deck.splitlines() if line.startswith("wavelength "))
        lo, hi = (float(x) for x in wl_line.split()[1:3])
        is_fine = "fine" in next(
            line for line in deck.splitlines() if line.startswith("mol_abs_param ")
        )
        wl = np.arange(lo, hi + 1.0, 1.0)
        toa = np.full(wl.shape, 0.2 if is_fine else 0.1)  # value encodes which model ran
        eg = np.full(wl.shape, 0.8)
        return wl, eg, toa

    monkeypatch.setattr(runner, "_run_uvspec", _fake_uvspec)

    geom, atmo = _scene()
    ds = runner._build_scene_lut(geom, atmo)
    wl = ds["wavelength"].values
    toa1 = ds["TOA_rho1"].values  # (1,1,1,1,1,aot,tcwv,wl)
    fine = (wl >= 500.0) & (wl < 600.0)
    base = ~fine
    assert np.allclose(toa1[..., fine], 0.2)  # fine segment placed correctly
    assert np.allclose(toa1[..., base], 0.1)  # base (coarse) segments placed correctly


def test_build_scene_lut_no_nan_holes_at_segment_edges(monkeypatch: pytest.MonkeyPatch) -> None:
    # Regression: reptran emits at representative wavelengths that can start just
    # ABOVE seg_lo (and end below seg_hi). Stitching must clamp at the edges, not
    # leave NaN holes (which would poison any band whose RSRF reaches a seam - and
    # the adaptive default segments every scene).
    cfg = LibRadtranAlgorithmConfig(
        mol_abs_param="reptran",
        wavelength_min_nm=400.0,
        wavelength_max_nm=700.0,
        scene_lut_max_nodes_per_axis=2,
        scene_lut_max_cases=8,
        mol_abs_regions=[(500.0, 600.0, "reptran fine")],
        adaptive_deep_water_fine=False,
    )
    runner = LibRadtranRunner(libradtran_config=cfg)
    monkeypatch.setattr(
        runner,
        "_ensure_paths",
        lambda: LibRadtranPaths(uvspec=Path("/x/uvspec"), data_dir=Path("/x/data")),
    )

    def _fake_uvspec(_paths, deck: str):
        wl_line = next(line for line in deck.splitlines() if line.startswith("wavelength "))
        lo, hi = (float(x) for x in wl_line.split()[1:3])
        wl = np.arange(lo + 2.0, hi - 1.0, 1.0)  # NOT flush with the segment edges
        return wl, np.full(wl.shape, 0.8), np.full(wl.shape, 0.3)

    monkeypatch.setattr(runner, "_run_uvspec", _fake_uvspec)
    ds = runner._build_scene_lut(*_scene())
    for var in ("TOA_rho1", "TOA_rho2", "Eg_rho1", "Eg_rho2"):
        assert np.isfinite(ds[var].values).all(), f"{var} has NaN holes at segment edges"


def test_wavelength_grid_no_overshoot_for_non_integer_step(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # np.arange(stop=max+step) could place a sample past wavelength_max_nm for a
    # non-integer step; the endpoint-safe grid must pin both ends exactly.
    cfg = LibRadtranAlgorithmConfig(
        mol_abs_param="reptran",
        wavelength_min_nm=400.0,
        wavelength_max_nm=401.0,
        wavelength_step_nm=0.1,
        scene_lut_max_nodes_per_axis=2,
        scene_lut_max_cases=8,
        adaptive_deep_water_fine=False,
    )
    runner = LibRadtranRunner(libradtran_config=cfg)
    monkeypatch.setattr(
        runner,
        "_ensure_paths",
        lambda: LibRadtranPaths(uvspec=Path("/x/uvspec"), data_dir=Path("/x/data")),
    )

    def _fake(_paths, deck: str):
        wl_line = next(line for line in deck.splitlines() if line.startswith("wavelength "))
        lo, hi = (float(x) for x in wl_line.split()[1:3])
        wl = np.linspace(lo, hi, 200)
        return wl, np.full(wl.shape, 0.8), np.full(wl.shape, 0.3)

    monkeypatch.setattr(runner, "_run_uvspec", _fake)
    wl = runner._build_scene_lut(*_scene())["wavelength"].values
    assert wl.min() == pytest.approx(400.0)
    assert wl.max() == pytest.approx(401.0)  # exactly the window max, no overshoot
    assert np.all(wl <= 401.0 + 1e-9)


def test_mol_abs_regions_validation() -> None:
    ok = LibRadtranAlgorithmConfig(
        wavelength_min_nm=340.0,
        wavelength_max_nm=2300.0,
        mol_abs_regions=[(1085.0, 1180.0, "reptran fine"), (890.0, 985.0, "reptran fine")],
    )
    # validator sorts by lo.
    assert [r[0] for r in ok.mol_abs_regions] == [890.0, 1085.0]
    with pytest.raises(ValueError, match="overlap"):
        LibRadtranAlgorithmConfig(
            wavelength_min_nm=340.0,
            wavelength_max_nm=2300.0,
            mol_abs_regions=[(900.0, 1000.0, "reptran fine"), (950.0, 1100.0, "reptran medium")],
        )
    with pytest.raises(ValueError, match="within"):
        LibRadtranAlgorithmConfig(
            wavelength_min_nm=400.0,
            wavelength_max_nm=2300.0,
            mol_abs_regions=[(300.0, 500.0, "reptran fine")],
        )
    with pytest.raises(ValueError, match="lo_nm < hi_nm"):
        LibRadtranAlgorithmConfig(
            wavelength_min_nm=400.0,
            wavelength_max_nm=2300.0,
            mol_abs_regions=[(900.0, 900.0, "reptran fine")],
        )


def test_warn_band_straddle(caplog: pytest.LogCaptureFixture) -> None:
    cfg = LibRadtranAlgorithmConfig(
        mol_abs_param="reptran",
        wavelength_min_nm=400.0,
        wavelength_max_nm=2300.0,
        mol_abs_regions=[(890.0, 985.0, "reptran fine")],
    )
    runner = LibRadtranRunner(libradtran_config=cfg)
    # A band whose support [855, 905] crosses the 890 nm boundary -> warn.
    straddling = SensorBand(
        name="BX", center_wavelength=880.0, bandwidth=25.0, resolution=20.0, band_index=0
    )
    with caplog.at_level("WARNING"):
        runner._warn_band_straddle(straddling)
    assert any("straddles" in r.message for r in caplog.records)
    # A band fully inside a window (B04 ~665) does not warn.
    caplog.clear()
    clear = SensorBand(
        name="B04", center_wavelength=665.0, bandwidth=30.0, resolution=10.0, band_index=3
    )
    with caplog.at_level("WARNING"):
        runner._warn_band_straddle(clear)
    assert not any("straddles" in r.message for r in caplog.records)


def test_deep_water_bands_constant_avoids_s2_surface_centers() -> None:
    # The recommended fine regions must not contain any S2 surface band centre
    # (else a per-region scheme would split that band across a seam).
    s2_centers = [443, 490, 560, 665, 705, 740, 783, 833, 865, 1614, 2202]
    for lo, hi in DEEP_WATER_H2O_BANDS_NM:
        assert not any(lo < c < hi for c in s2_centers), (lo, hi)


def _fake_spectrum() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    wl = np.array([490.0, 491.0, 492.0], dtype=np.float64)
    eg = np.array([0.80, 0.81, 0.82], dtype=np.float64)
    toa = np.array([0.10, 0.11, 0.12], dtype=np.float64)
    return wl, eg, toa


def test_run_cache_skips_uvspec_on_repeat_deck(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An identical deck is served from disk on the second call (and re-runs)."""
    runner = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig(run_cache_dir=tmp_path))
    assert runner._run_cache_dir == tmp_path
    paths = LibRadtranPaths(uvspec=tmp_path / "uvspec", data_dir=tmp_path / "data")
    calls = {"n": 0}

    def _fake_invoke(
        _paths: LibRadtranPaths, _deck: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        calls["n"] += 1
        return _fake_spectrum()

    monkeypatch.setattr(runner, "_invoke_uvspec", _fake_invoke)

    deck = "mol_abs_param reptran\nwavelength 490 492\n"
    wl1, eg1, toa1 = runner._run_uvspec(paths, deck)
    wl2, eg2, toa2 = runner._run_uvspec(paths, deck)

    assert calls["n"] == 1  # second call served from disk, not uvspec
    assert (runner._run_cache_hits, runner._run_cache_misses) == (1, 1)
    np.testing.assert_array_equal(wl1, wl2)
    np.testing.assert_array_equal(eg1, eg2)
    np.testing.assert_array_equal(toa1, toa2)

    # A fresh runner sharing the cache dir also hits (cross-run / cross-process).
    runner2 = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig(run_cache_dir=tmp_path))

    def _must_not_run(*_a: object, **_k: object) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        pytest.fail("uvspec must not run for a cached deck")

    monkeypatch.setattr(runner2, "_invoke_uvspec", _must_not_run)
    wl3, _eg3, _toa3 = runner2._run_uvspec(paths, deck)
    np.testing.assert_array_equal(wl1, wl3)
    assert (runner2._run_cache_hits, runner2._run_cache_misses) == (1, 0)


def test_run_cache_disabled_always_invokes_uvspec(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    runner = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig(run_cache_enabled=False))
    assert runner._run_cache_dir is None
    paths = LibRadtranPaths(uvspec=tmp_path / "uvspec", data_dir=tmp_path / "data")
    calls = {"n": 0}

    def _fake_invoke(
        _paths: LibRadtranPaths, _deck: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        calls["n"] += 1
        return _fake_spectrum()

    monkeypatch.setattr(runner, "_invoke_uvspec", _fake_invoke)
    runner._run_uvspec(paths, "same-deck")
    runner._run_uvspec(paths, "same-deck")
    assert calls["n"] == 2  # caching off => uvspec runs every time


def test_run_cache_distinct_decks_recompute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    runner = LibRadtranRunner(libradtran_config=LibRadtranAlgorithmConfig(run_cache_dir=tmp_path))
    paths = LibRadtranPaths(uvspec=tmp_path / "uvspec", data_dir=tmp_path / "data")
    calls = {"n": 0}

    def _fake_invoke(
        _paths: LibRadtranPaths, _deck: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        calls["n"] += 1
        return _fake_spectrum()

    monkeypatch.setattr(runner, "_invoke_uvspec", _fake_invoke)
    runner._run_uvspec(paths, "deck-A")
    runner._run_uvspec(paths, "deck-B")
    assert calls["n"] == 2  # distinct decks => distinct keys => both computed
    assert runner._run_cache_misses == 2
