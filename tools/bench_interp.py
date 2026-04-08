#!/usr/bin/env python
"""Benchmark: xr.DataArray.interp() vs scipy RegularGridInterpolator.

Compares wall-clock time and numerical agreement for:
  1. 2D (aot, tcwv) grids — the spectral LUT path after scene subsetting
  2. 5D (sza, vza, raa, aot, tcwv) grids — the compact LUT path

Both paths are used in SIAC's LUT backend.
"""

import time

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator


# ---------------------------------------------------------------------------
# Grid builders
# ---------------------------------------------------------------------------

def make_2d_grid(n_aot: int = 16, n_tcwv: int = 12) -> xr.DataArray:
    aot = np.linspace(0.0, 2.0, n_aot, dtype=np.float64)
    tcwv = np.linspace(0.0, 7.0, n_tcwv, dtype=np.float64)
    rng = np.random.default_rng(42)
    return xr.DataArray(
        rng.random((n_aot, n_tcwv)),
        dims=("aot", "tcwv"),
        coords={"aot": aot, "tcwv": tcwv},
    )


def make_5d_grid(
    n_sza: int = 13,
    n_vza: int = 7,
    n_raa: int = 10,
    n_aot: int = 16,
    n_tcwv: int = 12,
) -> xr.DataArray:
    """Compact-LUT shape: (sza, vza, raa, aot, tcwv)."""
    rng = np.random.default_rng(42)
    return xr.DataArray(
        rng.random((n_sza, n_vza, n_raa, n_aot, n_tcwv)),
        dims=("sza", "vza", "raa", "aot", "tcwv"),
        coords={
            "sza": np.linspace(0, 75, n_sza),
            "vza": np.linspace(0, 45, n_vza),
            "raa": np.linspace(0, 180, n_raa),
            "aot": np.linspace(0, 2, n_aot),
            "tcwv": np.linspace(0, 7, n_tcwv),
        },
    )


# ---------------------------------------------------------------------------
# Query-point builders
# ---------------------------------------------------------------------------

def make_query_points(n_points: int, grid: xr.DataArray) -> dict[str, xr.DataArray]:
    rng = np.random.default_rng(7)
    coords = {}
    for name in grid.dims:
        lo = float(grid.coords[name].min())
        hi = float(grid.coords[name].max())
        coords[name] = xr.DataArray(rng.uniform(lo, hi, n_points), dims=["point"])
    return coords


# ---------------------------------------------------------------------------
# Interpolation methods
# ---------------------------------------------------------------------------

def interp_xarray(grid: xr.DataArray, coords: dict[str, xr.DataArray]) -> np.ndarray:
    return grid.interp(**coords).values.astype(np.float32)


def interp_scipy(grid: xr.DataArray, coords: dict[str, xr.DataArray]) -> np.ndarray:
    dim_names = list(grid.dims)
    grid_axes = tuple(grid.coords[d].values.astype(np.float64) for d in dim_names)
    interp = RegularGridInterpolator(
        grid_axes,
        grid.values.astype(np.float64),
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    pts = np.column_stack([coords[d].values.astype(np.float64) for d in dim_names])
    return interp(pts).astype(np.float32)


# ---------------------------------------------------------------------------
# Timing harness
# ---------------------------------------------------------------------------

def bench(label: str, fn, grid, coords, n_repeats: int = 3):
    fn(grid, coords)  # warm-up
    times = []
    for _ in range(n_repeats):
        t0 = time.perf_counter()
        result = fn(grid, coords)
        times.append(time.perf_counter() - t0)
    best = min(times)
    print(f"  {label:25s}  best={best*1000:9.2f} ms  (n_repeats={n_repeats})")
    return result, best


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_suite(label: str, grid: xr.DataArray, point_counts: list[int]):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  Grid shape: {grid.shape}  dims={grid.dims}")
    print(f"{'='*60}")
    for n_points in point_counts:
        print(f"\n--- {n_points:,} query points ---")
        coords = make_query_points(n_points, grid)

        result_xa, t_xa = bench("xr.DataArray.interp()", interp_xarray, grid, coords)
        result_sp, t_sp = bench("RegularGridInterpolator", interp_scipy, grid, coords)

        max_diff = float(np.nanmax(np.abs(result_xa - result_sp)))
        speedup = t_xa / t_sp if t_sp > 0 else float("inf")
        print(f"  max |diff| = {max_diff:.2e},  speedup = {speedup:.1f}x")


def main():
    point_counts = [100, 1_000, 10_000, 48_400]

    # 2D — spectral LUT path (scene-subsetted grids)
    run_suite("2D spectral-LUT path  (aot × tcwv)", make_2d_grid(), point_counts)

    # 5D — compact LUT path
    run_suite("5D compact-LUT path  (sza × vza × raa × aot × tcwv)", make_5d_grid(), point_counts)

    # Simulate the full M3 query: 12 calls × 2D grid  (3 bands × 4 variables)
    print(f"\n{'='*60}")
    print("  Simulated M3 surface prior:  12 calls × 48,400 pts on 2D grid")
    print(f"{'='*60}")
    grid = make_2d_grid()
    coords = make_query_points(48_400, grid)

    t_xa_total = 0.0
    t_sp_total = 0.0
    for i in range(12):
        t0 = time.perf_counter()
        interp_xarray(grid, coords)
        t_xa_total += time.perf_counter() - t0

        t0 = time.perf_counter()
        interp_scipy(grid, coords)
        t_sp_total += time.perf_counter() - t0

    print(f"  xr.DataArray.interp()    total = {t_xa_total*1000:9.2f} ms")
    print(f"  RegularGridInterpolator  total = {t_sp_total*1000:9.2f} ms")
    print(f"  speedup = {t_xa_total / t_sp_total:.1f}x")


if __name__ == "__main__":
    main()

