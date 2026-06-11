"""Block-grid-search machinery for the multigrid aerosol solver.

Implements the no-Jacobian solver path of
``siac.algorithms.solver.multigrid.MultiGridSolver`` as module-level steps
with explicit parameters: observation/prior stacking and validation,
candidate RT-coefficient provider construction (including the joint-LUT
negotiation with capable backends), block-grid aggregation, and dispatch of
the candidate cost-cube evaluation kernels.

The Rust evaluation kernels are passed in by the caller rather than
imported here so they keep resolving through the ``multigrid`` module
namespace (tests monkeypatch them there).
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt

from siac.domain.protocols import rt_optional_capability
from siac.runtime import AtmosphericState, GeometryAngles

if TYPE_CHECKING:
    from siac.domain import SensorBand
    from siac.domain.protocols import RTModelBackend
    from siac.runtime import SurfacePrior

logger = logging.getLogger(__name__)

BoolArray: TypeAlias = npt.NDArray[np.bool_]
Float32Array: TypeAlias = npt.NDArray[np.float32]
# Mirrors the alias of the same name in multigrid.py.
FixedAtmosphericParameter: TypeAlias = Literal["none", "aot", "tcwv"]
# Invoked once per (aot, tcwv) candidate; returns the (xap, xbp, xcp) stacks.
CandidateCoeffProvider: TypeAlias = Callable[
    [float, float], tuple[np.ndarray, np.ndarray, np.ndarray]
]


def fixed_axis_from_prior(
    prior: np.ndarray,
    bounds: tuple[float, float],
) -> Float32Array:
    """Collapse a per-pixel prior to a single-node axis for a fixed parameter."""
    finite = prior[np.isfinite(prior)]
    value = float(np.mean(finite)) if finite.size else 0.5 * (float(bounds[0]) + float(bounds[1]))
    return cast(
        "Float32Array",
        np.array([np.clip(value, float(bounds[0]), float(bounds[1]))], dtype=np.float32),
    )


def _finite_fixed_prior_field(
    prior: np.ndarray,
    bounds: tuple[float, float],
) -> Float32Array:
    values = np.asarray(prior, dtype=np.float32)
    if np.all(np.isfinite(values)):
        return values.astype(np.float32, copy=True)

    fixed_axis = fixed_axis_from_prior(values, bounds)
    filled = values.astype(np.float32, copy=True)
    filled[~np.isfinite(filled)] = fixed_axis[0]
    return cast("Float32Array", filled)


def _subsample_geometry(geometry: GeometryAngles, step: int) -> GeometryAngles:
    """Return a subsampled copy of the geometry angles."""
    return GeometryAngles(
        sza=geometry.sza[::step, ::step],
        saa=geometry.saa[::step, ::step],
        vza=geometry.vza[::step, ::step],
        vaa=geometry.vaa[::step, ::step],
    )


def _subsample_atmo_state(state: AtmosphericState, step: int) -> AtmosphericState:
    """Return a subsampled copy of the atmospheric state."""
    return AtmosphericState(
        aot=state.aot[::step, ::step],
        tcwv=state.tcwv[::step, ::step],
        tco3=state.tco3[::step, ::step],
        aot_unc=state.aot_unc[::step, ::step],
        tcwv_unc=state.tcwv_unc[::step, ::step],
        tco3_unc=state.tco3_unc[::step, ::step],
        elevation=state.elevation[::step, ::step],
    )


def _empty_float32_like(template: xr.DataArray) -> xr.DataArray:
    return xr.DataArray(
        np.empty(template.shape, dtype=np.float32),
        dims=template.dims,
        coords=template.coords,
        attrs=template.attrs,
        name=template.name,
    )


def _assign_coeff_stack(
    target: Float32Array,
    band_index: int,
    values: xr.DataArray,
) -> None:
    coeff_values = np.asarray(values.values, dtype=np.float32)
    expected_shape = target[band_index].shape
    if coeff_values.shape != expected_shape:
        raise ValueError(
            f"RT coefficients for band {band_index} have shape {coeff_values.shape}, "
            f"expected {expected_shape}"
        )
    target[band_index] = coeff_values


def prepare_grid_search_observations(
    *,
    toa: xr.DataArray,
    surface_prior: SurfacePrior,
    n_bands: int,
    shape: tuple[int, int],
    min_boa_unc: float,
    aot_prior: Float32Array,
    tcwv_prior: Float32Array,
    aot_prior_unc: Float32Array,
    tcwv_prior_unc: Float32Array,
) -> tuple[Float32Array, BoolArray, Float32Array, Float32Array, BoolArray]:
    """Stack and validate the TOA/BOA arrays used by the grid search.

    Returns ``(toa_values, no_observation_mask, boa_prior, boa_unc,
    support_mask)`` where ``support_mask`` flags pixels with full
    observation support and finite surface/atmospheric priors.
    """
    toa_values = toa.values.astype(np.float32)
    if toa_values.ndim == 2:
        toa_values = toa_values[np.newaxis, ...]
    if toa_values.shape != (n_bands, *shape):
        raise ValueError(
            f"TOA shape {toa_values.shape} incompatible with {n_bands} bands and grid {shape}"
        )
    toa_values = np.ascontiguousarray(toa_values, dtype=np.float32)
    observation_band_mask = np.isfinite(toa_values) & (toa_values > 0.0) & (toa_values < 1.0)
    observation_support_mask = np.all(observation_band_mask, axis=0)
    no_observation_mask = ~np.any(observation_band_mask, axis=0)

    boa_prior = surface_prior.boa.values.astype(np.float32)
    if boa_prior.ndim == 2:
        boa_prior = np.broadcast_to(boa_prior, (n_bands, *shape))
    if boa_prior.ndim != 3 or boa_prior.shape[-2:] != shape:
        raise ValueError(
            f"BOA prior shape {boa_prior.shape} incompatible with {n_bands} bands and grid {shape}"
        )
    if boa_prior.shape[0] < n_bands:
        raise ValueError(f"BOA prior has {boa_prior.shape[0]} bands, needs at least {n_bands}")
    if boa_prior.shape[0] > n_bands:
        boa_prior = boa_prior[:n_bands]
    boa_prior = np.ascontiguousarray(boa_prior, dtype=np.float32)

    boa_unc = np.maximum(surface_prior.boa_unc.values.astype(np.float32), min_boa_unc)
    if boa_unc.ndim == 2:
        boa_unc = np.broadcast_to(boa_unc, (n_bands, *shape))
    if boa_unc.ndim != 3 or boa_unc.shape[-2:] != shape:
        raise ValueError(
            f"BOA uncertainty shape {boa_unc.shape} incompatible with {n_bands} bands and grid {shape}"
        )
    if boa_unc.shape[0] < n_bands:
        raise ValueError(f"BOA uncertainty has {boa_unc.shape[0]} bands, needs at least {n_bands}")
    if boa_unc.shape[0] > n_bands:
        boa_unc = boa_unc[:n_bands]
    boa_unc = np.ascontiguousarray(boa_unc, dtype=np.float32)
    prior_support_mask = (
        np.all(np.isfinite(boa_prior), axis=0)
        & np.all(np.isfinite(boa_unc), axis=0)
        & np.isfinite(aot_prior)
        & np.isfinite(tcwv_prior)
        & np.isfinite(aot_prior_unc)
        & np.isfinite(tcwv_prior_unc)
    )
    support_mask = observation_support_mask & prior_support_mask
    return toa_values, no_observation_mask, boa_prior, boa_unc, support_mask


def build_candidate_coeff_provider(
    *,
    rt_model: RTModelBackend,
    bands: list[SensorBand],
    geometry: GeometryAngles,
    atmo_prior: AtmosphericState,
    aot_axis: Float32Array,
    tcwv_axis: Float32Array,
    solve_aot: bool,
    solve_tcwv: bool,
    aot_bounds: tuple[float, float],
    tcwv_bounds: tuple[float, float],
    rt_sample_step: int,
    shape: tuple[int, int],
) -> CandidateCoeffProvider:
    """Build the per-candidate RT coefficient provider for the grid search.

    The provider fills shared ``(n_bands, ny, nx)`` coefficient stacks for
    one (aot, tcwv) candidate pair, sampling RT coefficients on the block
    grid when ``rt_sample_step > 1``.
    """
    n_bands = len(bands)
    coeff_spatial_shape = (
        shape
        if rt_sample_step <= 1
        else (
            (shape[0] + rt_sample_step - 1) // rt_sample_step,
            (shape[1] + rt_sample_step - 1) // rt_sample_step,
        )
    )

    xap_stack: Float32Array = np.empty((n_bands, *coeff_spatial_shape), dtype=np.float32)
    xbp_stack: Float32Array = np.empty((n_bands, *coeff_spatial_shape), dtype=np.float32)
    xcp_stack: Float32Array = np.empty((n_bands, *coeff_spatial_shape), dtype=np.float32)

    coeff_geometry = geometry
    coeff_atmo_template = atmo_prior
    if rt_sample_step > 1:
        coeff_geometry = _subsample_geometry(geometry, rt_sample_step)
        coeff_atmo_template = _subsample_atmo_state(atmo_prior, rt_sample_step)

    candidate_aot = _empty_float32_like(coeff_atmo_template.aot)
    candidate_tcwv = _empty_float32_like(coeff_atmo_template.tcwv)
    candidate_aot_values = cast("Float32Array", candidate_aot.data)
    candidate_tcwv_values = cast("Float32Array", candidate_tcwv.data)
    if not solve_aot:
        candidate_aot_values[...] = _finite_fixed_prior_field(
            coeff_atmo_template.aot.values,
            aot_bounds,
        )
    if not solve_tcwv:
        candidate_tcwv_values[...] = _finite_fixed_prior_field(
            coeff_atmo_template.tcwv.values,
            tcwv_bounds,
        )
    candidate_atmo = AtmosphericState(
        aot=candidate_aot,
        tcwv=candidate_tcwv,
        tco3=coeff_atmo_template.tco3,
        aot_unc=coeff_atmo_template.aot_unc,
        tcwv_unc=coeff_atmo_template.tcwv_unc,
        tco3_unc=coeff_atmo_template.tco3_unc,
        elevation=coeff_atmo_template.elevation,
    )

    # Try to build a single joint (aot × tcwv × geometry) LUT spanning the
    # entire grid-search range. The block-grid-search invokes the coeff
    # provider hundreds of times (N_aot × N_tcwv pairs); each call would
    # otherwise run a fresh 6S batch per band. The joint LUT amortises
    # those 6S calls across all candidates by precomputing one big LUT and
    # serving each candidate by interpolation instead. Because the LUT
    # nodes in the (aot, tcwv) axes coincide with the grid-search points,
    # the lookups are exact at the candidate values and only the geometric
    # dimensions are linearly interpolated — the same approximation the
    # per-candidate scene-LUT path already makes.
    #
    # IMPORTANT: pass ``coeff_atmo_template`` (the prior, with deterministic
    # per-pixel aot/tcwv values) — NOT ``candidate_atmo``, whose aot/tcwv
    # arrays are allocated with ``np.empty`` and may contain uninitialised
    # memory at this point (the values only get filled inside the provider
    # closure once a real candidate is in hand). Feeding uninitialised
    # values to the builder leaks NaNs into the LUT's ``valid_mask``
    # check, producing a non-deterministic valid pixel set across runs.
    # The aot/tcwv per-pixel arrays passed in here are only used to
    # populate the LUT's case_arrays for completeness; the actual aot/tcwv
    # axes the LUT spans come from the explicit ``aot_axis`` / ``tcwv_axis``
    # arguments below.
    #
    # When the backend doesn't support this optimization (non-6S backend,
    # or sixs.mode == "direct"), build_joint_grid_search_lut returns None
    # and we fall through to the original per-candidate compute path.
    joint_lut = None
    joint_lut_builder = rt_optional_capability(rt_model, "build_joint_grid_search_lut")
    if joint_lut_builder is not None:
        try:
            joint_lut = joint_lut_builder(
                geometry=coeff_geometry,
                atmo_state=coeff_atmo_template,
                aot_axis=aot_axis.astype(np.float64, copy=False),
                tcwv_axis=tcwv_axis.astype(np.float64, copy=False),
                bands=bands,
            )
        except Exception:
            logger.exception(
                "Failed to build joint grid-search LUT; falling back to the "
                "per-candidate scene-LUT compute path."
            )
            joint_lut = None

    if joint_lut is not None:
        logger.info(
            "Block-grid-search using joint LUT: %d aot × %d tcwv × %d bands "
            "(LUT case_count=%d, scene pixels=%d).",
            int(aot_axis.size),
            int(tcwv_axis.size),
            len(bands),
            joint_lut.plan.lut_case_count,
            joint_lut.plan.direct_case_count,
        )

        def _candidate_coeff_provider(
            aot_val: float, tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            if solve_aot:
                candidate_aot_values.fill(np.float32(aot_val))
            if solve_tcwv:
                candidate_tcwv_values.fill(np.float32(tcwv_val))
            band_outputs = joint_lut.evaluate(float(aot_val), float(tcwv_val))
            for ib, outputs in enumerate(band_outputs):
                _assign_coeff_stack(xap_stack, ib, outputs["xap"])
                _assign_coeff_stack(xbp_stack, ib, outputs["xbp"])
                _assign_coeff_stack(xcp_stack, ib, outputs["xcp"])
            return xap_stack, xbp_stack, xcp_stack
    else:

        def _candidate_coeff_provider(
            aot_val: float, tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            if solve_aot:
                candidate_aot_values.fill(np.float32(aot_val))
            if solve_tcwv:
                candidate_tcwv_values.fill(np.float32(tcwv_val))

            for ib, band in enumerate(bands):
                coeffs = rt_model.compute_coefficients(
                    coeff_geometry,
                    candidate_atmo,
                    band,
                    compute_jacobian=False,
                )
                _assign_coeff_stack(xap_stack, ib, coeffs.xap)
                _assign_coeff_stack(xbp_stack, ib, coeffs.xbp)
                _assign_coeff_stack(xcp_stack, ib, coeffs.xcp)
            return xap_stack, xbp_stack, xcp_stack

    return _candidate_coeff_provider


def aggregate_valid_counts(mask: np.ndarray, step: int) -> np.ndarray:
    """Count valid pixels contributing to each NxN block."""
    source = np.asarray(mask, dtype=bool)
    if step <= 1:
        return source.astype(np.int32, copy=False)

    by = (source.shape[0] + step - 1) // step
    bx = (source.shape[1] + step - 1) // step
    pad_y = by * step - source.shape[0]
    pad_x = bx * step - source.shape[1]
    padded = np.pad(source, ((0, pad_y), (0, pad_x)), mode="constant", constant_values=False)
    return padded.reshape(by, step, bx, step).sum(axis=(1, 3), dtype=np.int32)


def aggregate_block_pixel_counts(shape: tuple[int, int], step: int) -> np.ndarray:
    """Return the number of full-resolution pixels in each NxN block."""
    if step <= 1:
        return np.ones(shape, dtype=np.int32)

    by = (shape[0] + step - 1) // step
    bx = (shape[1] + step - 1) // step
    rows = np.minimum(step, shape[0] - np.arange(by, dtype=np.int32) * step)
    cols = np.minimum(step, shape[1] - np.arange(bx, dtype=np.int32) * step)
    return (rows[:, np.newaxis] * cols[np.newaxis, :]).astype(np.int32, copy=False)


def evaluate_candidate_cost_cube(
    *,
    coeff_provider: CandidateCoeffProvider,
    aot_axis: Float32Array,
    tcwv_axis: Float32Array,
    toa_values: Float32Array,
    boa_prior: Float32Array,
    boa_unc: Float32Array,
    band_weights: Float32Array,
    solve_valid_mask: BoolArray,
    aot_prior: Float32Array,
    tcwv_prior: Float32Array,
    aot_prior_unc: Float32Array,
    tcwv_prior_unc: Float32Array,
    block_size: int,
    fixed_parameter: FixedAtmosphericParameter,
    block_support_mask: np.ndarray,
    block_cost_cube_fn: Callable[..., Any],
    pixel_cost_cube_fn: Callable[..., Any],
) -> tuple[Float32Array, npt.NDArray[np.uint16], BoolArray, BoolArray]:
    """Evaluate the (aot × tcwv) candidate cost cube via the supplied kernel.

    Dispatches to ``block_cost_cube_fn`` when ``block_size > 1`` (one shared
    cost per NxN block) and ``pixel_cost_cube_fn`` otherwise, then masks
    unsupported blocks out of the cube. Returns ``(costs, obs_counts,
    block_valid_mask, refine_valid_mask)``.
    """
    block_valid_mask = block_support_mask.copy()
    if block_size > 1:
        costs_raw, obs_counts_raw, block_valid_mask_raw = block_cost_cube_fn(
            coeff_provider,
            aot_axis.astype(np.float32, copy=False),
            tcwv_axis.astype(np.float32, copy=False),
            toa_values,
            boa_prior,
            boa_unc,
            band_weights,
            solve_valid_mask.astype(bool, copy=False),
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
            block_size,
            fixed_parameter,
        )
        block_valid_mask = np.asarray(block_valid_mask_raw, dtype=bool) & block_support_mask
        refine_valid_mask = block_valid_mask.astype(bool, copy=False)
    else:
        costs_raw, obs_counts_raw = pixel_cost_cube_fn(
            coeff_provider,
            aot_axis.astype(np.float32, copy=False),
            tcwv_axis.astype(np.float32, copy=False),
            toa_values,
            boa_prior,
            boa_unc,
            band_weights,
            solve_valid_mask.astype(bool, copy=False),
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
            fixed_parameter,
        )
        refine_valid_mask = solve_valid_mask.astype(bool, copy=False)
    costs = np.asarray(costs_raw, dtype=np.float32)
    obs_counts = np.asarray(obs_counts_raw, dtype=np.uint16)
    if block_size > 1:
        costs = np.where(block_valid_mask[np.newaxis, np.newaxis, :, :], costs, np.inf).astype(
            np.float32,
            copy=False,
        )
        obs_counts = np.where(
            block_valid_mask[np.newaxis, np.newaxis, :, :],
            obs_counts,
            np.uint16(0),
        ).astype(np.uint16, copy=False)
    return costs, obs_counts, block_valid_mask, refine_valid_mask
