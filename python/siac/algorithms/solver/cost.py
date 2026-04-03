"""
Cost function components for aerosol retrieval.

This module implements the cost function and its gradient for the
multi-grid optimization of atmospheric parameters (AOT, TCWV).

The total cost function is:
    J = J_obs + J_prior + J_smooth

where:
    J_obs = Σ_bands Σ_pixels w_band * (boa_model - boa_prior)² / σ²_boa
    J_prior = (aot - aot_prior)²/σ²_aot + (tcwv - tcwv_prior)²/σ²_tcwv
    J_smooth = γ_aot * Σ φ_δ(∇aot) + γ_tcwv * Σ φ_δ(∇tcwv)

The smoothness term uses a Pseudo-Huber penalty on first-order finite
differences: φ_δ(a) = δ²(√(1 + (a/δ)²) − 1).  This behaves quadratically
for small gradients (noise suppression) and linearly for large gradients
(hotspot preservation).  The parameter δ controls the transition threshold.

The gradient dJ/d(aot, tcwv) is computed analytically using the chain rule
through the RT model Jacobians.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from scipy import sparse

from siac.domain.protocols import RTModelBackend

if TYPE_CHECKING:
    from siac.domain import SensorBand
    from siac.runtime import AtmosphericState, GeometryAngles, SurfacePrior

logger = logging.getLogger(__name__)


@dataclass
class CostFunctionConfig:
    """Configuration for the cost function."""

    # Smoothness regularization weights
    aot_gamma: float = 10.0
    tcwv_gamma: float = 5.0

    # Parameter bounds
    aot_min: float = 0.001
    aot_max: float = 2.5
    tcwv_min: float = 0.0
    tcwv_max: float = 8.0

    # Band weighting power (negative values weight shorter wavelengths more for AOT)
    band_weight_power: float = -1.6

    # Pseudo-Huber transition threshold for smoothness.
    # Gradients below delta are penalized quadratically (noise suppression);
    # gradients above delta are penalized linearly (hotspot preservation).
    # Units match the field: AOT per grid cell for aot, g/cm² per cell for tcwv.
    smoothness_delta: float = 0.02

    # Minimum uncertainty floor
    min_boa_unc: float = 0.001
    min_aot_unc: float = 0.02
    min_tcwv_unc: float = 0.1


class CostFunction:
    """
    Cost function for aerosol retrieval optimization.

    Implements the observation, prior, and smoothness cost terms
    with analytical gradients for use with L-BFGS-B optimization.
    """

    def __init__(
        self,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: RTModelBackend,
        bands: list[SensorBand],
        mask: xr.DataArray,
        config: CostFunctionConfig | None = None,
    ):
        """
        Initialize cost function.

        Args:
            toa: TOA reflectance (bands, y, x)
            surface_prior: Surface reflectance prior from BRDF
            geometry: Viewing geometry
            atmo_prior: Prior atmospheric state
            rt_model: Radiative transfer model backend
            bands: List of sensor bands to use
            mask: Valid pixel mask (True = valid)
            config: Cost function configuration
        """
        if not isinstance(rt_model, RTModelBackend):
            raise TypeError(
                f"rt_model must implement RTModelBackend protocol, "
                f"got {type(rt_model).__name__}"
            )
        self.toa = toa
        self.surface_prior = surface_prior
        self.geometry = geometry
        self.atmo_prior = atmo_prior
        self.rt_model = rt_model
        self.bands = bands
        self.mask = mask.values
        self.config = config or CostFunctionConfig()

        # Get dimensions
        self.shape = (mask.sizes.get("y", mask.shape[0]), mask.sizes.get("x", mask.shape[1]))
        self.n_bands = len(bands)
        self.n_pixels = np.prod(self.shape)

        # Setup prior arrays
        self._setup_priors()

        # Pre-allocate mutable DataArrays for the optimizer inner loop to
        # avoid reconstructing xarray objects on every L-BFGS-B iteration.
        self._aot_da = xr.DataArray(
            self.atmo_prior.aot.values.copy(),
            dims=self.atmo_prior.aot.dims,
            coords=self.atmo_prior.aot.coords,
        )
        self._tcwv_da = xr.DataArray(
            self.atmo_prior.tcwv.values.copy(),
            dims=self.atmo_prior.tcwv.dims,
            coords=self.atmo_prior.tcwv.coords,
        )

        # Setup band weights
        self._setup_band_weights()

        # Setup smoothness operator
        self._setup_smoothness()

    def _setup_priors(self) -> None:
        """Setup prior arrays and uncertainties."""
        # Prior values
        self.aot_prior = self.atmo_prior.aot.values.copy()
        self.tcwv_prior = self.atmo_prior.tcwv.values.copy()

        # Prior uncertainties (with minimum floor)
        self.aot_unc = np.maximum(
            self.atmo_prior.aot_unc.values, self.config.min_aot_unc
        )
        self.tcwv_unc = np.maximum(
            self.atmo_prior.tcwv_unc.values, self.config.min_tcwv_unc
        )

        # BOA prior and uncertainty
        self.boa_prior = self.surface_prior.boa.values
        self.boa_unc = np.maximum(
            self.surface_prior.boa_unc.values, self.config.min_boa_unc
        )

        # Ensure finite values – fall back to the configured minimum uncertainty
        # rather than an arbitrary 1.0 which is far too large for AOT/TCWV.
        self.aot_unc = np.where(
            np.isfinite(self.aot_unc), self.aot_unc, self.config.min_aot_unc
        )
        self.tcwv_unc = np.where(
            np.isfinite(self.tcwv_unc), self.tcwv_unc, self.config.min_tcwv_unc
        )
        self.boa_unc = np.where(
            np.isfinite(self.boa_unc), self.boa_unc, self.config.min_boa_unc
        )

    def _setup_band_weights(self) -> None:
        """Setup wavelength-based band weights."""
        wavelengths = np.array([b.center_wavelength for b in self.bands])
        # Convert to micrometers for weighting
        wl_um = wavelengths / 1000.0

        # Power-law weighting (shorter wavelengths weighted more for AOT)
        weights = wl_um ** self.config.band_weight_power
        self.band_weights = weights / weights.sum()

    def _setup_smoothness(self) -> None:
        """Setup smoothness regularization (no precomputation needed)."""

    def __call__(self, p: np.ndarray) -> tuple[float, np.ndarray]:
        """
        Compute cost function and gradient.

        Args:
            p: Parameter vector [aot.ravel(), tcwv.ravel()]

        Returns:
            Tuple of (cost, gradient)
        """
        return self.cost_and_gradient(p)

    def cost_and_gradient(self, p: np.ndarray) -> tuple[float, np.ndarray]:
        """
        Compute total cost function and gradient.

        Args:
            p: Parameter vector [aot.ravel(), tcwv.ravel()]

        Returns:
            Tuple of (cost, gradient)
        """
        # Unpack parameters
        aot, tcwv = self._unpack_params(p)

        # Compute observation cost
        j_obs, dj_obs_aot, dj_obs_tcwv = self._observation_cost(aot, tcwv)

        # Compute prior cost
        j_prior, dj_prior_aot, dj_prior_tcwv = self._prior_cost(aot, tcwv)

        # Compute smoothness cost
        j_smooth, dj_smooth_aot, dj_smooth_tcwv = self._smoothness_cost(aot, tcwv)

        # Total cost
        j_total = j_obs + j_prior + j_smooth

        # Total gradient
        dj_aot = dj_obs_aot + dj_prior_aot + dj_smooth_aot
        dj_tcwv = dj_obs_tcwv + dj_prior_tcwv + dj_smooth_tcwv

        gradient = self._pack_params(dj_aot, dj_tcwv)

        return j_total, gradient

    def _unpack_params(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Unpack parameter vector into AOT and TCWV arrays."""
        n = self.n_pixels
        aot = p[:n].reshape(self.shape)
        tcwv = p[n:2*n].reshape(self.shape)

        # Apply bounds
        aot = np.clip(aot, self.config.aot_min, self.config.aot_max)
        tcwv = np.clip(tcwv, self.config.tcwv_min, self.config.tcwv_max)

        return aot, tcwv

    def _pack_params(self, aot: np.ndarray, tcwv: np.ndarray) -> np.ndarray:
        """Pack AOT and TCWV arrays into parameter vector."""
        return np.concatenate([aot.ravel(), tcwv.ravel()])

    def _observation_cost(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Compute observation cost term.

        J_obs = Σ_bands Σ_pixels w_band * (boa_model - boa_prior)² / σ²_boa

        Returns:
            Tuple of (cost, d_cost/d_aot, d_cost/d_tcwv)
        """
        j_obs = 0.0
        dj_aot = np.zeros(self.shape)
        dj_tcwv = np.zeros(self.shape)

        # Update pre-allocated DataArrays in-place to avoid xarray overhead
        self._aot_da.values[:] = aot
        self._tcwv_da.values[:] = tcwv
        atmo_state = self.atmo_prior.with_updated_aot_tcwv(
            aot=self._aot_da,
            tcwv=self._tcwv_da,
        )

        # Compute RT coefficients for all bands (batch if supported)
        if hasattr(self.rt_model, 'compute_coefficients_multi'):
            all_coeffs = self.rt_model.compute_coefficients_multi(
                self.geometry, atmo_state, self.bands, compute_jacobian=True,
            )
        else:
            all_coeffs = [
                self.rt_model.compute_coefficients(
                    self.geometry, atmo_state, band, compute_jacobian=True,
                )
                for band in self.bands
            ]

        for i, (_band, coeffs) in enumerate(zip(self.bands, all_coeffs)):

            # Apply correction to get modeled BOA
            toa_band_da = self.toa.isel(band=i) if "band" in self.toa.dims else self.toa
            boa_model = coeffs.apply_correction(toa_band_da).values

            # Get prior BOA for this band
            if self.boa_prior.ndim == 3:
                boa_prior_band = self.boa_prior[i]
                boa_unc_band = self.boa_unc[i]
            else:
                boa_prior_band = self.boa_prior
                boa_unc_band = self.boa_unc
                if i == 0:
                    logger.warning(
                        "boa_prior is %dD; using the same surface prior for all %d bands.",
                        self.boa_prior.ndim,
                        len(self.bands),
                    )

            # Residual
            diff = boa_model - boa_prior_band
            weight = self.band_weights[i] / (boa_unc_band ** 2)

            # Mask invalid pixels
            valid = self.mask & np.isfinite(diff) & np.isfinite(weight)
            weight = np.where(valid, weight, 0.0)

            # Cost contribution
            j_band = 0.5 * np.sum(weight * diff ** 2)
            j_obs += j_band

            # Gradient via chain rule
            if coeffs.has_jacobian:
                d_boa_aot, d_boa_tcwv = coeffs.compute_boa_jacobian(
                    toa_band_da
                )

                dj_aot += weight * diff * d_boa_aot.values
                dj_tcwv += weight * diff * d_boa_tcwv.values

        return j_obs, dj_aot, dj_tcwv

    def _prior_cost(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Compute prior constraint cost term.

        J_prior = (aot - aot_prior)²/σ²_aot + (tcwv - tcwv_prior)²/σ²_tcwv

        Returns:
            Tuple of (cost, d_cost/d_aot, d_cost/d_tcwv)
        """
        # AOT prior cost
        diff_aot = aot - self.aot_prior
        weight_aot = 1.0 / (self.aot_unc ** 2)
        j_aot = 0.5 * np.sum(weight_aot * diff_aot ** 2)
        dj_aot = weight_aot * diff_aot

        # TCWV prior cost
        diff_tcwv = tcwv - self.tcwv_prior
        weight_tcwv = 1.0 / (self.tcwv_unc ** 2)
        j_tcwv = 0.5 * np.sum(weight_tcwv * diff_tcwv ** 2)
        dj_tcwv = weight_tcwv * diff_tcwv

        j_prior = j_aot + j_tcwv

        return j_prior, dj_aot, dj_tcwv

    # ------------------------------------------------------------------
    # Pseudo-Huber smoothness on finite differences
    # ------------------------------------------------------------------

    @staticmethod
    def _pseudo_huber_cost_grad(
        field: np.ndarray,
        gamma: float,
        delta: float,
    ) -> tuple[float, np.ndarray]:
        """Compute Pseudo-Huber smoothness cost and gradient for a 2-D field.

        The cost is  γ² Σ φ_δ(dx) + φ_δ(dy)
        where φ_δ(a) = δ²(√(1 + (a/δ)²) − 1)  (Pseudo-Huber loss).

        For |a| << δ this is ≈ a²/2  (quadratic, noise suppression).
        For |a| >> δ this is ≈ δ|a|  (linear, hotspot preservation).

        The gradient is computed via the adjoint (negative divergence) of the
        element-wise derivative  φ'(a) = a / √(1 + (a/δ)²).
        """
        # Forward finite differences with Neumann (zero-gradient) boundaries.
        dy = np.diff(field, axis=0)  # shape (ny-1, nx)
        dx = np.diff(field, axis=1)  # shape (ny, nx-1)

        delta2 = delta * delta

        # Pseudo-Huber: φ(a) = δ²(√(1 + (a/δ)²) − 1)
        r_dy = np.sqrt(1.0 + (dy * dy) / delta2)
        r_dx = np.sqrt(1.0 + (dx * dx) / delta2)

        cost = gamma * gamma * delta2 * (np.sum(r_dy - 1.0) + np.sum(r_dx - 1.0))

        # φ'(a) = a / √(1 + (a/δ)²)
        dphi_dy = dy / r_dy  # shape (ny-1, nx)
        dphi_dx = dx / r_dx  # shape (ny, nx-1)

        # Adjoint of forward-difference → gradient w.r.t. field.
        # ∂/∂f[i] Σ φ(f[i+1]-f[i]) = -φ'(f[i+1]-f[i]) + φ'(f[i]-f[i-1])
        grad = np.zeros_like(field)
        grad[:-1, :] -= dphi_dy   # -φ'(dy[i]) at source i
        grad[1:, :] += dphi_dy    # +φ'(dy[i]) at target i+1
        grad[:, :-1] -= dphi_dx   # -φ'(dx[j]) at source j
        grad[:, 1:] += dphi_dx    # +φ'(dx[j]) at target j+1

        grad *= gamma * gamma

        return cost, grad

    def _smoothness_cost(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Compute smoothness regularization cost term using Pseudo-Huber penalty.

        Uses first-order finite differences with a Pseudo-Huber loss that
        transitions from quadratic (smooth regions) to linear (sharp features).

        Returns:
            Tuple of (cost, d_cost/d_aot, d_cost/d_tcwv)
        """
        delta = self.config.smoothness_delta

        j_aot, dj_aot = self._pseudo_huber_cost_grad(
            aot, self.config.aot_gamma, delta,
        )
        j_tcwv, dj_tcwv = self._pseudo_huber_cost_grad(
            tcwv, self.config.tcwv_gamma, delta,
        )

        return j_aot + j_tcwv, dj_aot, dj_tcwv

    def observation_cost_only(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray]:
        """Compute only observation cost (for diagnostics)."""
        j_obs, dj_aot, dj_tcwv = self._observation_cost(aot, tcwv)
        return j_obs, self._pack_params(dj_aot, dj_tcwv)

    def prior_cost_only(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray]:
        """Compute only prior cost (for diagnostics)."""
        j_prior, dj_aot, dj_tcwv = self._prior_cost(aot, tcwv)
        return j_prior, self._pack_params(dj_aot, dj_tcwv)

    def smoothness_cost_only(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray]:
        """Compute only smoothness cost (for diagnostics)."""
        j_smooth, dj_aot, dj_tcwv = self._smoothness_cost(aot, tcwv)
        return j_smooth, self._pack_params(dj_aot, dj_tcwv)


def compute_laplacian_eigenvalues(nx: int, ny: int) -> np.ndarray:
    """
    Compute eigenvalues of discrete Laplacian for regularization.

    The eigenvalues in the DCT basis are:
        λ_{i,j} = 2 - 2*cos(πi/M) + 2 - 2*cos(πj/N)

    Args:
        nx: Grid width
        ny: Grid height

    Returns:
        2D array of eigenvalues
    """
    i = np.arange(ny)
    j = np.arange(nx)
    jj, ii = np.meshgrid(j, i)

    lambda_vals = 2 - 2 * np.cos(np.pi * ii / ny) + 2 - 2 * np.cos(np.pi * jj / nx)
    return lambda_vals


def apply_smoothness_filter(
    x: np.ndarray,
    gamma: float,
    _lambda_vals: np.ndarray | None = None,
    *,
    delta: float = 0.02,
    n_iter: int = 20,
) -> np.ndarray:
    """
    Apply edge-preserving smoothness filter via iterated Pseudo-Huber diffusion.

    Iteratively solves a weighted-Laplacian system where the weights are
    derived from the Pseudo-Huber penalty, approximating the proximal
    operator of the Pseudo-Huber smoothness term.

    For small gradients (|∇x| << delta) this behaves like Tikhonov smoothing.
    For large gradients (|∇x| >> delta) the smoothing is reduced, preserving
    sharp features such as aerosol hotspots.

    Args:
        x: Input 2-D array.
        gamma: Regularization strength.
        lambda_vals: Unused (kept for backward compatibility).
        delta: Pseudo-Huber transition threshold.
        n_iter: Number of diffusion iterations.

    Returns:
        Smoothed array.
    """
    z = x.astype(np.float64).copy()
    gamma2 = gamma * gamma
    delta2 = delta * delta

    # Stable step size: the data-fidelity Hessian contributes 1 per pixel,
    # and the smoothness Hessian contributes at most γ² * 4 (4-connected).
    # For stability we need τ * (1 + γ² * 4) < 1.
    tau = 1.0 / (1.0 + 4.0 * gamma2)

    for _ in range(n_iter):
        # Forward differences
        dy = np.diff(z, axis=0)
        dx = np.diff(z, axis=1)

        # Pseudo-Huber weight: w(a) = 1 / √(1 + (a/δ)²)
        # This down-weights large gradients (edges).
        wy = 1.0 / np.sqrt(1.0 + (dy * dy) / delta2)
        wx = 1.0 / np.sqrt(1.0 + (dx * dx) / delta2)

        # Smoothness gradient (adjoint of weighted forward-difference)
        sg = np.zeros_like(z)
        sg[:-1, :] -= wy * dy
        sg[1:, :] += wy * dy
        sg[:, :-1] -= wx * dx
        sg[:, 1:] += wx * dx

        # Gradient descent step: z ← z − τ * ((z − x) + γ² * sg)
        z -= tau * ((z - x) + gamma2 * sg)

    return z


def create_sparse_laplacian(nx: int, ny: int) -> sparse.csc_matrix:
    """
    Create sparse discrete Laplacian matrix.

    Uses Neumann boundary conditions (zero gradient at edges).

    Args:
        nx: Grid width
        ny: Grid height

    Returns:
        Sparse Laplacian matrix
    """
    n = nx * ny

    # Diagonal values
    main_diag = 4 * np.ones(n)

    # Off-diagonals for x neighbors
    off_diag_1 = -1 * np.ones(n)
    off_diag_1[ny - 1 :: ny] = 0  # No connection across row boundaries

    # Off-diagonals for y neighbors
    off_diag_ny = -1 * np.ones(n)

    # Boundary adjustments (Neumann: reduce diagonal at edges)
    indices = np.arange(n)
    rows = indices // ny
    cols = indices % ny
    main_diag[rows == 0] -= 1.0
    main_diag[rows == nx - 1] -= 1.0
    main_diag[cols == 0] -= 1.0
    main_diag[cols == ny - 1] -= 1.0

    # Build sparse matrix
    laplacian = sparse.diags(
        [main_diag, off_diag_1[:-1], off_diag_1[:-1], off_diag_ny[:-ny], off_diag_ny[:-ny]],
        [0, 1, -1, ny, -ny],
        format="csc",
    )

    return laplacian
