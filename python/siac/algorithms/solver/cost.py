"""
Cost function components for aerosol retrieval.

This module implements the cost function and its gradient for the
multi-grid optimization of atmospheric parameters (AOT, TCWV).

The total cost function is:
    J = J_obs + J_prior + J_smooth

where:
    J_obs = Σ_bands Σ_pixels w_band * (boa_model - boa_prior)² / σ²_boa
    J_prior = (aot - aot_prior)²/σ²_aot + (tcwv - tcwv_prior)²/σ²_tcwv
    J_smooth = γ_aot * ||∇aot||² + γ_tcwv * ||∇tcwv||²

The gradient dJ/d(aot, tcwv) is computed analytically using the chain rule
through the RT model Jacobians.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr
from scipy import sparse
from scipy.fftpack import dct, idct

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

    # Band weighting power (alpha = -2 weights shorter wavelengths more)
    band_weight_power: float = -2.0

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
        rt_model: Any,  # RTModelBackend protocol
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
        """Setup smoothness regularization operators."""
        ny, nx = self.shape

        # DCT-based smoothness using eigenvalue decomposition
        # This is more efficient than sparse matrix representation for large grids
        self._setup_dct_smoothness(nx, ny)

    def _setup_dct_smoothness(self, nx: int, ny: int) -> None:
        """Setup DCT-based smoothness operator."""
        # Create frequency grids for DCT
        c, r = np.mgrid[:ny, :nx]

        # Laplacian eigenvalues in DCT domain
        omega_x = np.pi * (c / ny)
        omega_y = np.pi * (r / nx)

        # Eigenvalues of discrete Laplacian
        self.lambda_smooth = 2 - 2 * np.cos(omega_x) + 2 - 2 * np.cos(omega_y)

        # Precompute filter for smoothness cost
        # Smoothness cost = gamma^2 * sum(lambda * |dct(x)|^2)
        # Gradient = 2 * gamma^2 * idct(lambda * dct(x))

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

        # Create temporary atmospheric state with current parameters
        atmo_state = self.atmo_prior.with_updated_aot_tcwv(
            aot=xr.DataArray(aot, dims=self.atmo_prior.aot.dims),
            tcwv=xr.DataArray(tcwv, dims=self.atmo_prior.tcwv.dims),
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

    def _smoothness_cost(
        self, aot: np.ndarray, tcwv: np.ndarray
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Compute smoothness regularization cost term.

        J_smooth = γ_aot * ||∇aot||² + γ_tcwv * ||∇tcwv||²

        Uses DCT-based implementation for efficiency.

        Returns:
            Tuple of (cost, d_cost/d_aot, d_cost/d_tcwv)
        """
        gamma_aot = self.config.aot_gamma
        gamma_tcwv = self.config.tcwv_gamma

        # AOT smoothness using DCT
        aot_dct = dct(dct(aot, axis=0, norm="ortho"), axis=1, norm="ortho")
        j_aot = 0.5 * gamma_aot ** 2 * np.sum(self.lambda_smooth * aot_dct ** 2)

        # Gradient in DCT domain then transform back
        dj_aot_dct = gamma_aot ** 2 * self.lambda_smooth * aot_dct
        dj_aot = idct(idct(dj_aot_dct, axis=1, norm="ortho"), axis=0, norm="ortho")

        # TCWV smoothness using DCT
        tcwv_dct = dct(dct(tcwv, axis=0, norm="ortho"), axis=1, norm="ortho")
        j_tcwv = 0.5 * gamma_tcwv ** 2 * np.sum(self.lambda_smooth * tcwv_dct ** 2)

        dj_tcwv_dct = gamma_tcwv ** 2 * self.lambda_smooth * tcwv_dct
        dj_tcwv = idct(idct(dj_tcwv_dct, axis=1, norm="ortho"), axis=0, norm="ortho")

        j_smooth = j_aot + j_tcwv

        return j_smooth, dj_aot, dj_tcwv

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
    x: np.ndarray, gamma: float, lambda_vals: np.ndarray
) -> np.ndarray:
    """
    Apply Tikhonov smoothness regularization via DCT.

    Solves: (I + γ² Δ) x_smooth = x

    In DCT domain: X_smooth = X / (1 + γ² Λ)

    Args:
        x: Input array
        gamma: Regularization strength
        lambda_vals: Laplacian eigenvalues

    Returns:
        Smoothed array
    """
    # Forward DCT
    x_dct = dct(dct(x, axis=0, norm="ortho"), axis=1, norm="ortho")

    # Apply filter in DCT domain
    filter_vals = 1.0 / (1.0 + gamma ** 2 * lambda_vals)
    x_dct_smooth = x_dct * filter_vals

    # Inverse DCT
    x_smooth = idct(idct(x_dct_smooth, axis=1, norm="ortho"), axis=0, norm="ortho")

    return x_smooth


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
    for i in range(n):
        row = i // ny
        col = i % ny

        if row == 0:
            main_diag[i] -= 1
        if row == nx - 1:
            main_diag[i] -= 1
        if col == 0:
            main_diag[i] -= 1
        if col == ny - 1:
            main_diag[i] -= 1

    # Build sparse matrix
    laplacian = sparse.diags(
        [main_diag, off_diag_1[:-1], off_diag_1[:-1], off_diag_ny[:-ny], off_diag_ny[:-ny]],
        [0, 1, -1, ny, -ny],
        format="csc",
    )

    return laplacian
