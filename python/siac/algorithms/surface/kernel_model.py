"""
Surface prior derivation using BRDF kernel model.

Computes surface reflectance prior from MODIS MCD43 BRDF parameters
using the Ross-Thick Li-Sparse kernel model, with optional PSF
convolution for scale matching.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr
from scipy import ndimage

from siac.algorithms.brdf.kernels import BRDFKernels, compute_reflectance
from siac.algorithms.surface.spectral_mapping import map_multispectral_reflectance
from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.algorithms.surface.spectral_mapping import SpectralMappingConfig
    from siac.domain import SensorBand


logger = logging.getLogger(__name__)
FloatArray = np.ndarray[Any, np.dtype[np.floating[Any]]]

# --- Dark-target surface-prior uncertainty inflation -----------------------
# Dark surfaces (dense forest, water-adjacent land) drive aerosol *over*-
# retrieval: a near-black blue prior forces the solver to attribute the observed
# blue TOA to aerosol, since increasing AOT is the only way to darken the modelled
# BOA enough to match the prior. The BRDF/MODIS prior is least reliable and the
# blue-band inversion most ill-conditioned exactly where the surface is darkest,
# so the prior should be trusted *less* there. We add a brightness-dependent floor
# to the per-band reflectance uncertainty:
#
#     sigma_eff = sqrt(sigma**2 + (gain * max(0, threshold - boa))**2)
#
# applied per band per pixel keyed on that band's own darkness. Bright pixels
# (boa >= threshold) are untouched; the effect concentrates on the blue band,
# which is both the darkest and the most AOT-sensitive.
#
# Default OFF (opt-in). A 66-site AERONET A/B lowers mean absolute error
# (0.084 -> 0.078) and fixes the worst dark over-retrievers (NEON_Bartlett
# 1.38 -> 1.15 flips within EE; Santarem and Chachoengsao pulled toward truth),
# but is within-EE-count-neutral (51/66 -> 50/66): keyed on brightness alone it
# cannot separate "dark AND over-retrieving" from "dark but already correct", so
# it over-corrects a few borderline dark scenes (ATTO-Campina 0.24 -> 0.11, into
# under). A discrepancy-keyed variant -- inflate only where the surface implies
# far more aerosol than the CAMS prior -- would target the failure cleanly; until
# that exists this stays a toggle, not a default. Threshold/gain are insensitive
# in [0.04, 0.06] / [0.6, 1.2]; the dark-pixel population dominates.
_APPLY_DARK_TARGET_UNC = False
_DARK_TARGET_REFLECTANCE_THRESHOLD = 0.06
_DARK_TARGET_UNC_GAIN = 0.6


def inflate_dark_target_uncertainty(boa: xr.DataArray, boa_unc: xr.DataArray) -> xr.DataArray:
    """Inflate surface-prior uncertainty for dark (low-reflectance) pixels.

    Returns ``boa_unc`` unchanged when the toggle is off. See the module note for
    the rationale and model. ``boa`` and ``boa_unc`` must share dims/shape.
    """
    if not _APPLY_DARK_TARGET_UNC:
        return boa_unc
    boa_v = np.asarray(boa.values, dtype=np.float32)
    unc_v = np.asarray(boa_unc.values, dtype=np.float32)
    deficit = np.maximum(np.float32(0.0), np.float32(_DARK_TARGET_REFLECTANCE_THRESHOLD) - boa_v)
    inflated = np.sqrt(unc_v**2 + (np.float32(_DARK_TARGET_UNC_GAIN) * deficit) ** 2)
    return boa_unc.copy(data=inflated.astype(unc_v.dtype))


# --- Visible surface-prior de-bias ----------------------------------------
# The MODIS-BRDF kernel prior is systematically darker than the real surface in
# the visible solver bands. Measured against the S2 BOA corrected at the known
# AERONET AOD (the true surface in the solver's own RT frame) over 66 sites:
# prior(B02) = -0.0133 + 0.928*ref, prior(B04) = -0.0076 + 0.950*ref, with the
# B02/B04 biases 0.91-correlated -- one common MODIS->S2 offset, not band noise.
# A dark prior makes the solver raise AOT to darken the modelled BOA down to it,
# i.e. it contributes to the over-retrieval that dominates the failure tail. We
# invert the per-band affine fit, corrected = (boa + a)/b.
#
# Default OFF -- a rigorous negative result. The de-bias correctly removes the
# *reflectance* offset (prior MAE 0.025 -> 0.018) and, at strength 0.3, nulls the
# *population* AOT bias (+0.022 -> -0.005). But it does NOT improve per-site
# retrieval: a 66-site A/B is net -2 within-EE at every strength, because the
# prior error is a small correctable bias (~0.021) sitting on LARGER irreducible
# per-site scatter (std ~0.028). A uniform correction nulls the mean but cannot
# touch the scatter, so it just relocates errors (fixes low-AOD over-retrievers,
# breaks others). SWIR cannot supply better per-site surface info (R^2 <= 0.24),
# so no global prior correction helps -- the MODIS-BRDF prior is at the achievable
# quality ceiling for the available information. Kept as an opt-in toggle.
_APPLY_PRIOR_DEBIAS = False
# band -> (a, b) such that the full de-bias is corrected = (boa + a) / b
_PRIOR_DEBIAS: dict[str, tuple[float, float]] = {
    "B02": (0.0133, 0.928),
    "B04": (0.0076, 0.950),
}
# Fraction of the affine de-bias to apply. The full correction (1.0) removes the
# reflectance bias but over-corrects the retrieved AOT, because the
# reflectance->AOT sensitivity is large at high AOD (flat cost cube): a uniform
# brightening drops high-AOD retrievals far more than low-AOD ones. A 66-site A/B
# at strength 1.0 flipped the population AOT bias +0.022 -> -0.051 (over -> under);
# strength 0.3 nulls it (-0.005). 0.3 is retained as the bias-nulling value should
# the toggle ever be enabled, though it does not improve within-EE (see above).
_PRIOR_DEBIAS_STRENGTH = 0.3


def debias_visible_prior(boa: xr.DataArray) -> xr.DataArray:
    """Affine de-bias of the visible surface prior (see module note).

    Returns ``boa`` unchanged when the toggle is off or it has no band dim.
    Bands without a calibration entry pass through (a=0, b=1).
    """
    if not _APPLY_PRIOR_DEBIAS or "band" not in getattr(boa, "dims", ()):
        return boa
    bands = [str(x) for x in boa.coords["band"].values]
    a = xr.DataArray(
        [_PRIOR_DEBIAS.get(b, (0.0, 1.0))[0] for b in bands],
        dims=["band"],
        coords={"band": boa.coords["band"]},
    )
    b = xr.DataArray(
        [_PRIOR_DEBIAS.get(bd, (0.0, 1.0))[1] for bd in bands],
        dims=["band"],
        coords={"band": boa.coords["band"]},
    )
    corrected = (boa + a) / b
    return boa + np.float32(_PRIOR_DEBIAS_STRENGTH) * (corrected - boa)


class KernelModelDeriver:
    """
    Derive surface reflectance prior from BRDF kernel weights.

    Uses the Ross-Thick Li-Sparse BRDF model to compute expected
    surface reflectance at the observation geometry.

    Args:
        psf_sigma_x: PSF standard deviation in x direction (pixels).
            If *None*, computed automatically from ``source_resolution_m``
            and ``target_resolution_m``.
        psf_sigma_y: PSF standard deviation in y direction (pixels).
            If *None*, uses the same value as *psf_sigma_x*.
        apply_psf: Whether to apply PSF convolution for scale matching.
        source_resolution_m: Spatial resolution of the BRDF source (metres).
        target_resolution_m: Spatial resolution of the target sensor (metres).
    """

    def __init__(
        self,
        psf_sigma_x: float | None = None,
        psf_sigma_y: float | None = None,
        apply_psf: bool = True,
        source_resolution_m: float = 500.0,
        target_resolution_m: float = 10.0,
    ) -> None:
        if psf_sigma_x is not None:
            self.psf_sigma_x = psf_sigma_x
            self.psf_sigma_y = psf_sigma_y if psf_sigma_y is not None else psf_sigma_x
        else:
            # Derive sigma (in target pixels) from the resolution ratio.
            # The PSF FWHM ≈ source_resolution / target_resolution, so
            # σ = (source / target) / (2 √(2 ln 2)) ≈ ratio / 2.355
            ratio = source_resolution_m / target_resolution_m
            sigma = ratio / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            self.psf_sigma_x = sigma
            self.psf_sigma_y = psf_sigma_y if psf_sigma_y is not None else sigma
        self.apply_psf = apply_psf
        # ``psf_sigma_x/y`` count pixels *at ``target_resolution_m``* — the
        # resolution the PSF was calibrated against (10 m for S2). The surface
        # prior, however, is built on the aerosol-retrieval grid, which is
        # usually coarser (default 120 m). ``_scale_psf_sigmas`` converts the
        # pixel sigma to the actual grid so the *physical* PSF footprint stays
        # constant; without it a 29.75 px blur at 120 m would span ~3.6 km
        # (~12x too much) and smear neighbouring dark targets into the prior.
        self.source_resolution_m = source_resolution_m
        self.target_resolution_m = target_resolution_m
        self._kernels = BRDFKernels(hb=2.0, br=1.0)

    def compute_surface_prior(
        self,
        brdf_weights: BRDFKernelWeights,
        geometry: GeometryAngles,
        psf_params: tuple[float, float] | None = None,
        *,
        source_bands: Sequence[SensorBand] | None = None,
        target_bands: Sequence[SensorBand] | None = None,
        spectral_library: SpectralMappingConfig | None = None,
        spectral_k_neighbors: int = 5,
        grid_resolution_m: float | None = None,
    ) -> SurfacePrior:
        """
        Compute surface reflectance prior from BRDF parameters.

        Args:
            brdf_weights: BRDF kernel coefficients (f0, f1, f2)
            geometry: Observation geometry (sza, vza, raa)
            psf_params: Optional (sigma_x, sigma_y) override for PSF
            grid_resolution_m: Pixel size (metres) of the grid the prior is
                computed on. Used to rescale the PSF sigma — calibrated in
                pixels at ``target_resolution_m`` — so the physical PSF
                footprint is preserved. ``None`` leaves the sigma unscaled
                (legacy single-resolution behaviour).

        Returns:
            SurfacePrior with BOA reflectance, uncertainty, and mask
        """
        if psf_params is not None:
            sigma_x, sigma_y = psf_params
        else:
            sigma_x, sigma_y = self.psf_sigma_x, self.psf_sigma_y
        sigma_x, sigma_y = self._scale_psf_sigmas(sigma_x, sigma_y, grid_resolution_m)

        # Compute kernel values at observation geometry
        k_vol, k_geo = self._kernels.compute(geometry.vza, geometry.sza, geometry.raa)
        k_vol_da = cast("xr.DataArray", k_vol)
        k_geo_da = cast("xr.DataArray", k_geo)
        k_vol_da, k_geo_da = self._align_kernels_to_brdf_grid(
            k_vol_da,
            k_geo_da,
            brdf_weights.f0,
        )

        # Compute reflectance per band
        boa = cast(
            "xr.DataArray",
            compute_reflectance(
                brdf_weights.f0,
                brdf_weights.f1,
                brdf_weights.f2,
                k_vol_da,
                k_geo_da,
            ),
        )

        # Compute uncertainty
        boa_unc = self._compute_uncertainty(
            brdf_weights,
            k_vol_da,
            k_geo_da,
        )

        if source_bands and target_bands:
            boa, mapped_unc = map_multispectral_reflectance(
                boa,
                source_bands=source_bands,
                target_bands=target_bands,
                source_uncertainty=boa_unc,
                spectral_library=spectral_library,
                k_neighbors=spectral_k_neighbors,
            )
            boa_unc = mapped_unc

        # Apply PSF convolution if needed
        if self.apply_psf and sigma_x > 0 and sigma_y > 0:
            boa = self._apply_psf(boa, sigma_x, sigma_y)
            boa_unc = self._apply_psf(boa_unc, sigma_x, sigma_y)

        # De-bias the visible prior (MODIS->S2 dark offset) before downstream use.
        boa = debias_visible_prior(boa)
        # Loosen the prior over dark targets to curb aerosol over-retrieval.
        boa_unc = inflate_dark_target_uncertainty(boa, boa_unc)

        # Create validity mask. The kernel-model prior is built from BRDF
        # kernels and should produce non-negative reflectance; any negative
        # value is a numerical artefact rather than acceptable correction
        # noise, so we use ``> 0`` here (not ``> BOA_VALID_MIN``).
        from siac.constants import BOA_VALID_MAX

        mask = (boa > 0) & (boa < BOA_VALID_MAX) & np.isfinite(boa)

        return SurfacePrior(
            boa=boa,
            boa_unc=boa_unc,
            kernels=brdf_weights,
            mask=mask,
        )

    def _align_kernels_to_brdf_grid(
        self,
        k_vol: xr.DataArray,
        k_geo: xr.DataArray,
        ref: xr.DataArray,
    ) -> tuple[xr.DataArray, xr.DataArray]:
        """Align geometry kernels to the BRDF grid to avoid empty coord intersections."""
        if not isinstance(k_vol, xr.DataArray) or not isinstance(k_geo, xr.DataArray):
            return k_vol, k_geo
        if not isinstance(ref, xr.DataArray):
            return k_vol, k_geo

        if "y" not in ref.dims or "x" not in ref.dims:
            return k_vol, k_geo

        target_y = ref.coords.get("y")
        target_x = ref.coords.get("x")
        if target_y is None or target_x is None:
            return k_vol, k_geo

        # Already aligned.
        if (
            "y" in k_vol.dims
            and "x" in k_vol.dims
            and k_vol.sizes.get("y") == ref.sizes.get("y")
            and k_vol.sizes.get("x") == ref.sizes.get("x")
        ):
            return k_vol, k_geo

        try:
            return (
                k_vol.interp(y=target_y, x=target_x, method="linear"),
                k_geo.interp(y=target_y, x=target_x, method="linear"),
            )
        except (ValueError, KeyError, RuntimeError):
            # Narrowed from ``except Exception`` (REVIEW.md §2.1, §3.5
            # kernel_model.py:189). xarray.interp raises ValueError for
            # non-monotonic / unaligned coords, KeyError when a dim is
            # missing, and rasterio-backed reproject can raise
            # RuntimeError on shape mismatch — those are the "fallback
            # to shape-only resize" cases. Other exception classes
            # (e.g. MemoryError) should propagate.
            target_shape = (int(ref.sizes["y"]), int(ref.sizes["x"]))
            return (
                self._resize_kernel_grid(k_vol, target_shape, target_y, target_x),
                self._resize_kernel_grid(k_geo, target_shape, target_y, target_x),
            )

    @staticmethod
    def _resize_kernel_grid(
        da: xr.DataArray,
        target_shape: tuple[int, int],
        target_y: xr.DataArray,
        target_x: xr.DataArray,
    ) -> xr.DataArray:
        """Resize a 2-D kernel grid to target_shape while preserving output coords."""
        src = np.asarray(da.values, dtype=np.float32)
        h_out, w_out = target_shape

        if src.ndim != 2:
            return da
        if src.shape == target_shape:
            return xr.DataArray(src, dims=["y", "x"], coords={"y": target_y, "x": target_x})
        if src.shape[0] == 0 or src.shape[1] == 0:
            empty_out: FloatArray = np.full(target_shape, np.nan, dtype=np.float32)
            return xr.DataArray(empty_out, dims=["y", "x"], coords={"y": target_y, "x": target_x})

        zoom_y = h_out / src.shape[0]
        zoom_x = w_out / src.shape[1]
        out: FloatArray = np.asarray(ndimage.zoom(src, (zoom_y, zoom_x), order=1), dtype=np.float32)
        out = out[:h_out, :w_out]

        if out.shape != target_shape:
            padded: FloatArray = np.full(target_shape, np.nan, dtype=np.float32)
            h = min(out.shape[0], h_out)
            w = min(out.shape[1], w_out)
            padded[:h, :w] = out[:h, :w]
            out = padded

        return xr.DataArray(out, dims=["y", "x"], coords={"y": target_y, "x": target_x})

    def _compute_uncertainty(
        self,
        brdf_weights: BRDFKernelWeights,
        k_vol: xr.DataArray,
        k_geo: xr.DataArray,
    ) -> xr.DataArray:
        """
        Compute uncertainty in surface reflectance.

        Propagates BRDF parameter uncertainties through the kernel model.
        """
        return brdf_weights.compute_reflectance_uncertainty(k_vol, k_geo)

    def _scale_psf_sigmas(
        self,
        sigma_x: float,
        sigma_y: float,
        grid_resolution_m: float | None,
    ) -> tuple[float, float]:
        """Rescale pixel sigmas from ``target_resolution_m`` to the grid resolution.

        ``psf_sigma_x/y`` count pixels at ``self.target_resolution_m``. On a grid
        whose pixels are ``grid_resolution_m`` metres, the same physical blur
        spans ``sigma * target_resolution_m / grid_resolution_m`` pixels. When no
        grid resolution is supplied (or it is non-positive / non-finite) the
        sigmas pass through unchanged, preserving legacy behaviour.
        """
        if (
            grid_resolution_m is None
            or not np.isfinite(grid_resolution_m)
            or grid_resolution_m <= 0.0
        ):
            return sigma_x, sigma_y
        scale = self.target_resolution_m / grid_resolution_m
        return sigma_x * scale, sigma_y * scale

    def _apply_psf(
        self,
        data: xr.DataArray,
        sigma_x: float,
        sigma_y: float,
    ) -> xr.DataArray:
        """
        Apply Gaussian PSF convolution for scale matching.

        This accounts for the difference between the MODIS BRDF resolution
        (500m) and the target satellite resolution (e.g., 10m for S2).
        """
        if "band" in data.dims:
            # Batch all bands in one call with sigma=0 along the band axis.
            result = self._convolve_nd(data.values, sigma_x, sigma_y, band_axis=0)
            return xr.DataArray(result, dims=data.dims, coords=data.coords)
        else:
            result = self._convolve_2d(data.values, sigma_x, sigma_y)
            return xr.DataArray(result, dims=data.dims, coords=data.coords)

    def _convolve_nd(
        self,
        data: FloatArray,
        sigma_x: float,
        sigma_y: float,
        band_axis: int = 0,
    ) -> FloatArray:
        """Apply 2D Gaussian convolution to each band in a 3-D cube.

        Uses in-place operations to reduce peak memory from ~5x to ~3x input size.
        """
        data_arr = np.asarray(data, dtype=np.float32)
        mask = np.isfinite(data_arr)

        # Fill NaN in-place on a working copy (avoids separate data_filled allocation).
        result: FloatArray = np.where(mask, data_arr, 0.0).astype(np.float32, copy=False)

        # Build per-axis sigma: 0 for band axis, sigma_y/sigma_x for spatial.
        sigma_nd = [0.0] * data_arr.ndim
        spatial = [i for i in range(data_arr.ndim) if i != band_axis]
        sigma_nd[spatial[0]] = sigma_y
        sigma_nd[spatial[1]] = sigma_x

        # Convolve data in-place via output parameter.
        ndimage.gaussian_filter(result, sigma=sigma_nd, mode="reflect", output=result)
        # Convolve mask (reuse mask array to avoid extra allocation).
        norm = mask.astype(np.float32, copy=False)
        ndimage.gaussian_filter(norm, sigma=sigma_nd, mode="reflect", output=norm)

        # Normalize in-place, then mask invalid regions.
        with np.errstate(divide="ignore", invalid="ignore"):
            np.maximum(norm, 1.0e-10, out=norm)
            np.divide(result, norm, out=result)
        result[norm < (1.0e-10 + 0.01)] = np.nan
        return result

    def _convolve_2d(
        self,
        data: FloatArray,
        sigma_x: float,
        sigma_y: float,
    ) -> FloatArray:
        """Apply 2D Gaussian convolution.

        Uses in-place operations to reduce peak memory from ~5x to ~3x input size.
        """
        data_arr = np.asarray(data, dtype=np.float32)
        mask = np.isfinite(data_arr)

        # Fill NaN in-place on a working copy.
        result: FloatArray = np.where(mask, data_arr, 0.0).astype(np.float32, copy=False)

        # Convolve data in-place via output parameter.
        ndimage.gaussian_filter(result, sigma=[sigma_y, sigma_x], mode="reflect", output=result)
        # Convolve mask (reuse array to avoid extra allocation).
        norm = mask.astype(np.float32, copy=False)
        ndimage.gaussian_filter(norm, sigma=[sigma_y, sigma_x], mode="reflect", output=norm)

        # Normalize in-place, then mask invalid regions.
        with np.errstate(divide="ignore", invalid="ignore"):
            np.maximum(norm, 1.0e-10, out=norm)
            np.divide(result, norm, out=result)
        result[norm < (1.0e-10 + 0.01)] = np.nan
        return result


class PSFConvolver:
    """
    PSF convolution for scale matching using DCT.

    More efficient than direct convolution for large kernels.
    Uses discrete cosine transform (DCT) for frequency-domain
    convolution.
    """

    def __init__(self, sigma_x: float = 29.75, sigma_y: float = 39.0) -> None:
        self.sigma_x = sigma_x
        self.sigma_y = sigma_y

    def convolve(self, data: np.ndarray | xr.DataArray) -> np.ndarray | xr.DataArray:
        """
        Apply PSF convolution using DCT.

        Args:
            data: 2D array to convolve

        Returns:
            Convolved array
        """
        if isinstance(data, xr.DataArray):
            arr = np.asarray(data.values, dtype=np.float32)
            result = self._dct_convolve(arr)
            return xr.DataArray(result, dims=data.dims, coords=data.coords)

        arr = np.asarray(data, dtype=np.float32)
        return self._dct_convolve(arr)

    def _dct_convolve(self, data: FloatArray) -> FloatArray:
        """DCT-based Gaussian convolution."""
        from scipy.fftpack import dct, idct

        data_arr = np.asarray(data, dtype=np.float32)
        height, width = data_arr.shape

        # Handle NaN
        mask = np.isfinite(data_arr)
        data_filled = np.where(mask, data_arr, 0.0).astype(np.float32, copy=False)

        # Gaussian in frequency domain
        u = np.arange(height, dtype=np.float32) / max(height, 1)
        v = np.arange(width, dtype=np.float32) / max(width, 1)

        gx = np.exp(-2 * np.pi**2 * self.sigma_x**2 * (0.5 * u) ** 2).astype(np.float32)
        gy = np.exp(-2 * np.pi**2 * self.sigma_y**2 * (0.5 * v) ** 2).astype(np.float32)
        kernel: FloatArray = np.asarray(np.outer(gx, gy), dtype=np.float32)

        # Forward DCT
        data_dct = np.asarray(
            dct(dct(data_filled, axis=0, norm="ortho"), axis=1, norm="ortho"),
            dtype=np.float32,
        )
        mask_dct = np.asarray(
            dct(dct(mask.astype(np.float32), axis=0, norm="ortho"), axis=1, norm="ortho"),
            dtype=np.float32,
        )

        # Multiply in frequency domain
        convolved_dct = data_dct * kernel
        norm_dct = mask_dct * kernel

        # Inverse DCT
        convolved = np.asarray(
            idct(idct(convolved_dct, axis=1, norm="ortho"), axis=0, norm="ortho"),
            dtype=np.float32,
        )
        norm = np.asarray(
            idct(idct(norm_dct, axis=1, norm="ortho"), axis=0, norm="ortho"),
            dtype=np.float32,
        )

        # Normalize
        with np.errstate(divide="ignore", invalid="ignore"):
            result: FloatArray = np.asarray(
                convolved / np.maximum(norm, 1.0e-10),
                dtype=np.float32,
            )

        restored: FloatArray = np.asarray(np.where(norm > 0.01, result, np.nan), dtype=np.float32)
        return restored
