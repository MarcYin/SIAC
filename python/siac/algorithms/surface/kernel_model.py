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
    ) -> SurfacePrior:
        """
        Compute surface reflectance prior from BRDF parameters.

        Args:
            brdf_weights: BRDF kernel coefficients (f0, f1, f2)
            geometry: Observation geometry (sza, vza, raa)
            psf_params: Optional (sigma_x, sigma_y) override for PSF

        Returns:
            SurfacePrior with BOA reflectance, uncertainty, and mask
        """
        if psf_params is not None:
            sigma_x, sigma_y = psf_params
        else:
            sigma_x, sigma_y = self.psf_sigma_x, self.psf_sigma_y

        # Compute kernel values at observation geometry
        k_vol, k_geo = self._kernels.compute(
            geometry.vza, geometry.sza, geometry.raa
        )
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

        # Create validity mask
        mask = (boa > 0) & (boa < 1.5) & np.isfinite(boa)

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
            "y" in k_vol.dims and "x" in k_vol.dims
            and k_vol.sizes.get("y") == ref.sizes.get("y")
            and k_vol.sizes.get("x") == ref.sizes.get("x")
        ):
            return k_vol, k_geo

        try:
            return (
                k_vol.interp(y=target_y, x=target_x, method="linear"),
                k_geo.interp(y=target_y, x=target_x, method="linear"),
            )
        except Exception:
            # Fallback to shape-only resize when coordinates are not monotonic/aligned.
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
        """Apply 2D Gaussian convolution to each band in a 3-D cube."""
        data_arr = np.asarray(data, dtype=np.float32)
        mask = np.isfinite(data_arr)
        data_filled = np.where(mask, data_arr, 0.0).astype(np.float32, copy=False)

        # Build per-axis sigma: 0 for band axis, sigma_y/sigma_x for spatial.
        sigma_nd = [0.0] * data_arr.ndim
        sigma_nd[band_axis] = 0.0
        # Spatial axes are the remaining two (in order: y, x).
        spatial = [i for i in range(data_arr.ndim) if i != band_axis]
        sigma_nd[spatial[0]] = sigma_y
        sigma_nd[spatial[1]] = sigma_x

        convolved = np.asarray(
            ndimage.gaussian_filter(data_filled, sigma=sigma_nd, mode="reflect"),
            dtype=np.float32,
        )
        norm = np.asarray(
            ndimage.gaussian_filter(mask.astype(np.float32), sigma=sigma_nd, mode="reflect"),
            dtype=np.float32,
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            result: FloatArray = np.asarray(convolved / np.maximum(norm, 1.0e-10), dtype=np.float32)
        restored: FloatArray = np.asarray(np.where(norm > 0.01, result, np.nan), dtype=np.float32)
        return restored

    def _convolve_2d(
        self,
        data: FloatArray,
        sigma_x: float,
        sigma_y: float,
    ) -> FloatArray:
        """Apply 2D Gaussian convolution."""
        data_arr = np.asarray(data, dtype=np.float32)
        # Handle NaN values
        mask = np.isfinite(data_arr)
        data_filled = np.where(mask, data_arr, 0.0).astype(np.float32, copy=False)

        # Convolve data and mask
        convolved = np.asarray(
            ndimage.gaussian_filter(
                data_filled, sigma=[sigma_y, sigma_x], mode="reflect"
            ),
            dtype=np.float32,
        )
        norm = np.asarray(
            ndimage.gaussian_filter(
                mask.astype(np.float32), sigma=[sigma_y, sigma_x], mode="reflect"
            ),
            dtype=np.float32,
        )

        # Normalize to account for missing data
        with np.errstate(divide="ignore", invalid="ignore"):
            result: FloatArray = np.asarray(
                convolved / np.maximum(norm, 1.0e-10),
                dtype=np.float32,
            )

        # Restore NaN where no valid data
        restored: FloatArray = np.asarray(np.where(norm > 0.01, result, np.nan), dtype=np.float32)
        return restored


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
