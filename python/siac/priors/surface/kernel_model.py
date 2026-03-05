"""
Surface prior derivation using BRDF kernel model.

Computes surface reflectance prior from MODIS MCD43 BRDF parameters
using the Ross-Thick Li-Sparse kernel model, with optional PSF
convolution for scale matching.
"""

from __future__ import annotations

import logging

import numpy as np
import xarray as xr
from scipy import ndimage

from siac.core.types import BRDFKernelWeights, GeometryAngles, SurfacePrior
from siac.priors.brdf.kernels import BRDFKernels, compute_reflectance

logger = logging.getLogger(__name__)


class KernelModelDeriver:
    """
    Derive surface reflectance prior from BRDF kernel weights.

    Uses the Ross-Thick Li-Sparse BRDF model to compute expected
    surface reflectance at the observation geometry.

    Args:
        psf_sigma_x: PSF standard deviation in x direction (pixels)
        psf_sigma_y: PSF standard deviation in y direction (pixels)
        apply_psf: Whether to apply PSF convolution for scale matching
    """

    def __init__(
        self,
        psf_sigma_x: float = 29.75,
        psf_sigma_y: float = 39.0,
        apply_psf: bool = True,
    ):
        self.psf_sigma_x = psf_sigma_x
        self.psf_sigma_y = psf_sigma_y
        self.apply_psf = apply_psf
        self._kernels = BRDFKernels(hb=2.0, br=1.0)

    def compute_surface_prior(
        self,
        brdf_weights: BRDFKernelWeights,
        geometry: GeometryAngles,
        psf_params: tuple[float, float] | None = None,
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
        k_vol, k_geo = self._align_kernels_to_brdf_grid(k_vol, k_geo, brdf_weights.f0)

        # Compute reflectance per band
        boa = compute_reflectance(
            brdf_weights.f0,
            brdf_weights.f1,
            brdf_weights.f2,
            k_vol,
            k_geo,
        )

        # Compute uncertainty
        boa_unc = self._compute_uncertainty(
            brdf_weights, k_vol, k_geo
        )

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
            out = np.full(target_shape, np.nan, dtype=np.float32)
            return xr.DataArray(out, dims=["y", "x"], coords={"y": target_y, "x": target_x})

        zoom_y = h_out / src.shape[0]
        zoom_x = w_out / src.shape[1]
        out = ndimage.zoom(src, (zoom_y, zoom_x), order=1)
        out = out[:h_out, :w_out]

        if out.shape != target_shape:
            padded = np.full(target_shape, np.nan, dtype=np.float32)
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
        # Variance from each component
        var_f0 = brdf_weights.f0_unc ** 2
        var_f1 = (k_vol * brdf_weights.f1_unc) ** 2
        var_f2 = (k_geo * brdf_weights.f2_unc) ** 2

        # Total variance
        var_total = var_f0 + var_f1 + var_f2

        return np.sqrt(var_total)

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
            # Process each band separately
            result_list = []
            for band in data.band.values:
                band_data = data.sel(band=band)
                convolved = self._convolve_2d(band_data.values, sigma_x, sigma_y)
                result_list.append(convolved)

            result = np.stack(result_list, axis=0)
            return xr.DataArray(
                result,
                dims=data.dims,
                coords=data.coords,
            )
        else:
            result = self._convolve_2d(data.values, sigma_x, sigma_y)
            return xr.DataArray(
                result,
                dims=data.dims,
                coords=data.coords,
            )

    def _convolve_2d(
        self,
        data: np.ndarray,
        sigma_x: float,
        sigma_y: float,
    ) -> np.ndarray:
        """Apply 2D Gaussian convolution."""
        # Handle NaN values
        mask = np.isfinite(data)
        data_filled = np.where(mask, data, 0.0)

        # Convolve data and mask
        convolved = ndimage.gaussian_filter(
            data_filled, sigma=[sigma_y, sigma_x], mode="reflect"
        )
        norm = ndimage.gaussian_filter(
            mask.astype(float), sigma=[sigma_y, sigma_x], mode="reflect"
        )

        # Normalize to account for missing data
        with np.errstate(divide="ignore", invalid="ignore"):
            result = convolved / np.maximum(norm, 1e-10)

        # Restore NaN where no valid data
        result = np.where(norm > 0.01, result, np.nan)

        return result


class PSFConvolver:
    """
    PSF convolution for scale matching using DCT.

    More efficient than direct convolution for large kernels.
    Uses discrete cosine transform (DCT) for frequency-domain
    convolution.
    """

    def __init__(self, sigma_x: float = 29.75, sigma_y: float = 39.0):
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
        is_xarray = isinstance(data, xr.DataArray)
        if is_xarray:
            arr = data.values
            template = data
        else:
            arr = np.asarray(data)
            template = None

        # Apply DCT-based convolution
        result = self._dct_convolve(arr)

        if is_xarray:
            return xr.DataArray(result, dims=template.dims, coords=template.coords)
        return result

    def _dct_convolve(self, data: np.ndarray) -> np.ndarray:
        """DCT-based Gaussian convolution."""
        from scipy.fftpack import dct, idct

        height, width = data.shape

        # Handle NaN
        mask = np.isfinite(data)
        data_filled = np.where(mask, data, 0.0)

        # Gaussian in frequency domain
        u = np.arange(height) / height
        v = np.arange(width) / width

        gx = np.exp(-2 * np.pi**2 * self.sigma_x**2 * (0.5 * u) ** 2)
        gy = np.exp(-2 * np.pi**2 * self.sigma_y**2 * (0.5 * v) ** 2)
        kernel = np.outer(gx, gy)

        # Forward DCT
        data_dct = dct(dct(data_filled, axis=0, norm="ortho"), axis=1, norm="ortho")
        mask_dct = dct(dct(mask.astype(float), axis=0, norm="ortho"), axis=1, norm="ortho")

        # Multiply in frequency domain
        convolved_dct = data_dct * kernel
        norm_dct = mask_dct * kernel

        # Inverse DCT
        convolved = idct(idct(convolved_dct, axis=1, norm="ortho"), axis=0, norm="ortho")
        norm = idct(idct(norm_dct, axis=1, norm="ortho"), axis=0, norm="ortho")

        # Normalize
        with np.errstate(divide="ignore", invalid="ignore"):
            result = convolved / np.maximum(norm, 1e-10)

        result = np.where(norm > 0.01, result, np.nan)

        return result
