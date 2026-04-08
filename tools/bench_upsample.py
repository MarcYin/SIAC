"""Benchmark angle grid upsampling methods."""
import numpy as np
import time

src = np.random.rand(23, 23).astype(np.float32)
h_out, w_out = 10980, 10980

# Method 1: scipy.ndimage.zoom (current)
from scipy.ndimage import zoom
_ = zoom(src, (h_out / 23, w_out / 23), order=1)  # warm up
t0 = time.perf_counter()
for _ in range(4):
    out = zoom(src, (h_out / 23, w_out / 23), order=1)
print(f"scipy.zoom: {(time.perf_counter() - t0) / 4:.3f}s per call, shape={out.shape}")

# Method 2: PIL/Pillow BILINEAR
from PIL import Image
img = Image.fromarray(src, mode="F")
_ = np.array(img.resize((w_out, h_out), Image.BILINEAR))  # warm up
t0 = time.perf_counter()
for _ in range(4):
    img = Image.fromarray(src, mode="F")
    out2 = np.array(img.resize((w_out, h_out), Image.BILINEAR))
print(f"PIL resize: {(time.perf_counter() - t0) / 4:.3f}s per call, shape={out2.shape}")

# Method 3: scipy.interpolate.RegularGridInterpolator
from scipy.interpolate import RegularGridInterpolator
yy_src = np.arange(src.shape[0], dtype=np.float64)
xx_src = np.arange(src.shape[1], dtype=np.float64)
interp = RegularGridInterpolator((yy_src, xx_src), src.astype(np.float64), method="linear", bounds_error=False, fill_value=None)
yy_out = np.linspace(0, src.shape[0] - 1, h_out)
xx_out = np.linspace(0, src.shape[1] - 1, w_out)
yy_grid, xx_grid = np.meshgrid(yy_out, xx_out, indexing="ij")
pts = np.column_stack([yy_grid.ravel(), xx_grid.ravel()])
_ = interp(pts).reshape(h_out, w_out)  # warm up
t0 = time.perf_counter()
for _ in range(2):
    out3 = interp(pts).reshape(h_out, w_out).astype(np.float32)
print(f"RegularGridInterp: {(time.perf_counter() - t0) / 2:.3f}s per call, shape={out3.shape}")

# Method 4: map_coordinates with precomputed coords
from scipy.ndimage import map_coordinates
yy_coords = np.linspace(0, src.shape[0] - 1, h_out, dtype=np.float64)
xx_coords = np.linspace(0, src.shape[1] - 1, w_out, dtype=np.float64)
yy_g, xx_g = np.meshgrid(yy_coords, xx_coords, indexing="ij")
coords = np.array([yy_g, xx_g])
_ = map_coordinates(src, coords, order=1, mode="nearest")  # warm up
t0 = time.perf_counter()
for _ in range(4):
    out4 = map_coordinates(src, coords, order=1, mode="nearest").astype(np.float32)
print(f"map_coordinates: {(time.perf_counter() - t0) / 4:.3f}s per call, shape={out4.shape}")

# Check max difference
print(f"\nMax diff PIL vs scipy.zoom: {np.max(np.abs(out.astype(np.float32) - out2)):.6f}")
print(f"Max diff map_coords vs scipy.zoom: {np.max(np.abs(out.astype(np.float32) - out4)):.6f}")
