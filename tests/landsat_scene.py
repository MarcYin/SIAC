"""Synthetic Landsat Collection-2 Level-1 files shared by unit and seam tests."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np
import rasterio
from rasterio.transform import from_origin

from siac.adapters.satellite.landsat import parse_mtl_text

if TYPE_CHECKING:
    from pathlib import Path

PRODUCT_ID = "LC09_L1TP_199024_20240615_20240615_02_T1"
HEIGHT, WIDTH = 40, 50
ORIGIN_X, ORIGIN_Y = 600000.0, 5700000.0
CRS = "EPSG:32631"
MULT, ADD = 2.0e-5, -0.1
SZA_DEG = 35.0

# QA_PIXEL test values (Collection-2 bit layout).
QA_FILL = 1 << 0
QA_CLEAR = 1 << 6
QA_CLOUD = (1 << 3) | (1 << 1)
QA_SHADOW = 1 << 4
QA_CIRRUS = 1 << 2


def _write(path: Path, data: np.ndarray, dtype: str) -> None:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype=dtype,
        crs=CRS,
        transform=from_origin(ORIGIN_X, ORIGIN_Y, 30.0, 30.0),
    ) as dst:
        dst.write(data.astype(dtype), 1)


def _mtl_text(spacecraft: str, level: str, *, with_angle_names: bool) -> str:
    files = "\n".join(f'    FILE_NAME_BAND_{n} = "{PRODUCT_ID}_B{n}.TIF"' for n in range(1, 8))
    angles = ""
    if with_angle_names:
        angles = (
            f'    FILE_NAME_ANGLE_SOLAR_ZENITH_BAND_4 = "{PRODUCT_ID}_SZA.TIF"\n'
            f'    FILE_NAME_ANGLE_SOLAR_AZIMUTH_BAND_4 = "{PRODUCT_ID}_SAA.TIF"\n'
            f'    FILE_NAME_ANGLE_SENSOR_ZENITH_BAND_4 = "{PRODUCT_ID}_VZA.TIF"\n'
            f'    FILE_NAME_ANGLE_SENSOR_AZIMUTH_BAND_4 = "{PRODUCT_ID}_VAA.TIF"\n'
            f'    FILE_NAME_QUALITY_L1_PIXEL = "{PRODUCT_ID}_QA_PIXEL.TIF"\n'
        )
    rescale = "\n".join(
        f"    REFLECTANCE_MULT_BAND_{n} = {MULT:.4E}\n    REFLECTANCE_ADD_BAND_{n} = {ADD:.6f}"
        for n in range(1, 10)
    )
    return f"""GROUP = LANDSAT_METADATA_FILE
  GROUP = PRODUCT_CONTENTS
    LANDSAT_PRODUCT_ID = "{PRODUCT_ID}"
    COLLECTION_NUMBER = 02
    COLLECTION_CATEGORY = "T1"
    PROCESSING_LEVEL = "{level}"
{files}
{angles}  END_GROUP = PRODUCT_CONTENTS
  GROUP = IMAGE_ATTRIBUTES
    SPACECRAFT_ID = "{spacecraft}"
    SENSOR_ID = "OLI_TIRS"
    WRS_PATH = 199
    WRS_ROW = 24
    DATE_ACQUIRED = 2024-06-15
    SCENE_CENTER_TIME = "10:47:12.3456789Z"
    CLOUD_COVER = 12.5
    SUN_AZIMUTH = 150.0
    SUN_ELEVATION = {90.0 - SZA_DEG}
    EARTH_SUN_DISTANCE = 1.0158
  END_GROUP = IMAGE_ATTRIBUTES
  GROUP = LEVEL1_RADIOMETRIC_RESCALING
{rescale}
  END_GROUP = LEVEL1_RADIOMETRIC_RESCALING
END_GROUP = LANDSAT_METADATA_FILE
END
"""


def _qa_grid() -> np.ndarray:
    qa = np.full((HEIGHT, WIDTH), QA_CLEAR, dtype=np.uint16)
    qa[:5, :] = QA_FILL
    qa[10:15, 10:20] = QA_CLOUD
    qa[20:25, 10:20] = QA_SHADOW
    qa[30:35, 30:40] = QA_CIRRUS
    return qa


def _make_scene(
    root: Path,
    *,
    spacecraft: str = "LANDSAT_9",
    level: str = "L1TP",
    with_angles: bool = True,
    with_angle_names: bool = True,
    mtl_format: str = "txt",
) -> Path:
    scene = root / PRODUCT_ID
    scene.mkdir()
    fill_rows = slice(0, 5)
    for n in range(1, 8):
        dn = np.full((HEIGHT, WIDTH), 10000 + 1000 * n, dtype=np.uint16)
        dn[fill_rows, :] = 0
        _write(scene / f"{PRODUCT_ID}_B{n}.TIF", dn, "uint16")
    if with_angles:
        sza = np.full((HEIGHT, WIDTH), SZA_DEG * 100)
        saa = np.full((HEIGHT, WIDTH), 150.0 * 100)
        vza = np.tile(np.linspace(-700, 700, WIDTH), (HEIGHT, 1))
        vza = np.abs(vza)
        # Azimuths are stored in [-180, 180] deg; alternate -179 / +179 deg inside each block.
        vaa = np.where(np.arange(WIDTH)[None, :] % 2 == 0, -17900, 17900) * np.ones((HEIGHT, 1))
        for name, grid in (("SZA", sza), ("SAA", saa), ("VZA", vza), ("VAA", vaa)):
            grid = grid.copy()
            grid[fill_rows, :] = 0
            _write(scene / f"{PRODUCT_ID}_{name}.TIF", grid, "int16")
    _write(scene / f"{PRODUCT_ID}_QA_PIXEL.TIF", _qa_grid(), "uint16")

    text = _mtl_text(spacecraft, level, with_angle_names=with_angle_names and with_angles)
    if mtl_format == "json":
        payload = json.dumps(parse_mtl_text(text))
        (scene / f"{PRODUCT_ID}_MTL.json").write_text(payload, encoding="utf-8")
    else:
        (scene / f"{PRODUCT_ID}_MTL.txt").write_text(text, encoding="utf-8")
    return scene
