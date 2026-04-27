"""
Layer 3 — Grid assembler (M4) unit tests.
"""

import dataclasses
import warnings

import numpy as np
import pytest
import xarray as xr
from rasterio.enums import Resampling as RasterioResampling

import siac.algorithms.grid.assembler as assembler_mod
from siac.algorithms.grid.assembler import assemble_grids
from siac.catalog import SENTINEL2A_CONFIG
from siac.config.schema import SharpTransitionFilterConfig
from siac.domain import SensorBand, SensorConfig
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.validation import validate_solver_input_bundle


@pytest.fixture
def large_obs_bundle():
    """ObservationBundle at 64x64 (simulating ~10 m native) with 4 bands."""
    shape = (64, 64)
    config = SensorConfig(
        sensor_id="MOCK",
        satellite_id="TEST",
        bands=(
            SensorBand("B01", 443.0, 20.0, 60.0, 0),
            SensorBand("B02", 490.0, 65.0, 10.0, 1),
            SensorBand("B03", 560.0, 35.0, 10.0, 2),
            SensorBand("B04", 665.0, 30.0, 10.0, 3),
        ),
    )
    from datetime import datetime

    toa = xr.Dataset(
        {
            "B01": xr.DataArray(
                np.random.RandomState(10).uniform(0.05, 0.3, shape).astype(np.float32),
                dims=["y", "x"],
            ),
            "B02": xr.DataArray(
                np.random.RandomState(11).uniform(0.05, 0.3, shape).astype(np.float32),
                dims=["y", "x"],
            ),
            "B03": xr.DataArray(
                np.random.RandomState(12).uniform(0.05, 0.3, shape).astype(np.float32),
                dims=["y", "x"],
            ),
            "B04": xr.DataArray(
                np.random.RandomState(13).uniform(0.05, 0.3, shape).astype(np.float32),
                dims=["y", "x"],
            ),
        }
    )
    geometry = GeometryAngles(
        sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
        saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
    )
    cloud = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
    return ObservationBundle(
        toa=toa,
        geometry=geometry,
        cloud_mask=cloud,
        sensor_config=config,
        metadata={"observation_time": datetime(2023, 7, 15, 10, 30)},
        crs="EPSG:32632",
        bounds=(300000.0, 5500000.0, 300640.0, 5500640.0),
    )


@pytest.fixture
def large_atmo():
    """AtmosphericState at 64x64."""
    shape = (64, 64)
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
    )


@pytest.fixture
def large_surface():
    """SurfacePrior at 64x64."""
    shape = (64, 64)
    brdf = BRDFKernelWeights(
        f0=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        f1=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
        f2=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        f0_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
        f1_unc=xr.DataArray(np.full(shape, 0.005), dims=["y", "x"]),
        f2_unc=xr.DataArray(np.full(shape, 0.002), dims=["y", "x"]),
    )
    return SurfacePrior(
        boa=xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
        boa_unc=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
        kernels=brdf,
        mask=xr.DataArray(np.ones(shape, dtype=bool), dims=["y", "x"]),
    )


class TestAssembleGrids:
    def test_output_type(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        result = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        assert isinstance(result, SolverInputBundle)

    def test_spatial_alignment(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """All raster fields should share the same (y, x) shape."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        # toa has (band, y, x)
        toa_spatial = sib.toa.shape[1:]
        assert sib.geometry.sza.shape == toa_spatial
        assert sib.cloud_mask.shape == toa_spatial
        assert sib.atmo_prior.aot.shape == toa_spatial
        assert sib.surface_prior.boa.shape == toa_spatial

    def test_aux_resolution(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        sib = assemble_grids(
            large_obs_bundle, large_atmo, large_surface, mock_rt_model, aux_resolution_m=500.0
        )
        assert sib.aux_resolution_m == 500.0

    def test_aux_resolution_falls_back_to_solver_grid_for_legacy_callers(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ):
        sib = assemble_grids(
            large_obs_bundle, large_atmo, large_surface, mock_rt_model, aux_resolution_m=320.0
        )

        assert sib.toa.shape[1:] == (2, 2)
        assert sib.aerosol_resolution_m == 320.0

    def test_band_selection(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Non-MSI sensors should keep the wavelength-based aerosol defaults."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        assert [b.name for b in sib.bands] == ["B01", "B02"]

    def test_band_selection_includes_requested_stage_bands_in_sensor_order(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            solver_band_names=("B04", "B02", "B01"),
        )

        assert [band.name for band in sib.bands] == ["B01", "B02", "B04"]
        assert list(np.asarray(sib.toa.coords["band"].values).tolist()) == ["B01", "B02", "B04"]

    def test_surface_prior_is_aligned_to_solver_band_names(self, large_atmo, mock_rt_model):
        shape = (4, 4)
        toa = xr.Dataset(
            {
                "B01": xr.DataArray(np.full(shape, 0.08, dtype=np.float32), dims=["y", "x"]),
                "B02": xr.DataArray(np.full(shape, 0.12, dtype=np.float32), dims=["y", "x"]),
                "B04": xr.DataArray(np.full(shape, 0.18, dtype=np.float32), dims=["y", "x"]),
            }
        )
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5, dtype=np.float32), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5, dtype=np.float32), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5, dtype=np.float32), dims=["y", "x"]),
        )
        obs = ObservationBundle(
            toa=toa,
            geometry=geometry,
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
            sensor_config=SensorConfig(
                sensor_id="MSI",
                satellite_id="S2A",
                bands=tuple(SENTINEL2A_CONFIG.get_band(name) for name in ("B01", "B02", "B04")),
            ),
            metadata={"observation_time": None},
            crs="EPSG:32632",
            bounds=(300000.0, 5500000.0, 300040.0, 5500040.0),
        )
        surface = SurfacePrior(
            boa=xr.DataArray(
                np.stack(
                    [
                        np.full(shape, 0.44, dtype=np.float32),
                        np.full(shape, 0.11, dtype=np.float32),
                        np.full(shape, 0.22, dtype=np.float32),
                    ],
                    axis=0,
                ),
                dims=["band", "y", "x"],
                coords={"band": ["B04", "B01", "B02"]},
            ),
            boa_unc=xr.DataArray(
                np.stack(
                    [
                        np.full(shape, 0.04, dtype=np.float32),
                        np.full(shape, 0.01, dtype=np.float32),
                        np.full(shape, 0.02, dtype=np.float32),
                    ],
                    axis=0,
                ),
                dims=["band", "y", "x"],
                coords={"band": ["B04", "B01", "B02"]},
            ),
            kernels=None,
            mask=xr.DataArray(
                np.ones((3, *shape), dtype=bool),
                dims=["band", "y", "x"],
                coords={"band": ["B04", "B01", "B02"]},
            ),
        )

        sib = assemble_grids(obs, large_atmo, surface, mock_rt_model, aerosol_resolution_m=10.0)

        assert [band.name for band in sib.bands] == ["B02", "B04"]
        assert list(sib.surface_prior.boa.coords["band"].values) == ["B02", "B04"]
        np.testing.assert_allclose(sib.surface_prior.boa.sel(band="B02").values, 0.22)
        np.testing.assert_allclose(sib.surface_prior.boa.sel(band="B04").values, 0.44)

    def test_cloud_mask_conservative(
        self, large_obs_bundle, large_atmo, large_surface, mock_rt_model
    ):
        """Any native pixel that is cloud → aux pixel is cloud."""
        import dataclasses

        shape = (64, 64)
        cloud = np.zeros(shape, dtype=bool)
        cloud[10:20, 10:20] = True
        obs = dataclasses.replace(
            large_obs_bundle,
            cloud_mask=xr.DataArray(cloud, dims=["y", "x"]),
        )
        sib = assemble_grids(obs, large_atmo, large_surface, mock_rt_model)
        # The cloud region should still have True values after resampling
        assert sib.cloud_mask.values.any(), "Cloud region should be preserved"

    def test_water_mask_is_loaded_and_aggregated_to_solver_grid(
        self,
        monkeypatch,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ):
        native = xr.DataArray(np.zeros((64, 64), dtype=bool), dims=["y", "x"])
        native.values[10:20, 10:20] = True

        captured: dict[str, object] = {}

        def _fake_load(bounds, crs, *, source=None, cache_dir=None, session=None):  # noqa: ANN001
            del session
            captured["bounds"] = bounds
            captured["crs"] = crs
            captured["source"] = source
            captured["cache_dir"] = cache_dir
            return native

        monkeypatch.setattr(assembler_mod, "load_water_mask_subset", _fake_load)

        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            water_mask_path="https://example.com/landWater2020.vrt",
            water_mask_cache_dir="/tmp/water-mask-cache",
        )

        assert sib.water_mask is not None
        assert sib.water_mask.values.any()
        assert captured["source"] == "https://example.com/landWater2020.vrt"
        assert captured["cache_dir"] == "/tmp/water-mask-cache"

    def test_water_mask_buffer_is_applied_before_resampling(
        self,
        monkeypatch,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ):
        native = xr.DataArray(np.zeros((64, 64), dtype=bool), dims=["y", "x"])
        native.values[20, 20] = True

        def _fake_load(_bounds, _crs, *, source=None, cache_dir=None, session=None):  # noqa: ANN001
            del source, cache_dir, session
            return native

        monkeypatch.setattr(assembler_mod, "load_water_mask_subset", _fake_load)

        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=10.0,
            water_mask_path="https://example.com/landWater2020.vrt",
            water_mask_buffer_pixels=1,
        )

        expected = assembler_mod._dilate_mask(native.values, 1)
        np.testing.assert_array_equal(sib.water_mask.values, expected)

    def test_water_mask_is_geospatially_reprojected_to_solver_grid(
        self,
        monkeypatch,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        import rioxarray  # noqa: F401
        from rasterio.transform import from_bounds
        from rasterio.warp import transform_bounds

        lon_min, lat_min, lon_max, lat_max = transform_bounds(
            large_obs_bundle.crs,
            "EPSG:4326",
            *large_obs_bundle.bounds,
        )
        height = width = 64
        x = lon_min + (np.arange(width, dtype=np.float64) + 0.5) * ((lon_max - lon_min) / width)
        y = lat_max - (np.arange(height, dtype=np.float64) + 0.5) * ((lat_max - lat_min) / height)
        native = xr.DataArray(
            np.zeros((height, width), dtype=bool),
            dims=["y", "x"],
            coords={"x": x, "y": y},
        )
        native.values[:, : width // 2] = True
        native = native.rio.set_spatial_dims(x_dim="x", y_dim="y")
        native = native.rio.write_crs("EPSG:4326")
        native = native.rio.write_transform(
            from_bounds(lon_min, lat_min, lon_max, lat_max, width, height)
        )

        def _fake_load(_bounds, _crs, *, source=None, cache_dir=None, session=None):  # noqa: ANN001
            del _bounds, _crs, source, cache_dir, session
            return native

        monkeypatch.setattr(assembler_mod, "load_water_mask_subset", _fake_load)

        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=320.0,
            water_mask_path="https://example.com/landWater2020.vrt",
        )

        target_template = assembler_mod._build_target_template(
            large_obs_bundle.bounds,
            large_obs_bundle.crs,
            320.0,
        )
        expected_fraction = assembler_mod._resample_da(
            native.astype(np.float32),
            target_template.shape,
            "area",
            template=target_template,
        )
        expected = expected_fraction.values > 0.0

        np.testing.assert_array_equal(sib.water_mask.values, expected)

    def test_shape_only_mask_remap_requires_matching_crs_and_bounds(
        self,
        large_obs_bundle,
    ) -> None:
        import rioxarray  # noqa: F401
        from rasterio.transform import from_bounds
        from rasterio.warp import transform_bounds

        target_template = assembler_mod._build_target_template(
            large_obs_bundle.bounds,
            large_obs_bundle.crs,
            320.0,
        )

        lon_min, lat_min, lon_max, lat_max = transform_bounds(
            large_obs_bundle.crs,
            "EPSG:4326",
            *large_obs_bundle.bounds,
        )
        wgs84_x = lon_min + (np.arange(64, dtype=np.float64) + 0.5) * ((lon_max - lon_min) / 64.0)
        wgs84_y = lat_max - (np.arange(64, dtype=np.float64) + 0.5) * ((lat_max - lat_min) / 64.0)
        wgs84_mask = xr.DataArray(
            np.zeros((64, 64), dtype=bool),
            dims=["y", "x"],
            coords={"x": wgs84_x, "y": wgs84_y},
        )
        wgs84_mask = wgs84_mask.rio.set_spatial_dims(x_dim="x", y_dim="y")
        wgs84_mask = wgs84_mask.rio.write_crs("EPSG:4326")
        wgs84_mask = wgs84_mask.rio.write_transform(
            from_bounds(lon_min, lat_min, lon_max, lat_max, 64, 64)
        )

        assert not assembler_mod._shape_only_mask_remap_is_safe(wgs84_mask, target_template)

        utm_xmin, utm_ymin, utm_xmax, utm_ymax = large_obs_bundle.bounds
        utm_x = utm_xmin + (np.arange(64, dtype=np.float64) + 0.5) * ((utm_xmax - utm_xmin) / 64.0)
        utm_y = utm_ymax - (np.arange(64, dtype=np.float64) + 0.5) * ((utm_ymax - utm_ymin) / 64.0)
        utm_mask = xr.DataArray(
            np.zeros((64, 64), dtype=bool),
            dims=["y", "x"],
            coords={"x": utm_x, "y": utm_y},
        )
        utm_mask = utm_mask.rio.set_spatial_dims(x_dim="x", y_dim="y")
        utm_mask = utm_mask.rio.write_crs(large_obs_bundle.crs)
        utm_mask = utm_mask.rio.write_transform(
            from_bounds(utm_xmin, utm_ymin, utm_xmax, utm_ymax, 64, 64)
        )

        assert assembler_mod._shape_only_mask_remap_is_safe(utm_mask, target_template)

    def test_cloud_mask_geospatial_fallback_when_source_is_not_aligned(
        self,
        large_obs_bundle,
    ) -> None:
        import rioxarray  # noqa: F401
        from rasterio.transform import from_bounds
        from rasterio.warp import transform_bounds

        target_template = assembler_mod._build_target_template(
            large_obs_bundle.bounds,
            large_obs_bundle.crs,
            320.0,
        )

        lon_min, lat_min, lon_max, lat_max = transform_bounds(
            large_obs_bundle.crs,
            "EPSG:4326",
            *large_obs_bundle.bounds,
        )
        x = lon_min + (np.arange(64, dtype=np.float64) + 0.5) * ((lon_max - lon_min) / 64.0)
        y = lat_max - (np.arange(64, dtype=np.float64) + 0.5) * ((lat_max - lat_min) / 64.0)
        mask = xr.DataArray(
            np.zeros((64, 64), dtype=bool),
            dims=["y", "x"],
            coords={"x": x, "y": y},
        )
        mask.values[:, :32] = True
        mask = mask.rio.set_spatial_dims(x_dim="x", y_dim="y")
        mask = mask.rio.write_crs("EPSG:4326")
        mask = mask.rio.write_transform(from_bounds(lon_min, lat_min, lon_max, lat_max, 64, 64))

        resampled = assembler_mod._resample_cloud_mask(
            mask,
            target_template.shape,
            template=target_template,
        )
        expected_fraction = assembler_mod._resample_da(
            mask.astype(np.float32),
            target_template.shape,
            "area",
            template=target_template,
        )
        expected = expected_fraction.values > 0.0

        np.testing.assert_array_equal(resampled.values, expected)

    def test_cloud_mask_preserves_non_divisible_edge_clouds(self):
        mask = xr.DataArray(np.zeros((5, 5), dtype=bool), dims=["y", "x"])
        mask.values[-1, -1] = True

        resampled = assembler_mod._resample_cloud_mask(mask, (2, 2))

        expected = np.array([[False, False], [False, True]])
        np.testing.assert_array_equal(resampled.values, expected)

    def test_identity_same_res(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """When aerosol_resolution matches the scene pixel size, output ≈ input."""
        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=10.0,
        )
        # Should have same shape as input (64x64)
        assert sib.toa.shape[1:] == (64, 64)

    def test_downsampling_reduces_shape(
        self, large_obs_bundle, large_atmo, large_surface, mock_rt_model
    ):
        """Resampling to a coarser aerosol grid should produce a smaller grid."""
        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=500.0,
        )
        # 64 px @ 10m -> ~1.3 px @ 500m -> at least (1, 1)
        assert sib.toa.shape[1] <= 64
        assert sib.toa.shape[2] <= 64

    @pytest.mark.parametrize(
        ("feature_kind", "expected_cells"),
        [
            ("bright_point", [(0, 0)]),
            ("dark_point", [(0, 0)]),
            ("road", [(0, 0), (1, 0), (2, 0), (3, 0)]),
        ],
    )
    def test_sharp_transition_filter_excludes_solver_cells_for_native_targets(
        self,
        feature_kind,
        expected_cells,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        shape = (64, 64)
        toa = {
            name: xr.DataArray(
                np.full(shape, 0.18, dtype=np.float32),
                dims=["y", "x"],
            )
            for name in large_obs_bundle.toa.data_vars
        }
        if feature_kind == "bright_point":
            toa["B01"].values[8, 8] = 0.95
            toa["B02"].values[8, 8] = 0.95
        elif feature_kind == "dark_point":
            toa["B01"].values[8, 8] = 0.005
            toa["B02"].values[8, 8] = 0.005
        elif feature_kind == "road":
            toa["B01"].values[:, 8] = 0.95
            toa["B02"].values[:, 8] = 0.95
        else:
            raise AssertionError(f"Unhandled feature kind {feature_kind}")

        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset(toa),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        )

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=160.0,
            sharp_transition_filter=SharpTransitionFilterConfig(
                enabled=True,
                blur_kernel_pixels_native=31,
                residual_threshold_uint8=12,
                dilation_pixels=0,
                solver_cell_fraction_threshold=0.0,
            ),
        )

        exclusion_mask = sib.sharp_transition_mask
        assert exclusion_mask is not None
        assert exclusion_mask.shape == sib.cloud_mask.shape
        for y_idx, x_idx in expected_cells:
            assert bool(exclusion_mask.values[y_idx, x_idx])

    @pytest.mark.parametrize(
        ("feature_kind", "expected_cells"),
        [
            ("native_road", [(0, 0), (1, 0), (2, 0), (3, 0)]),
            ("native_point", [(0, 0)]),
        ],
    )
    def test_sharp_transition_filter_native_local_metrics_cover_realistic_subgrid_features(
        self,
        feature_kind,
        expected_cells,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        shape = (64, 64)
        toa = {
            name: xr.DataArray(
                np.full(shape, 0.18, dtype=np.float32),
                dims=["y", "x"],
            )
            for name in large_obs_bundle.toa.data_vars
        }
        if feature_kind == "native_road":
            toa["B01"].values[:, 8:12] = 0.30
            toa["B02"].values[:, 8:12] = 0.30
        elif feature_kind == "native_point":
            toa["B01"].values[8, 8] = 0.40
            toa["B02"].values[8, 8] = 0.40
        else:
            raise AssertionError(f"Unhandled feature kind {feature_kind}")

        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset(toa),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        )

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=160.0,
            sharp_transition_filter=SharpTransitionFilterConfig(
                enabled=True,
                blur_kernel_pixels_native=31,
                residual_threshold_uint8=12,
                dilation_pixels=0,
                solver_cell_fraction_threshold=0.0,
            ),
        )

        exclusion_mask = sib.sharp_transition_mask
        assert exclusion_mask is not None
        assert exclusion_mask.shape == sib.cloud_mask.shape
        for y_idx, x_idx in expected_cells:
            assert bool(exclusion_mask.values[y_idx, x_idx])

    def test_sharp_transition_filter_preserves_homogeneous_solver_cells(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        shape = (64, 64)
        toa = {
            name: xr.DataArray(
                np.full(shape, 0.18, dtype=np.float32),
                dims=["y", "x"],
            )
            for name in large_obs_bundle.toa.data_vars
        }
        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset(toa),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        )

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=160.0,
            sharp_transition_filter=SharpTransitionFilterConfig(
                enabled=True,
                blur_kernel_pixels_native=31,
                residual_threshold_uint8=12,
                dilation_pixels=0,
                solver_cell_fraction_threshold=0.0,
            ),
        )

        exclusion_mask = sib.sharp_transition_mask
        assert exclusion_mask is not None
        assert not bool(exclusion_mask.values.any())

    def test_sharp_transition_filter_marks_step_edge_boundary_cells(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        shape = (64, 64)
        toa = {
            name: xr.DataArray(
                np.full(shape, 0.18, dtype=np.float32),
                dims=["y", "x"],
            )
            for name in large_obs_bundle.toa.data_vars
        }
        for band in toa.values():
            band.values[:, 32:] = 0.30

        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset(toa),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
        )

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=160.0,
            sharp_transition_filter=SharpTransitionFilterConfig(
                enabled=True,
                blur_kernel_pixels_native=31,
                residual_threshold_uint8=12,
                dilation_pixels=0,
                solver_cell_fraction_threshold=0.0,
            ),
        )

        exclusion_mask = sib.sharp_transition_mask
        assert exclusion_mask is not None
        assert bool(exclusion_mask.values[:, 1].all())
        assert bool(exclusion_mask.values[:, 2].all())
        assert not bool(exclusion_mask.values[:, 0].any())
        assert not bool(exclusion_mask.values[:, 3].any())

    def test_sharp_transition_filter_cloud_buffer_blocks_halo_detections(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        shape = (64, 64)
        toa = {
            name: xr.DataArray(
                np.full(shape, 0.18, dtype=np.float32),
                dims=["y", "x"],
            )
            for name in large_obs_bundle.toa.data_vars
        }
        for band in toa.values():
            band.values[8, 8] = 0.95

        cloud_mask = xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"])
        cloud_mask.values[6:11, 6:11] = True
        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset(toa),
            cloud_mask=cloud_mask,
        )

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=160.0,
            sharp_transition_filter=SharpTransitionFilterConfig(
                enabled=True,
                blur_kernel_pixels_native=31,
                residual_threshold_uint8=12,
                dilation_pixels=0,
                cloud_buffer_pixels=2,
                solver_cell_fraction_threshold=0.0,
            ),
        )

        exclusion_mask = sib.sharp_transition_mask
        assert exclusion_mask is not None
        assert not bool(exclusion_mask.values.any())

    def test_toa_downsampling_uses_gdal_average_when_georeferenced(
        self,
        monkeypatch: pytest.MonkeyPatch,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ) -> None:
        import rioxarray  # noqa: F401
        from affine import Affine
        from rioxarray.raster_array import RasterArray

        x = 300000.0 + (np.arange(64, dtype=np.float64) + 0.5) * 10.0
        y = 5500640.0 - (np.arange(64, dtype=np.float64) + 0.5) * 10.0
        transform = Affine(10.0, 0.0, 300000.0, 0.0, -10.0, 5500640.0)

        def _georef(da: xr.DataArray) -> xr.DataArray:
            out = da.assign_coords({"x": x, "y": y})
            out = out.rio.set_spatial_dims(x_dim="x", y_dim="y")
            out = out.rio.write_crs("EPSG:32632")
            return out.rio.write_transform(transform)

        obs = dataclasses.replace(
            large_obs_bundle,
            toa=xr.Dataset(
                {name: _georef(data) for name, data in large_obs_bundle.toa.data_vars.items()}
            ),
        )

        calls: list[RasterioResampling] = []
        original_reproject_match = RasterArray.reproject_match

        def _fake_reproject_match(self, target, *, resampling, **kwargs):  # type: ignore[no-untyped-def]
            calls.append(resampling)
            return original_reproject_match(self, target, resampling=resampling, **kwargs)

        monkeypatch.setattr(RasterArray, "reproject_match", _fake_reproject_match)

        sib = assemble_grids(
            obs,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=320.0,
        )

        assert sib.toa.shape[1:] == (2, 2)
        assert calls == [RasterioResampling.average, RasterioResampling.average]

    def test_georeferenced_bilinear_and_nearest_dispatch_use_gdal(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        import rioxarray  # noqa: F401
        from affine import Affine
        from rioxarray.raster_array import RasterArray

        x = 300000.0 + (np.arange(4, dtype=np.float64) + 0.5) * 10.0
        y = 5500040.0 - (np.arange(4, dtype=np.float64) + 0.5) * 10.0
        transform = Affine(10.0, 0.0, 300000.0, 0.0, -10.0, 5500040.0)

        source = xr.DataArray(
            np.arange(16, dtype=np.float32).reshape(4, 4),
            dims=["y", "x"],
            coords={"x": x, "y": y},
        )
        source = source.rio.set_spatial_dims(x_dim="x", y_dim="y")
        source = source.rio.write_crs("EPSG:32632")
        source = source.rio.write_transform(transform)

        template = xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=["y", "x"],
            coords={
                "x": 300000.0 + (np.arange(2, dtype=np.float64) + 0.5) * 20.0,
                "y": 5500040.0 - (np.arange(2, dtype=np.float64) + 0.5) * 20.0,
            },
        )
        template = template.rio.set_spatial_dims(x_dim="x", y_dim="y")
        template = template.rio.write_crs("EPSG:32632")
        template = template.rio.write_transform(Affine(20.0, 0.0, 300000.0, 0.0, -20.0, 5500040.0))

        calls: list[RasterioResampling] = []
        original_reproject_match = RasterArray.reproject_match

        def _fake_reproject_match(self, target, *, resampling, **kwargs):  # type: ignore[no-untyped-def]
            calls.append(resampling)
            return original_reproject_match(self, target, resampling=resampling, **kwargs)

        monkeypatch.setattr(RasterArray, "reproject_match", _fake_reproject_match)

        bilinear = assembler_mod._resample_da(source, (2, 2), method="bilinear", template=template)
        nearest = assembler_mod._resample_da(source, (2, 2), method="nearest", template=template)

        assert bilinear.shape == (2, 2)
        assert nearest.shape == (2, 2)
        assert calls == [RasterioResampling.bilinear, RasterioResampling.nearest]

    def test_resample_da_falls_back_when_gdal_reports_missing_georeferencing(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        import rioxarray  # noqa: F401
        from affine import Affine
        from rasterio.errors import NotGeoreferencedWarning
        from rioxarray.raster_array import RasterArray

        source = xr.DataArray(
            np.arange(16, dtype=np.float32).reshape(4, 4),
            dims=["y", "x"],
            coords={
                "x": 300000.0 + (np.arange(4, dtype=np.float64) + 0.5) * 10.0,
                "y": 5500040.0 - (np.arange(4, dtype=np.float64) + 0.5) * 10.0,
            },
        )
        source = source.rio.set_spatial_dims(x_dim="x", y_dim="y")
        source = source.rio.write_crs("EPSG:32632")

        template = xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=["y", "x"],
            coords={
                "x": 300000.0 + (np.arange(2, dtype=np.float64) + 0.5) * 20.0,
                "y": 5500040.0 - (np.arange(2, dtype=np.float64) + 0.5) * 20.0,
            },
        )
        template = template.rio.set_spatial_dims(x_dim="x", y_dim="y")
        template = template.rio.write_crs("EPSG:32632")
        template = template.rio.write_transform(Affine(20.0, 0.0, 300000.0, 0.0, -20.0, 5500040.0))

        def _warn_reproject_match(self, target, *, resampling, **kwargs):  # type: ignore[no-untyped-def]
            del self, target, resampling, kwargs
            warnings.warn(
                "Dataset has no geotransform, gcps, or rpcs. The identity matrix will be returned.",
                NotGeoreferencedWarning,
                stacklevel=1,
            )
            raise AssertionError("GDAL resampling should fall back after the warning.")

        monkeypatch.setattr(RasterArray, "reproject_match", _warn_reproject_match)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = assembler_mod._resample_da(source, (2, 2), method="area", template=template)

        assert out.shape == (2, 2)
        assert not any(isinstance(item.message, NotGeoreferencedWarning) for item in caught)

    def test_geographic_source_dims_dispatch_use_gdal(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        import rioxarray  # noqa: F401
        from affine import Affine
        from rioxarray.raster_array import RasterArray

        source = xr.DataArray(
            np.arange(16, dtype=np.float32).reshape(4, 4),
            dims=["latitude", "longitude"],
            coords={
                "longitude": 8.0 + (np.arange(4, dtype=np.float64) + 0.5) * 0.1,
                "latitude": 50.4 - (np.arange(4, dtype=np.float64) + 0.5) * 0.1,
            },
        )
        source = source.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
        source = source.rio.write_crs("EPSG:4326")
        source = source.rio.write_transform(Affine(0.1, 0.0, 8.0, 0.0, -0.1, 50.4))

        template = xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=["y", "x"],
            coords={
                "x": 8.0 + (np.arange(2, dtype=np.float64) + 0.5) * 0.2,
                "y": 50.4 - (np.arange(2, dtype=np.float64) + 0.5) * 0.2,
            },
        )
        template = template.rio.set_spatial_dims(x_dim="x", y_dim="y")
        template = template.rio.write_crs("EPSG:4326")
        template = template.rio.write_transform(Affine(0.2, 0.0, 8.0, 0.0, -0.2, 50.4))

        calls: list[RasterioResampling] = []
        original_reproject_match = RasterArray.reproject_match

        def _fake_reproject_match(self, target, *, resampling, **kwargs):  # type: ignore[no-untyped-def]
            calls.append(resampling)
            return original_reproject_match(self, target, resampling=resampling, **kwargs)

        monkeypatch.setattr(RasterArray, "reproject_match", _fake_reproject_match)

        bilinear = assembler_mod._resample_da(source, (2, 2), method="bilinear", template=template)
        nearest = assembler_mod._resample_da(source, (2, 2), method="nearest", template=template)

        assert bilinear.dims == ("y", "x")
        assert nearest.dims == ("y", "x")
        assert calls == [RasterioResampling.bilinear, RasterioResampling.nearest]

    def test_aerosol_resolution_controls_grid_and_georeference(
        self,
        large_obs_bundle,
        large_atmo,
        large_surface,
        mock_rt_model,
    ):
        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            large_surface,
            mock_rt_model,
            aerosol_resolution_m=320.0,
        )

        assert sib.toa.shape[1:] == (2, 2)
        assert sib.atmo_prior.aot.shape == (2, 2)
        assert sib.atmo_prior.aot.rio.crs is not None
        assert sib.atmo_prior.aot.rio.crs.to_string() == "EPSG:32632"
        assert sib.atmo_prior.aot.coords["x"].values.tolist() == pytest.approx([300160.0, 300480.0])
        assert sib.atmo_prior.aot.coords["y"].values.tolist() == pytest.approx(
            [5500480.0, 5500160.0]
        )

    def test_passes_validation(self, large_obs_bundle, large_atmo, large_surface, mock_rt_model):
        """Output should pass validate_solver_input_bundle()."""
        sib = assemble_grids(large_obs_bundle, large_atmo, large_surface, mock_rt_model)
        validate_solver_input_bundle(sib)  # should not raise

    def test_surface_prior_with_band_dimension(
        self, large_obs_bundle, large_atmo, large_surface, mock_rt_model
    ):
        """Assembler should handle banded SurfacePrior arrays from real BRDF providers."""
        shape = large_surface.boa.shape
        banded_boa = xr.DataArray(
            np.stack(
                [
                    np.full(shape, 0.11, dtype=np.float32),
                    np.full(shape, 0.13, dtype=np.float32),
                ],
                axis=0,
            ),
            dims=["band", "y", "x"],
            coords={"band": [1, 2]},
        )
        banded_unc = xr.full_like(banded_boa, 0.02)
        banded_mask = xr.DataArray(
            np.ones((2, *shape), dtype=bool),
            dims=["band", "y", "x"],
            coords={"band": [1, 2]},
        )
        banded_surface = SurfacePrior(
            boa=banded_boa,
            boa_unc=banded_unc,
            kernels=large_surface.kernels,
            mask=banded_mask,
        )

        sib = assemble_grids(large_obs_bundle, large_atmo, banded_surface, mock_rt_model)
        assert sib.surface_prior.boa.ndim == 3
        assert sib.surface_prior.boa.shape[0] == 2
        assert sib.surface_prior.mask.ndim == 2

    def test_surface_prior_at_solver_resolution_is_not_degraded(
        self,
        large_obs_bundle,
        large_atmo,
        mock_rt_model,
    ):
        """When surface prior is already at aerosol_resolution, data passes through unchanged.

        M3 queries the surface prior at aerosol_resolution using the same
        _build_target_template(bounds, crs, resolution) as M4.  When grids
        match, _resample_da detects the shared template and returns data as-is.
        This test verifies no numeric degradation from a redundant resample.
        """
        from siac.algorithms.grid.assembler import _build_target_template

        resolution_m = 320.0
        template = _build_target_template(
            large_obs_bundle.bounds,
            large_obs_bundle.crs,
            resolution_m,
        )
        target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))

        rng = np.random.RandomState(42)
        boa_vals = rng.uniform(0.05, 0.3, target_shape).astype(np.float32)
        unc_vals = rng.uniform(0.01, 0.05, target_shape).astype(np.float32)
        mask_vals = np.ones(target_shape, dtype=bool)
        mask_vals[0, 0] = False  # one invalid pixel

        from siac.geo.resample import copy_spatial_metadata_like

        boa = copy_spatial_metadata_like(
            xr.DataArray(boa_vals, dims=("y", "x"), coords=template.coords),
            template,
        )
        unc = copy_spatial_metadata_like(
            xr.DataArray(unc_vals, dims=("y", "x"), coords=template.coords),
            template,
        )
        mask = copy_spatial_metadata_like(
            xr.DataArray(mask_vals, dims=("y", "x"), coords=template.coords),
            template,
        )

        prior = SurfacePrior(boa=boa, boa_unc=unc, kernels=None, mask=mask)

        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            prior,
            mock_rt_model,
            aerosol_resolution_m=resolution_m,
        )

        np.testing.assert_array_equal(sib.surface_prior.boa.values, boa_vals)
        np.testing.assert_array_equal(sib.surface_prior.boa_unc.values, unc_vals)
        np.testing.assert_array_equal(sib.surface_prior.mask.values, mask_vals)

    def test_legacy_aux_resolution_override_resamples_surface_prior(
        self,
        large_obs_bundle,
        large_atmo,
        mock_rt_model,
    ):
        """When aux_resolution != 500 and aerosol_resolution == 120, M4 uses aux_resolution.

        This is the one scenario where M3 (querying at aerosol_resolution=120)
        and M4 (gridding at aux_resolution) can differ, producing a real resample.
        """
        from siac.algorithms.grid.assembler import _build_target_template

        # Surface prior generated at 120m (M3 default)
        m3_template = _build_target_template(
            large_obs_bundle.bounds,
            large_obs_bundle.crs,
            120.0,
        )
        m3_shape = (int(m3_template.sizes["y"]), int(m3_template.sizes["x"]))

        from siac.geo.resample import copy_spatial_metadata_like

        boa = copy_spatial_metadata_like(
            xr.DataArray(
                np.full(m3_shape, 0.15, dtype=np.float32),
                dims=("y", "x"),
                coords=m3_template.coords,
            ),
            m3_template,
        )
        unc = copy_spatial_metadata_like(
            xr.DataArray(
                np.full(m3_shape, 0.02, dtype=np.float32),
                dims=("y", "x"),
                coords=m3_template.coords,
            ),
            m3_template,
        )
        mask = copy_spatial_metadata_like(
            xr.DataArray(
                np.ones(m3_shape, dtype=bool),
                dims=("y", "x"),
                coords=m3_template.coords,
            ),
            m3_template,
        )
        prior = SurfacePrior(boa=boa, boa_unc=unc, kernels=None, mask=mask)

        # Trigger legacy override: aerosol_resolution=120 + aux_resolution=320
        sib = assemble_grids(
            large_obs_bundle,
            large_atmo,
            prior,
            mock_rt_model,
            aux_resolution_m=320.0,
            aerosol_resolution_m=120.0,
        )

        # Grid should be at 320m (legacy override), not 120m
        m4_template = _build_target_template(
            large_obs_bundle.bounds,
            large_obs_bundle.crs,
            320.0,
        )
        expected_shape = (int(m4_template.sizes["y"]), int(m4_template.sizes["x"]))
        assert sib.toa.shape[1:] == expected_shape
        assert sib.surface_prior.boa.shape == expected_shape
        assert sib.aerosol_resolution_m == pytest.approx(320.0)
