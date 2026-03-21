# Outputs

SIAC returns a typed `CorrectionResult` in Python and can also write outputs to disk through the configured output writer.

## Main artifacts

| Artifact | Meaning |
| --- | --- |
| `boa` | BOA reflectance dataset |
| `boa_unc` | Optional BOA uncertainty dataset |
| `aot` | Retrieved aerosol optical thickness |
| `tcwv` | Retrieved total column water vapour |
| `cloud_mask` | Cloud mask used by the workflow |

## Typical directory layout

When output writing is enabled, the layout usually looks like:

```text
output-dir/
├── boa/
├── boa_unc/
├── auxiliary/
└── quicklook.tif
```

The exact artifacts depend on `output.defaults`.

## Output formats

SIAC supports:

- GeoTIFF
- COG
- NetCDF
- Zarr

## Key output settings

| Setting | Effect |
| --- | --- |
| `output.defaults.format` | Selects raster or container format |
| `output.defaults.include_uncertainty` | Enables BOA uncertainty output |
| `output.defaults.include_auxiliary` | Writes AOT, TCWV, and cloud mask |
| `output.defaults.include_rgb` | Writes RGB quicklook when compatible bands are present |
| `output.defaults.boa_dtype` | Controls BOA encoding dtype |
| `output.defaults.compression` | Controls output compression |

## Interpreting the result in Python

```python
result.boa
result.boa_unc
result.aot
result.tcwv
result.cloud_mask
```

The result object is the cleanest way to inspect outputs inside notebooks, scripts, or downstream workflows before or instead of writing them to disk.
