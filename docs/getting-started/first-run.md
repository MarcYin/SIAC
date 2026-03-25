# First Run

This page shows the shortest useful SIAC runs: one through the CLI and one through Python. Both focus on the current primary workflow, Sentinel-2.

## CLI example

```bash
siac process-s2 T31UDQ_20240101 --output-path ./outputs/example
```

This command:

1. resolves the query into a local SAFE path or remote backend lookup
2. loads configuration defaults or your configured local defaults
3. runs the M1-M6 pipeline
4. writes outputs under `./outputs/example`

If you want to constrain the run to a smaller area:

```bash
siac process-s2 S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433 \
  --output-path ./outputs/example \
  --aoi-bbox 300000 1900000 400000 2000000 \
  --aoi-crs EPSG:32650
```

The AOI above covers the full `T50QLD` tile used by the current Sentinel-2C trial scene.

## Python example

```python
from siac import SIAC
from siac.config import SIACConfig

config = SIACConfig(sensor="s2")
siac = SIAC(config)
result = siac.process("/path/to/S2_SAFE/")

print(result.boa)
print(float(result.aot.mean()))
```

This path is useful when you want to inspect the typed result directly instead of relying only on files written to disk.

## What you should expect

Depending on output settings, SIAC writes:

- BOA reflectance rasters
- optional BOA uncertainty rasters
- AOT and TCWV auxiliary outputs
- a cloud mask
- optional quicklook imagery

See [Outputs](../user-guide/outputs.md) for the full layout.

## Before your first remote run

If your input is not already a local SAFE path, configure the credentials needed by your chosen backend and providers. The defaults and environment overlay behavior are explained in [Configuration Basics](../user-guide/configuration-basics.md).
