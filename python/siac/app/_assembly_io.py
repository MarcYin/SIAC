"""Output and grid-assembler assembly helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from siac.adapters.output import ConfiguredOutputWriter

if TYPE_CHECKING:
    from siac.workflows.pipeline import GridAssemblerFn


def resolve_output_writer(config: Any) -> ConfiguredOutputWriter:
    return ConfiguredOutputWriter(config.output.defaults)


def resolve_grid_assembler() -> GridAssemblerFn:
    from siac.algorithms.grid.assembler import assemble_grids

    return assemble_grids
