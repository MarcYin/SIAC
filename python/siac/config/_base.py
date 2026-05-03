"""Shared base model for SIAC configuration sections."""

from __future__ import annotations

from pydantic import BaseModel


class SIACBaseModel(BaseModel):
    model_config = {
        "extra": "forbid",
        "validate_default": True,
    }
