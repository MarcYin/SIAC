"""Cloud and cloud-shadow detection utilities."""

from siac.cloud.mapping import apply_class_mapping, validate_class_mapping
from siac.cloud.mask import build_cloud_classes, classes_to_bool_mask
from siac.cloud.providers import OmniCloudMaskProvider

__all__ = [
    "OmniCloudMaskProvider",
    "apply_class_mapping",
    "validate_class_mapping",
    "build_cloud_classes",
    "classes_to_bool_mask",
]
