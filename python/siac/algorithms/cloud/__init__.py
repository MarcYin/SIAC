"""Cloud and cloud-shadow detection utilities."""

from siac.algorithms.cloud.mapping import apply_class_mapping as apply_class_mapping
from siac.algorithms.cloud.mapping import validate_class_mapping as validate_class_mapping
from siac.algorithms.cloud.mask import build_cloud_classes as build_cloud_classes
from siac.algorithms.cloud.mask import classes_to_bool_mask as classes_to_bool_mask
from siac.algorithms.cloud.providers import OmniCloudMaskProvider as OmniCloudMaskProvider
