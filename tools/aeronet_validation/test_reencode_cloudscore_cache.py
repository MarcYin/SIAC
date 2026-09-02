"""Tests for the uint8 clear-score archive encoding."""

from __future__ import annotations

import numpy as np
import pytest
from tools.aeronet_validation.reencode_cloudscore_cache import decode, encode
from tools.aeronet_validation.build_cloudscore_index_local import (
    CLEAR_SCORE_NODATA,
    CLEAR_SCORE_SCALE,
)


def test_a_perfectly_clear_pixel_is_not_confused_with_nodata() -> None:
    """The reason the scale is 254 rather than 255.

    Scaling by 255 puts cs=1.0 on the nodata sentinel, which silently converts
    the clearest pixels in a scene into gaps -- and gaps cannot win a pixel, so
    the mosaic would quietly avoid exactly the best observations.
    """

    encoded = encode(np.array([[1.0, np.nan]]))
    assert encoded[0, 0] != CLEAR_SCORE_NODATA
    assert encoded[0, 1] == CLEAR_SCORE_NODATA
    assert np.isnan(decode(encoded)[0, 1])
    assert decode(encoded)[0, 0] == pytest.approx(1.0)


def test_round_trip_stays_within_half_a_quantisation_step() -> None:
    values = np.linspace(0.0, 1.0, 501).reshape(1, -1)
    error = np.abs(decode(encode(values)) - values)
    assert error.max() <= 0.5 / CLEAR_SCORE_SCALE + 1e-12


def test_out_of_range_scores_are_clipped_not_wrapped() -> None:
    # uint8 overflow would turn a slightly-over-1.0 score into a near-zero one.
    encoded = encode(np.array([[-0.2, 1.4]]))
    assert decode(encoded)[0, 0] == pytest.approx(0.0)
    assert decode(encoded)[0, 1] == pytest.approx(1.0)


def test_nodata_survives_a_round_trip_as_nan_not_a_score() -> None:
    values = np.array([[np.nan, 0.5, np.inf]])
    decoded = decode(encode(values))
    assert np.isnan(decoded[0, 0])
    assert np.isnan(decoded[0, 2])
    assert decoded[0, 1] == pytest.approx(0.5, abs=1 / CLEAR_SCORE_SCALE)
