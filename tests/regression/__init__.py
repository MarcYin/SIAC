"""Regression tests comparing against golden outputs.

Tests here must be decorated with ``@pytest.mark.regression`` and compare
numerical outputs against stored golden fixtures to detect behavior drift.
See ``docs/TEST_PLAN.md`` (Layer 8) for the conventions and
``docs/developer/dev-setup.md`` for how the suite integrates with CI.

Run with ``pytest -m regression`` (deselected from the fast suite).
"""
