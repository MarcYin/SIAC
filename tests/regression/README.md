# Regression tests

This package is reserved for regression tests that compare numerical
outputs against stored golden fixtures to detect behavior drift.

## Status

The package is currently **empty** — no regression tests live here yet.
It is kept in the tree so the `regression` marker (declared in
`pyproject.toml`'s `[tool.pytest.ini_options]`) has a documented home and
so future contributors have an obvious place to add golden-output tests.

## Conventions

- Every test in this directory MUST be decorated with
  `@pytest.mark.regression` (or set
  `pytestmark = pytest.mark.regression` at the module level).
- The `regression` marker is **deselected from the fast suite**. Run
  these tests explicitly with `pytest -m regression`.
- Golden fixtures should live under a `golden/` subdirectory and be
  versioned alongside the test that consumes them.

See `docs/TEST_PLAN.md` (Layer 8) and `docs/developer/dev-setup.md` for
the full conventions and CI integration.
