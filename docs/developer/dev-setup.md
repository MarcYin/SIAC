# Developer Setup

SIAC development is centered on the Pixi workspace in the repository root. The
workspace installs Python, Rust, geospatial dependencies, test tooling, and the
developer tasks used by CI.

## Prerequisites

- Git
- Pixi
- a supported Python version from the workspace
- Rust toolchain through the Pixi environment

The repository already declares the default development environment in
`pixi.toml`, so the normal path is to use Pixi rather than manage a separate
system Python and Rust installation by hand.

## Bootstrap

Install the environment and editable package:

```bash
pixi install
pixi run bootstrap
```

Build the native extension when you need the Rust-backed code paths:

```bash
pixi run build-rust
```

Run that build before `pixi run test`, `pixi run test-no-cov`, or
`pixi run coverage` in a fresh environment; those suites import `siac._rust`
during collection.

## Common Commands

| Goal | Command |
| --- | --- |
| install editable package | `pixi run bootstrap` |
| run full test suite | `pixi run test` |
| run fast test slice | `pixi run test-fast` |
| run tests without coverage | `pixi run test-no-cov` |
| run coverage gate | `pixi run coverage` |
| run lint | `pixi run lint` |
| format code | `pixi run format` |
| type-check main scoped targets | `pixi run typecheck-scoped` |
| build Rust extension | `pixi run build-rust` |

## Expected Local Workflow

1. `pixi install`
2. `pixi run bootstrap`
3. `pixi run build-rust`
4. make changes
5. `pixi run lint`
6. `pixi run typecheck-scoped`
7. `pixi run test-fast` or a narrower targeted test set
8. `pixi run test` before landing a broader change
9. `pixi run coverage` when you need to verify the enforced threshold
10. rerun `pixi run build-rust` after changing the native extension

## Test Layout

| Test area | Purpose |
| --- | --- |
| `tests/unit/` | package-level behavior, config, adapters, algorithms, workflow boundaries |
| `tests/integration/` | end-to-end slices and backend-oriented workflow coverage |
| `tests/regression/` | golden-output style protection when present |
| `tests/benchmarks/` | performance checks |

The CI pipeline currently runs unit tests, an integration slice, Rust build
smoke, wheel install smoke checks, and the coverage gate. If you change public
entry points, configuration behavior, or runtime orchestration, update tests in
the matching layer.

The repository also has a dedicated `Native 6S Smoke` GitHub Actions workflow.
That workflow validates the compiled `sixs` backend on Linux with the `rt6s`
Pixi environment, which pins Python 3.11 and `setuptools < 60` so the workflow
can force the compatible F2PY `distutils` backend for this legacy extension. It
downloads a fresh upstream 6SV2.1 source tree, builds the native module, and
runs a minimal array-backed correction case. When the smoke lane fails it now
persists the build diagnostics under `tmp/rt6s_ci_smoke/diagnostics/`, prints
the relevant log tails in the workflow output, and uploads the full
`tmp/rt6s_ci_smoke` tree as a workflow artifact so Linux-only builder failures
can be inspected after the run. Changes under
`python/siac/algorithms/rt/direct/`, `python/siac/config/schema.py`,
`python/siac/sixs_outputs.py`, and the native 6S tools/workflow files should be
validated against that workflow before landing.

## Docs Build Commands

The documentation site is intended to compile in CI and GitHub Pages using the
same local command. Prefer the project task when it exists:

```bash
pixi run docs-build
```

If you need to run the underlying command directly, use:

```bash
mkdocs build --strict
```

Use the direct command when debugging the docs toolchain or when you need to
confirm the same build outside a task wrapper.

## Working With The Layered Architecture

Use the package boundaries intentionally:

- add or change public user-facing behavior in `api/`
- add configuration schema and resolution changes in `config/`
- add assembly or registry behavior in `app/`
- add stage orchestration in `workflows/`
- add remote/backend integrations in `adapters/`
- add scientific or numerical behavior in `algorithms/`
- add shared runtime payloads in `runtime/`

Avoid collapsing responsibilities into a single layer. The codebase is easier to
extend when config, orchestration, data access, and numerical logic remain
separate.
