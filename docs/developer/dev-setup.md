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

## Common Commands

| Goal | Command |
| --- | --- |
| install editable package | `pixi run bootstrap` |
| run full test suite | `pixi run test` |
| run fast test slice | `pixi run test-fast` |
| run tests without coverage | `pixi run test-no-cov` |
| run lint | `pixi run lint` |
| format code | `pixi run format` |
| type-check main scoped targets | `pixi run typecheck-scoped` |
| build Rust extension | `pixi run build-rust` |

## Expected Local Workflow

1. `pixi install`
2. `pixi run bootstrap`
3. make changes
4. `pixi run lint`
5. `pixi run typecheck-scoped`
6. `pixi run test-fast` or a narrower targeted test set
7. `pixi run build-rust` if your change depends on native code

## Test Layout

| Test area | Purpose |
| --- | --- |
| `tests/unit/` | package-level behavior, config, adapters, algorithms, workflow boundaries |
| `tests/integration/` | end-to-end slices and backend-oriented workflow coverage |
| `tests/regression/` | golden-output style protection when present |
| `tests/benchmarks/` | performance checks |

The CI pipeline currently runs unit tests, an integration slice, Rust build
smoke, and wheel install smoke checks. If you change public entry points,
configuration behavior, or runtime orchestration, update tests in the matching
layer.

## Docs Build Commands

The documentation site is intended to compile in CI and GitHub Pages using the
same local command. Prefer the project task when it exists:

```bash
pixi run docs-build
```

If you need a direct fallback, use:

```bash
sphinx-build -M html docs docs/_build -W
```

Use the fallback when debugging the docs toolchain or before a dedicated task is
available locally.

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
