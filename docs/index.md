# SIAC Documentation

Sensor-Invariant Atmospheric Correction (SIAC) is a modular atmospheric correction framework for medium-resolution optical satellite data. The current codebase is organized around a small public API, a typed configuration system, a staged runtime pipeline, and pluggable adapters and algorithms for data access, priors, radiative transfer, and output writing.

```{mermaid}
flowchart LR
    User["Users and scripts"] --> API["API / CLI"]
    API --> Config["Config load + resolve"]
    Config --> App["Request coercion + assembly"]
    App --> Workflow["Workflows + pipeline"]
    Workflow --> Adapters["Adapters"]
    Workflow --> Algorithms["Algorithms"]
    Workflow --> Runtime["Runtime payloads"]
    Workflow --> Storage["Storage / outputs"]
    Algorithms --> Rust["Rust acceleration"]
```

## Choose a path

- **First-time user**: start with installation, first run, and running SIAC.
- **Scientific user**: start with the overview, theory, priors/mapping, and paper companion.
- **Operator/engineer**: start with configuration basics, outputs, troubleshooting, and developer setup.
- **Developer**: start with the package map, execution flow, data flow, and extension guide.

## Documentation map

```{toctree}
:maxdepth: 2
:caption: Concepts

concepts/overview
concepts/glossary
```

```{toctree}
:maxdepth: 2
:caption: Getting Started

getting-started/installation
getting-started/first-run
```

```{toctree}
:maxdepth: 2
:caption: User Guide

user-guide/configuration-basics
user-guide/running-siac
user-guide/outputs
user-guide/troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: Science

science/theory
science/priors-and-mapping
science/paper-companion
science/validation-and-limitations
```

```{toctree}
:maxdepth: 2
:caption: Reference

reference/configuration
reference/public-api
```

```{toctree}
:maxdepth: 2
:caption: Architecture

architecture/package-map
architecture/execution-flow
architecture/data-flow
architecture/rust-integration
```

```{toctree}
:maxdepth: 2
:caption: Developer

developer/dev-setup
developer/extending-siac
```
