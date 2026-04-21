# Glossary

## Core remote-sensing terms

**TOA**
: Top-of-atmosphere reflectance or radiance measured by the sensor before atmospheric correction.

**BOA**
: Bottom-of-atmosphere reflectance, the corrected surface reflectance target product.

**AOT**
: Aerosol optical thickness, a key atmospheric parameter retrieved or constrained during correction.

**TCWV**
: Total column water vapour.

**TCO3**
: Total column ozone.

**BRDF**
: Bidirectional reflectance distribution function. In SIAC, BRDF products help provide a surface prior.

**ARD**
: Analysis Ready Data. SIAC aims to produce correction outputs suitable for downstream analysis.

## SIAC-specific workflow terms

**AOI**
: Area of interest. In SIAC this may come from bounds, GeoJSON, or raster-backed geometry.

**SAFE**
: Sentinel-2 packaging format used for local scene input.

**Prior**
: A probabilistic or constraint input used before solving. SIAC uses atmospheric and surface priors.

**Spectral mapping**
: Mapping reflectance information between differing sensor spectral response functions.

**Spatial mapping / ePSF**
: Scale matching between fine-resolution observations and coarser prior products using an effective point spread function.

**Emulator**
: A fast learned approximation to radiative transfer model outputs.

**LUT**
: Lookup table for radiative transfer coefficients.

**ExecutionPlan**
: The fully assembled runtime plan that binds configuration, auth, inputs, and concrete callables for a run.

**ObservationBundle**
: The runtime payload produced by preprocessing, containing TOA data, geometry, cloud mask, metadata, and sensor description.

**CorrectionResult**
: The final typed result containing BOA outputs, auxiliary fields, and diagnostics.

## Solver and algorithm terms

**Multi-grid solver**
: The default solver strategy in SIAC. Retrieval starts on a coarse grid and progressively refines to the target aerosol resolution. Each level is solved with L-BFGS-B optimization. See `algorithms.solver.multigrid`.

**DCT regularization**
: Discrete Cosine Transform–based smoothness penalty applied to the atmospheric state during retrieval. Encourages spatially smooth AOT and TCWV fields while preserving large-scale gradients. See the smoothness term in `algorithms.solver.cost`.

**Band weight power (alpha)**
: Exponent controlling wavelength-dependent weighting in the observation cost. Configured as `algorithms.solver.alpha` (default `-1.6`). Higher absolute values give more weight to shorter wavelengths. Propagated through `MultiGridConfig.band_weight_power`.

**RAA**
: Relative Azimuth Angle, computed as `|VAA − SAA| mod 2π`. Used by RT models and BRDF kernels to parameterize the scattering geometry.

**RTModelBackend**
: The structural protocol (`domain.protocols.RTModelBackend`) that all radiative transfer backends must satisfy. Declares `simulate_toa`, `compute_coefficients`, and `supported_parameters`. Backends include the emulator, LUT, and native 6S adapters.

**SolverInputBundle**
: The typed payload consumed by the solver, containing resampled TOA, geometry, cloud mask, priors, RT model, and band metadata assembled onto the solver grid.
