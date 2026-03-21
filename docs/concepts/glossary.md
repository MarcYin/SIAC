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
