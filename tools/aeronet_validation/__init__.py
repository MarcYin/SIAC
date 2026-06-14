"""AERONET validation experiment for SIAC surface-prior approaches.

Pipeline stages (run via ``python -m tools.aeronet_validation.cli``):

1. ``fetch-aeronet``  - download global AERONET V3 AOD measurements.
2. ``matchup``        - find Sentinel-2 acquisitions coincident with AERONET.
3. ``run``            - SIAC aerosol retrieval per matchup x surface-prior approach.
4. ``compare``        - score retrieved AOT against AERONET AOD at 550 nm.
"""
