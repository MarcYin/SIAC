.. SIAC documentation master file, created by
   sphinx-quickstart on Thu Nov  8 18:01:40 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MULTIPLY - Sensor Invariant Atmospheric Correction (SIAC)
=========================================================

This atmospheric correction method uses MODIS MCD43 BRDF product to get
a coarse resolution simulation of earth surface. A model based on MODIS
PSF is built to deal with the scale differences between MODIS and other
sensors, and linear spectral mapping is used to map between different
sensors spectrally. We uses the ECMWF
`CAMS <http://apps.ecmwf.int/datasets/data/cams-nrealtime/levtype=sfc/>`__
prediction as a prior for the atmospheric states, coupling with 6S model
to solve for the atmospheric parameters, then the solved atmospheric
parameters are used to correct the TOA reflectances. The whole system is
built under Bayesian theory and the uncertainty is propagated through
the whole system. Since we do not rely on specific bands' relationship
to estimate the atmospheric states, but instead a more generic and
consistent way of inversion those parameters. The code can be downloaded
from `SIAC <https://github.com/MarcYin/Atmospheric_correction>`__ github
directly and futrher updates will make it more independent and can be
installed on different machines.



Development of this code has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 687320, under project H2020 MULTIPLY.
f this code has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 687320, under project `H2020 MULTIPLY <http://www.multiply-h2020.eu/>`_.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   01_Introduction.rst
   02_Installation.rst
   03_Quickstart.rst
   04_Data_access_and_adaption_for_Landsat5.rst
   05_Annex_technical.rst

.. toctree::
    :maxdepth: 2
    :caption: anything else

    Authors <authors>
    Changelog <changes>

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
