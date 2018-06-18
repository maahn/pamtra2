.. pamtra2 documentation master file, created by
   sphinx-quickstart on Sun Jun 17 16:56:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pamtra2's documentation!
###################################


Included Modules
****************

* pamtra2 - The core package.

All other packages are provided in `/libs`. They are installed as part of pamtra2 in `pamtra2.libs`, but can be also installed separately with the provided `setup.py` files:

* pamgasabs - Estimate gaseous absorption
* pyPamtraRadarSimulator - Simulate a radar Doppler spectrum
* pyPamtraRadarMoments - Estimate the moments of the radar Doppler spectrum
* refractiveIndex - general purpose refractive index of ice, water and Effective Medium Approximations
* singleScattering - general purpose electromagnetic single-scattering module. Includes several scattering methods




.. toctree::
   :numbered:
   :glob:

   install
   tutorial
   pamtra2
   libs
   examples



.. automodule:: pamtra2
   :members:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
