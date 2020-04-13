.. pamtra2 documentation master file, created by
   sphinx-quickstart on Sun Jun 17 16:56:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pamtra2 - Passive and Active Microwave TRANsfer 
##########################################################


Why pamtra2?
************

The original PAMTRA (still available here https://github.com/igmk/pamtra) was written in Fortran with an optional Python wrapper. We decided to start developing pamtra2, because we want to

* have easier to maintain code
* make it easier to extend with the potential to add eventually additional forward operators (e.g. lidar) 
* make installation less painful
* have the model core in pure Python
* use Fortran, C, Cython or Numba for expensive parts


Included Libraries
******************

Besides the core python package `pamtra2` , additonal are libraries provided in `/libs`. They are installed as part of pamtra2 in `pamtra2.libs`, but can be also installed separately with the provided `setup.py` files:

* meteo_si - Meteo SI is a collection of functions for atmospheric sciences following the SI-convention (unless stated otherwise)
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
   examples




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
