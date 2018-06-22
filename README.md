[![Documentation Status](//readthedocs.org/projects/pamtra2/badge/?version=latest)](https://pamtra2.readthedocs.io/en/latest/?badge=latest)


Pamtra2
#######

EARLY ALPHA VERSION! NOT FEATURE COMPLETE!

Included Modules
================

* pamtra2 - The core package.

All other packages are provided in `/libs`. They are installed as part of pamtra2 in `pamtra2.libs`, but can be also installed separately with the provided `setup.py` files:

* meteo_si - Meteo SI is a collection of functions for atmospheric sciences following the 
SI-convention (unless stated otherwise)
* pamgasabs - Estimate gaseous absorption
* pyPamtraRadarSimulator - Simulate a radar Doppler spectrum
* pyPamtraRadarMoments - Estimate the moments of the radar Doppler spectrum
* refractiveIndex - general purpose refractive index of ice, water and Effective Medium Approximations
* singleScattering - general purpose electromagnetic single-scattering module. Includes several scattering methods

Installation
============

First install all required python packages ::

    numpy
    scipy
    xarray
    dask
    numba
    cython

with e.g. `pip install ...` or `conda install`. Then, you can install pamtra2 
with ::

    python setup.py install [flag]

the flag --user will install the package under user's home and it is useful if you do not have root privileges. If you are actively developing, :: 

    python setup.py develop [flag]

is recommended. In case you want to run the tests, you also need :: 

    pytest

For building the documentation, please install :: 

    sphinx
    nbsphinx
    ipykernel
    pandoc
    nbconvert

