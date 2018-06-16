
================================================
Welcome to pyPamtraRadarSimulator's documentation!
================================================

:Release: |version|
:Date: |today|

pyPamtraRadarSimulator is the Pamtra Radar Simulator. It simulates a Doppler spectrum based on the backscattering of arbritrary hydrometeors..

Note that pyPamtraRadarSimulator only supports Python version >=3.5.

Installation
============

Dependencies on Ubuntu
**********************

On Ubuntu, install the following packages in advance in order to compile the Fortran Pamtra version::

    sudo apt-get install gcc gfortran liblapack-dev git gitk

The model is tested with gcc version 7.1.0.

To get the Python version, the following packages are additionally required::

    sudo apt-get install python-pandas python-numpy python-matplotlib python-scipy python-netcdf python-pip

The xarray and numba python packages (for exapmples only) might not be in the repositories, install with::

    sudo pip install numba xarray


Dependencies on Mac OS X
************************

On Mac OS X, it is recommended to use brew (http://brew.sh) to install gfortran (via gcc) and netcdf ::

    brew install gcc git
    brew install netcdf --enable-fortran

For the Python version, it is recommended not to use OS X's default python version,
but to install an independent one, e.g. with brew or conda
(https://www.continuum.io/downloads). 
In addition, the following packages are required::

    pip install numpy scipy matplotlib xarray numba


Get model from git repository
*****************************
The version control system git (http://git-scm.com/) is used to keep track of the code. Get a copy of the model with::

    git clone https://github.com/maahn/pyPamtraRadarSimulator

The very basics of git can be found here https://try.github.io/levels/1/challenges/1


Build pyPamtraRadarSimulator
**************
Simply type ::

  python setup.py build

Then, the python routines can be installed with ::

  python setup.py install

if you do not have root permissions you can also use instead of the last line::

    python setup.py install --user


You can start using pyPamtraRadarSimulator in python with ::

  import pyPamtraRadarSimulator

See examples for more information.



Provides Routines
*****************
.. automodule:: pyPamtraRadarSimulator.core
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
