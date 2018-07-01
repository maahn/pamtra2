
Installation
============

Note that Python 3.6 is required for pamtra2. A few libraries are required. For Ubuntu Trusty (16.04), do ::

    sudo apt-get install gcc gfortran liblapack-dev libfftw3-dev 

For mac OS, make sure the developer tools are installed ::

    xcode-select --install

and use a package manager like https://brew.sh/ to install ::

    brew install lapack fftw gcc

or conda ::

   conda install fftw lapack gcc

In order to install the python dependencies, do ::

    pip install numpy scipy xarray json dask numba cython netcdf4

or ::

    conda install numpy scipy xarray json dask numba cython netcdf4

if you are a conda user. Note that for pip, you likely have to use either `sudo`
or the `--user` flag. 

Then, you can install pamtra2 
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

