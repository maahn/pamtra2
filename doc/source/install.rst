
Installation
============

First install all required python packages ::

    numpy
    scipy
    xarray
    json
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

