
.. module:: pygasabs


====================================
Welcome to pygasabs's documentation!
====================================

:Release: |version|
:Date: |today|

pygasabs is the Pamtra's module to estimate geaseous absorption in the microwave range. 


Installation
============

Dependencies
************

Python 3 and Numpy are required


Get model from git repository
*****************************
The version control system git (http://git-scm.com/) is used to keep track of the code. Get a copy of the model with::

    git clone https://github.com/maahn/pamtra2/pygasabs

The very basics of git can be found here https://try.github.io/levels/1/challenges/1


Build pygasabs
**************
Simply type ::

  python setup.py build

Then, the python routines can be installed with ::

  python setup.py install

if you do not have root permissions you can also use instead of the last line::

    python setup.py install --user


You can start using pamRaMo in python with ::

  import pygasabs

See examples for more information.


Provided Routines
*****************

.. autosummary::
   :toctree: generated/

.. currentmodule:: pygasabs.core

.. autofunction:: calculate_gas_absorption

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

