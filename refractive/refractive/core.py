# -*- coding: utf-8 -*-
""" refractive module core submodule

This module provides a list of ice and water refractive index models to compute 
the dielectric properties of ice according to the requested frequencies and
temeperatures. The module is completed with some Effective Medium approximation
functions to compute the refractive index of inhomogeneous mixtures of materials
which are directly used to compute the dielectric properties of snow as a
dilution of ice in air.

This core file loads submodules and provide handy functions to
consistently call ice, water or snow refractive index modules

Example
-------
    $ python
    >>> import refractive
    >>> refractive.n(temperatures, frequencies, **kwargs)

and returns the complex refractive index of ice at the requested
Temperature [Kelvin] and Frequency [Hz]

Notes
-----
    It is possible to call the functions implemented in this module using
    nd-arrays. The function arguments must either have exactly the same
    shape allowing element-wise application of the functions or one of
    the two must be a scalar which will be spread across the nd computations

    Frequencies and Temperatures are always mandatory arguments as name of the
    substance, but specific algorithms requires special additional arguments to
    be passed in order to proceed (for instance snow density must be defined).
    The functions check for argument consistency and raise AttributeError if
    a wrong list of attributes is passed.

All of the argument quatities must be provided in SI units: Temperatures in 
Kelvin, frequencies in Hz, densities in kg/m3. The specific called algorithm 
check for arguments values to be within the limits of validity of the dielectric
model and raises ValueError in case they are not respected

"""

from . import ice
from . import snow
from . import water
from . import utilities

import numpy as np

substances_list = ['ice','water','snow']
#argument_list = ['temperatures','frequencies','densities','model','model_mix',
#                   'model_ice']
#argument_dict = {}
#for attr in argument_list:
#    argument_dict[attr] = None

def n(substance, temperatures, frequencies,**kwargs):# model=None, model_ice=None, model_mix=None, densities=None):
    """Complex index of refraction of the requested substance according to the 
        requested specifications

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [Kelvin]
    frequencies : float
        nd array of frequencies [Hz]
    **kwargs : additional arguments to be passed to the requested model

    Returns
    -------
    nd - complex
        Refractive index of the requested substance using the requested options

    Raises
    ------
    AttributeError
        If an uncorrect list of arguments is passed

    """
    return np.sqrt(eps(substance,temperatures,frequencies,**kwargs))

def eps(substance, temperatures, frequencies,**kwargs):# model=None, model_ice=None, model_mix=None, densities=None):
    """Complex relative dielectric permittivity of the requested substance 
        according to the requested specifications

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [Kelvin]
    frequencies : float
        nd array of frequencies [Hz]
    **kwargs : additional arguments to be passed to the requested model

    Returns
    -------
    nd - complex
        Refractive index of the requested substance using the requested options

    Raises
    ------
    AttributeError
        If an uncorrect list of arguments is passed

    """

#    for attr in argument_list:
#        if attr in kwargs:
#            self.__dict__[attr] = kwargs[attr]

    if (substance == 'ice'):
        if 'model' in kwargs.keys():
            model = kwargs['model']
        return ice.eps(temperatures,frequencies,**kwargs)#model=model)
    elif (substance == 'water'):
        if 'model' in kwargs.keys():
            model = kwargs['model']
        return water.eps(temperatures,frequencies,**kwargs)#model=model)
    elif (substance == 'snow'):
        return snow.eps(temperatures,frequencies,**kwargs)#densities=densities,model_mix=model_mix,model_ice=model_ice)
    else:
        raise AttributeError("I do not recognize the " + str(substance) + "as a " + 
        "valid substance I can only compute dielectric properties of " + str(substances_list))
