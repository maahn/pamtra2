# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot)
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)


from __future__ import absolute_import, division, print_function
import numpy as np

# from .due import due, Doi

import meteo_si.constants
import meteo_si.humidity


__all__ = ["kelvin_2_celsius", "celsius_to_kelvin", "T_virt_rh", "T_virt_q"]


'''
Functions to deal with temperature.

'''


def kelvin_2_celsius(T):
    """
    Convert the temperature from Kelvin to Celsius.

    Parameters
    ----------
    T
       Temperature in Kelvin.

    Returns
    -------

    C
        Temperature in Celsius.

    """

    return T + meteo_si.constants.Tnull


def celsius_to_kelvin(C):
    """
    Convert the temperature from Celsius to Kelvin.

    Parameters
    ----------

    C
        Temperature in Celsius.

    Returns
    -------
    T
       Temperature in Kelvin.
    """

    return C - meteo_si.constants.Tnull


def T_virt_rh(T, rh, p):
    '''
    Calculate the virtual temperature from air temperature,
    pressure, and relative humidity.

    Parameters
    ----------
    T
        Temperature in in K
    rh
        relative humidity in Pa/Pa
    p
        pressure in Pa

    Returns
    -------
    T_virt
        Virtual temperature in K
    '''

    with np.errstate(divide='ignore', invalid='ignore'):
        if np.any(rh > 5):
            raise TypeError("rh must not be in %")

    return T_virt_q(T, meteo_si.humidity.rh2q(rh, T, p))


def T_virt_q(T, q):
    '''
    Calculate the virtual temperature from air temperature and specific
    humidity.

    Parameters
    ----------
    T
        Temperature in in K
    q
        specific humidity in kg/kg

    Returns
    -------
    T_virt
        Virtual temperature in K
    '''
    return T + T * (meteo_si.constants.Rvapor/meteo_si.constants.Rair-1) * q
