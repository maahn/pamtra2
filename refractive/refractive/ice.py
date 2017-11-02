# -*- coding: utf-8 -*-
""" refractive.ice module.

This module provides a list of ice refractive index models to compute the
dielectric properties of ice according to the requested frequencies and
temeperatures.
The module can be also used as a standalone python script.

Example
-------
The python script is callable as

    $ python ice.py Temperature Frequency

and returns the complex refractive index of ice at the requested
Temperature [Kelvin] and Frequency [GHz]

Notes
-----
    It is possible to call the functions implemented in this module using
    nd-arrays. The function arguments must either have exactly the same
    shape allowing element-wise application of the functions or one of
    the two must be a scalar which will be spread across the nd computations

Temperatures should be provided in Kelvin and frequencies in GHz
A very basic non-negative value check is performed on the data

"""

import numpy as np

def warren_brandt_2008(temperatures,frequencies):
    print("Not implemented yet falling back to Matzler 2006")
    return matzler_2006(temperatures,frequencies)

def matzler_2006(temperatures,frequencies):
    """Ice complex relative dielectric constant according to Matzler (2006)
    "Microwave remote sensing" TODO: correct book title

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]

    Returns
    -------
    nd - complex
        Relative dielectric constant of ice at the requested frequencies and temperatures

    Raises
    ------
    ValueError
        If a negative frequency or temperature is passed as an argument

    """

    if (frequencies < 0).any():
        raise ValueError('A negative frequency value has been passed')
    if (temperatures < 0).any():
        raise ValueError('A negative temperature value has been passed')

    B1 = 0.0207
    b = 335.
    B2 = 1.16e-11
    c = 299792458.

    eps1  = 3.1884+(temperatures-273)*9.1e-4
    theta = 300./temperatures-1.
    alpha =(0.00504+0.0062*theta)*np.exp(-22.1*theta)
    deltabeta=np.exp(-9.963+0.0372*(temperatures-273.16))
    betaM = B1*np.exp(b/temperatures)/(temperatures*((np.exp(b/temperatures)-1)*(np.exp(b/temperatures)-1)))+B2*frequencies*frequencies
    beta  = betaM+deltabeta
    eps2  = alpha/frequencies + beta*frequencies
    return eps1 + 1j*eps2

#######################################################################################################

def eps(temperatures,frequencies,model="Matzler_2006"):
    """Ice complex relative dielectric constant according to the requested model

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]
    model : string
        dielectric model name default to Matzler (2006)

    Returns
    -------
    nd - complex
        Relative dielectric constant of ice at the requested frequencies and temperatures

    Raises
    ------
    ValueError
        If a negative frequency or temperature is passed as an argument

    """
    if (model == "Matzler_2006"):
        return matzler_2006(temperatures,frequencies)
    else:
        print("I do not recognize the ice refractive index specification, falling back to Matzler 2006")
        return matzler_2006(temperatures,frequencies)

def n(temperatures,frequencies,model="Matzler_2006"):
    """Ice complex refractive index according to the requested model

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]
    model : string
        dielectric model name default to Matzler (2006)

    Returns
    -------
    nd - complex
        Refractive index of ice at the requested frequencies and temperatures

    Raises
    ------
    ValueError
        If a negative frequency or temperature is passed as an argument

    """
    return np.sqrt(eps(temperatures,frequencies,model))

#######################################################################################################

if __name__ == "__main__":
    import sys
    n(float(sys.argv[1]),float(sys.argv[2]))
