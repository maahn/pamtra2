""" refractive.ice module.

This module provides a list of water refractive index models to compute the
dielectric properties of water according to the requested frequencies and
temeperatures.
The module can be also used as a standalone python script.

Example
-------
The python script is callable as

    $ python water.py Temperature Frequency

and returns the complex refractive index of water at the requested
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

def ellison(temperatures,frequencies):
    """Water complex relative dielectric constant according to Ellison (2005)
    "..." TODO: put the extensive correct reference here

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

    a0 = 5.7230
    a1 = 2.2379e-2
    a2 = -7.1237e-4
    a3 = 5.0478
    a4 = -7.0315e-2
    a5 = 6.0059e-4
    a6 = 3.6143
    a7 = 2.8841e-2
    a8 = 1.3652e-1
    a9 = 1.4825e-3
    a10 = 2.4166e-4

    T = temperatures-273.15
    es=(37088.6-82.168*T)/(421.854+T)
    einf=a6+a7*T
    e1=a0+T*(a1+T*a2)              #a0+a1*T+a2*T*T
    ni1=(45.0+T)/(a3+T*(a4+T*a5))  #(a3+a4*T+a5*T*T)
    ni2=(45.0+T)/(a8+T*(a9+T*a10)) #(a8+a9*T+a10*T*T)
    A1=frequencies/ni1
    A2=frequencies/ni2
    eps1=(es-e1)/(1+A1*A1)+(e1-einf)/(1+A2*A2)+einf
    eps2=(es*A1-e1*A1)/(1+A1*A1)+(e1*A2-einf*A2)/(1+A2*A2)
    return eps1 + 1j*eps2

def pamtra_water(Temperatures,Frequencies):
    return ellison(np.array(Temperatures),np.array(Frequencies))
# PLACEHOLDER FOR WHAT PAMTRA IS CURRENTLY COMPUTING

#######################################################################################################

def eps(Temperatures,Frequencies,model="ellison"):
    """Water complex relative dielectric constant according to the requested model

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]
    model : string
        dielectric model name default to Ellison (2005)

    Returns
    -------
    nd - complex
        Relative dielectric constant of water for the requested frequencies and temperatures

    Raises
    ------
    ValueError
        If a negative frequency or temperature is passed as an argument

    """
    if (model == "ellison"):
        return ellison(np.array(Temperatures),np.array(Frequencies))
    else:
        print("I do not recognize the ice refractive index specification, falling back to ellison")
        return ellison(np.array(Temperatures),np.array(Frequencies))

def n(Temperatures,Frequencies,model="ellison"):
    """Water complex refractive index according to the requested model

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]
    model : string
        dielectric model name default to Ellison (2005)

    Returns
    -------
    nd - complex
        Refractive index of water for the requested frequencies and temperatures

    Raises
    ------
    ValueError
        If a negative frequency or temperature is passed as an argument

    """
    return np.sqrt(eps(np.array(Temperatures),np.array(Frequencies),model))

#######################################################################################################

if __name__ == "__main__":
    import sys
    n(float(sys.argv[1]),float(sys.argv[2]))
