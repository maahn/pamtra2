""" refractive.utilities module
This module provides a short list of utilities and complementary functions
for the refractive index module.

Basic conversion from refractive index to dielectric permittivity
(and viceversa) is implemented.
The module also provides a conversion function from dielectric permittivity to
radar dielectric factor K2 which is of great importance in radar applications

"""

import numpy as np

eps2n = lambda eps: np.sqrt(eps)

n2eps = lambda n: n*n

def K(eps):
    """ Rayleigh complex dielectric factor
    This is basically the K complex factor that defines the Radar dielectric
    factor |K|**2. It is useful in Rayleigh theory to define absorption cross
    section from its imaginary part
    
    Parameters
    ----------
    eps : complex
        nd array of complex relative dielectric constants

    Returns
    -------
    nd - float
        Rayleigh complex dielectric factor K
    """
    return (eps-1.0)/(eps+2.0)

def K2(eps):
    """ Radar dielectric factor |K|**2

    Parameters
    ----------
    eps : complex
        nd array of complex relative dielectric constants

    Returns
    -------
    nd - float
        Radar dielectric factor |K|**2 real

    """
    K_complex = (eps-1.0)/(eps+2.0)
    return (K_complex*K_complex.conj()).real
