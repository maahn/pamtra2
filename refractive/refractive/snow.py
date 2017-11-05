"""
This module computes snow dielectric properties as a homogeneous mixture of ice and air
or maybe even other stuff ...


"""

import numpy as np
from . import ice 
from . import mixing
# from . import water

#Ice density in [mg/mm^3] [g/cm^3] [kg/dm^3]
ice_density = 0.9167


def n(temperatures,frequencies,densities,model_mix='bruggeman',model_ice='Matzler_2006'):
    """ Effective refractive index of snow according to the specified models for ice
        dielectric properties, effective medium approximation function and effective
        density of the snowflake

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]
    densities: float
        nd array of effective densities [g/cm**3]
    model_mix : string
        Effective Medium Approximation model name default to bruggeman
    model_ice : string
        dielectric model name default to Matzler (2006)
        
    Returns
    -------
    nd - complex
        Refractive index of snow at the requested frequencies and temperatures

    """
    return np.sqrt(eps(temperatures,frequencies,densities,model_mix=model_mix,model_ice=model_ice))

def eps(temperatures,frequencies,densities,model_mix='bruggeman',model_ice='Matzler_2006'):
    """ Effective complex relative dielectric constant of snow according to the specified
        models for ice dielectric properties, effective medium approximation function and
        effective density of the snowflake

    Parameters
    ----------
    temperatures : float
        nd array of temperatures [kelvin]
    frequencies : float
        nd array of frequencies [GHz]
    densities: float
        nd array of effective densities [g/cm**3]
    model_mix : string
        Effective Medium Approximation model name default to bruggeman
    model_ice : string
        dielectric model name default to Matzler (2006)
        
    Returns
    -------
    nd - complex
        Relative dielectric constant of snow at the requested frequencies and temperatures

    """

    fraction = densities/ice_density
    eps_ice = ice.eps(temperatures,frequencies)
    eps_air = complex(1.0,0.0)+0.0*eps_ice
    return mixing.eps([eps_ice,eps_air],[fraction,1.0-fraction],model=model_mix)
