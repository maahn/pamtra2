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
    return np.sqrt(eps(temperatures,frequencies,densities,model_mix=model_mix,model_ice=model_ice))

def eps(temperatures,frequencies,densities,model_mix='bruggeman',model_ice='Matzler_2006'):
    fraction = densities/ice_density
    eps_ice = ice.eps(temperatures,frequencies)
    eps_air = complex(1.0,0.0)+0.0*eps_ice
    return mixing.eps([eps_ice,eps_air],[fraction,1.0-fraction],model=model_mix)
