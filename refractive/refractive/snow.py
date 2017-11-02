"""
This module computes snow dielectric properties as a homogeneous mixture of ice and air
or maybe even other stuff ...


"""

from . import ice mixing #water

#Ice density in g/cm^3
ice_density = 0.9167


def n(density):
    return mixing.n(density)

def eps(density):
    return mixing.eps(density)

