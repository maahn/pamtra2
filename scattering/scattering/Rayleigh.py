""" scattering.Rayleigh.py
Rayleigh spherical scatterer object and member functions

"""
import numpy as np

from .scatterer import scatterer

class Rayleigh(scatterer):
    """
    This is Rayleigh for spheres
    """
    def __init__(self):
        scatterer.__init__(self)
        print('I am a Rayleigh instance')
        prefactor = np.pi**5 * self.K2 / self.wavelength**4
        self.back_spec = prefactor*self.diameter**6

#  K2 = np.asarray(K2)
#  diameter = np.asarray(diameter)
#back_spec =  prefactor[:,np.newaxis] * diameter**6
