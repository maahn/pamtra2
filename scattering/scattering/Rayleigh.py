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
        self.geometric_cross_section = np.pi*self.diameter*self.diameter*0.25
        
        prefactor = np.pi**5 * self.K2 / self.wavelength**4
        self.back_spec = prefactor*self.diameter**6
