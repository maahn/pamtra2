""" scattering.Rayleigh.py
Rayleigh spherical scatterer object and member functions

"""
import numpy as np
import sys

from .scatterer import scatterer

try:
    from refractive import utilities as ref_utils
except:
    sys.path.append('../../refractive/')
    from refractive import utilities as ref_utils

class Rayleigh(scatterer):
    """
    This is Rayleigh for spheres
    """
    def __init__(self,
                 diameter = 1.0e-3,
                 frequency = None,
                 wavelength = None,
                 refractive_index=None,
                 dielectric_permittivity=None,
                 theta_inc = 0.0,
                 phi_inc = 0.0,
                 theta_sca = 0.0,
                 phi_sca = 0.0):
        scatterer.__init__(self,
                           diameter = diameter,
                           frequency = frequency,
                           refractive_index=refractive_index,
                           dielectric_permittivity=dielectric_permittivity,
                           theta_inc = theta_inc,
                           phi_inc = phi_inc,
                           theta_sca = theta_sca,
                           phi_sca = phi_sca)

        self.geometric_cross_section = np.pi*self.diameter*self.diameter*0.25
        self.K = ref_utils.K(self.dielectric_permittivity)
        
        self.Cabs = 4.0*self.size_parameter*self.K.imag*self.geometric_cross_section
        self.Csca = 8.0*self.size_parameter**4.0*self.K2*self.geometric_cross_section
        self.Cext = self.Cabs + self.Csca
        self.Cbck = 4.0*self.size_parameter**4*self.K2*self.geometric_cross_section
        
        self.S1 = -1.5j*self.size_parameter**3.0
        self.S2 = self.S1*np.cos(self.scatt_angle)
        
        
        #prefactor = np.pi**5 * self.K2 / self.wavelength**4
        #self.back_spec = prefactor*self.diameter**6
