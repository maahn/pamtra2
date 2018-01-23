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
                 diameter = 1.,
                 frequency = None,
                 wavelength = None,
                 refractive_index=None,
                 dielectric_permittivity=None,):
        scatterer.__init__(self,
                           diameter = diameter,
                           frequency = frequency,
                           refractive_index=refractive_index,
                           dielectric_permittivity=dielectric_permittivity)
        print('I am a Rayleigh instance')
        self.geometric_cross_section = np.pi*self.diameter*self.diameter*0.25
        self.K = ref_utils.K(self.dielectric_permittivity)
        
        alpha = 4.*np.pi*(self.diameter*0.5)**3.0*self.K
        print(self.K,alpha, self.K2)
        self.S1 = -0.25j*alpha*self.wavenumber**3.0/np.pi
        self.theta = np.pi # TODO make it not fixed to backscattering
        self.S2 = self.S1*np.cos(self.theta)
        
        self.Cabs = self.wavenumber*alpha.imag
        self.Csca = self.wavenumber**4.0*(alpha*alpha.conj()).real/(6.0*np.pi)
        self.Cext = self.Csca + self.Cabs
        self.Cbck = 4.0*self.size_parameter**4*self.K2*self.geometric_cross_section
        
        prefactor = np.pi**5 * self.K2 / self.wavelength**4
        self.back_spec = prefactor*self.diameter**6
