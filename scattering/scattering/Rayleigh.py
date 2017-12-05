""" scattering.Rayleigh.py
Rayleigh spherical scatterer object and member functions

"""
import numpy as np

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
    def __init__(self):
        scatterer.__init__(self)
        print('I am a Rayleigh instance')
        self.geometric_cross_section = np.pi*self.diameter*self.diameter*0.25
        self.K = ref_utils.K(self.eps)
        
        alpha = 4.*np.pi*(self.diameter*0.5)**3.0*self.K
        self.S1 = -0.25j*alpha*self.k**3.0/np.pi
        self.theta = np.pi # TODO make it not fixed to backscattering
        self.S2 = S1*np.cos(theta)
        
        self.Cabs = self.k*alpha.imag
        self.Csca = self.k**4.0*(alpha*alpha.conj).real/(6.0*np.pi)
        self.Cext = self.Csca + self.Cabs
        self.Cbck = 4.0*self.x**4*self.K2*self.geometric_cross_section
        
        prefactor = np.pi**5 * self.K2 / self.wavelength**4
        self.back_spec = prefactor*self.diameter**6
