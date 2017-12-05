""" scattering.scatterer.py

This module implements the scatterer class and its scattering model specific subclasses.

"""

import numpy as np
import sys
from . import scattering_utilities

try:
    from refractive import utilities as ref_utils
except:
    sys.path.append('../../refractive/')
    from refractive import utilities as ref_utils
    
light_speed = 299792458.

class scatterer(object):
    """ Parent scatterer class from which every scattering method inherits
    
    """
    def __init__(self,
                 diameter = 1.,
                 frequency = 1.,
                 refractive_index=None,
                 dielectric_permittivity=None,
                 ):
        self.diameter = diameter
        self.frequency = frequency
        self.wavelength = light_speed/frequency
        self.wavenumber = 2.0*np.pi/self.wavelength 
        self.size_parameter = scattering_utilities.size_parameter(0.5*self.diameter,self.wavelength)
        
        self.set_dielectric_properties(refractive_index,dielectric_permittivity)
        
        
        
        
    def set_dielectric_properties(self,refractive_index,dielectric_permittivity):
        """
        Convenient setter of the dielectric properties of the scatterer instance
        """
        if (refractive_index is None):
            if (dielectric_permittivity is None):
                self.refractive_index = np.nan
                self.dielectric_permittivity = np.nan
                self.K2 = np.nan
            else:
                self.dielectric_permittivity = dielectric_permittivity
                self.refractive_index = ref_utils.eps2n(self.dielectric_permittivity)
                self.K2 = ref_utils.K2(self.dielectric_permittivity)
        elif (dielectric_permittivity is None):
            self.refractive_index = refractive_index
            self.dielectric_permittivity = ref_utils.n2eps(self.refractive_index)
            self.K2 = ref_utils.K2(self.dielectric_permittivity)
        else:
            raise AttributeError('Both dielectric permittivity and refractive index have been defined')

class Mie(scatterer):
    def __init__(self):
        scatterer.__init__(self)
        print('I am a Mie instance')
        raise NotImplementedError('Mie is not implemented yet')

class T_Matrix(scatterer):
    def __init__(self):
        scatterer.__init__(self)
        print('I am a T_Matrix instance')
        raise NotImplementedError(scattering_model_name + ' is not implemented yet')

class Liu_DB(scatterer):
    def __init__(self):
        scatterer.__init__(self)
        print('I am a Liu_DB instance')
        raise NotImplementedError(scattering_model_name + ' is not implemented yet')

class Hong_DB(scatterer):
    def __init__(self):
        scatterer.__init__(self)
        print('I am a Hong_DB instance')
        raise NotImplementedError(scattering_model_name + ' is not implemented yet')

class Leinonen_DB(scatterer):
    def __init__(self):
        scatterer.__init__(self)
        print('I am a Leinonen_DB instance')
        raise NotImplementedError(scattering_model_name + ' is not implemented yet')
