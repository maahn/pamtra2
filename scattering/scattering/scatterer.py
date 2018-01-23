""" scattering.scatterer.py

This module implements the scatterer class and its scattering model specific subclasses.

"""

import numpy as np
import sys
from . import scattering_utilities as scatt_utils

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
                 frequency = None,
                 wavelength = None,
                 refractive_index=None,
                 dielectric_permittivity=None,
                 theta_inc = 0.0,
                 phi_inc = 0.0,
                 theta_sca = 0.0,
                 phi_sca = 0.0
                 ):
        self.diameter = diameter
        
        self.set_electromagnetic_wave(wavelength, frequency)
        self.wavenumber = 2.0*np.pi/self.wavelength
        self.size_parameter = scatt_utils.size_parameter(0.5*self.diameter,self.wavelength)
        
        self.set_dielectric_properties(refractive_index,dielectric_permittivity)
        
        self.set_scattering_geometry([theta_inc,phi_inc,theta_sca,phi_sca])
       
        
    def set_electromagnetic_wave(self, wavelength,frequency):
        """ Convenient setter of the properties of the incoming electromagnetic wave
        
        """
        if (wavelength is None):
            if (frequency is None):
                raise AttributeError('Either frequency or wavelength must be set')
            else:
                self.frequency = frequency
                self.wavelength = light_speed/frequency
        elif (frequency is None):
            self.wavelength = wavelength
            self.frequency = light_speed/wavelength
        else:
            raise AttributeError('Both frequency and wavelength have been defined')
                
        
    def set_dielectric_properties(self, refractive_index, dielectric_permittivity):
        """ Convenient setter of the dielectric properties of the scatterer instance
        """
        if (refractive_index is None):
            if (dielectric_permittivity is None):
                raise AttributeError('Dielectric permittivity or refractive index should be defined')
            else:
                self.dielectric_permittivity = np.array(dielectric_permittivity)
                self.refractive_index = ref_utils.eps2n(self.dielectric_permittivity)
                self.K2 = ref_utils.K2(self.dielectric_permittivity)
        elif (dielectric_permittivity is None):
            self.refractive_index = np.array(refractive_index)
            self.dielectric_permittivity = ref_utils.n2eps(self.refractive_index)
            self.K2 = ref_utils.K2(self.dielectric_permittivity)
        else:
            raise AttributeError('Both dielectric permittivity and refractive index have been defined')
            
    def set_scattering_geometry(self,geometry):
        """ Convenient setter of the scattering geometry that takes as input
        a 4-element array containing all 4 incident and scattering angles
            
        """
        self.theta_inc = geometry[0]
        self.phi_inc = geometry[1]
        self.theta_sca = geometry[2]
        self.phi_sca = geometry[3]
        self.scatt_angle = scatt_utils.scattering_angle(self.theta_inc,
                                                        self.phi_inc,
                                                        self.theta_sca,
                                                        self.phi_sca)



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
