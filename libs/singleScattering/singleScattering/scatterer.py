""" scattering.scatterer.py

    Copyright (C) 2017 - 2018 Davide Ori dori@uni-koeln.de
    Institute for Geophysics and Meteorology - University of Cologne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


This module implements the scatterer class which is the parent class for
all of the specific scattering models.
The parent scatterer class is intended to provide an abstraction layer or a
common interface to the various scattering models.

"""

import numpy as np
import sys
from . import scattering_utilities as scatt_utils

try:
    import pamtra2.libs.refractiveIndex.utilities as ref_utils
except:
    sys.path.append('../../refractiveIndex/')
    from refractiveIndex import utilities as ref_utils
    
light_speed = 299792458.

class Scatterer(object):
    """ Parent Scatterer class from which every scattering method inherits
    
    Attributes:
        diameter: Equivalent size of the target in meters. The user should know
            how to properly define it.
        frequency, wavelength: The frequency [Hz] or the wavelength [meters] of
            incident light, see set_electromagnetic_wave function for details.
        refractive_index, dielectric_permittivity: The refractive index or the
            relative dielectric permittivity of the scattering target, see
            set_dielectric_properties for further details
        alpha, beta, gamma: Euler angles of the particle orientation... 
                            ... still do not know if I need all of them
        theta_inc, theta_sca:
        phi_inc, phi_sca:
    
    """
    def __init__(self,
                 diameter = 1.,
                 frequency = None,
                 wavelength = None,
                 refractive_index = None,
                 dielectric_permittivity = None,
                 alpha = 0.0,
                 beta = 0.0,
                 theta_inc = 0.0,
                 theta_sca = 0.0,
                 phi_inc = 0.0,
                 phi_sca = 0.0
                 ):
        self.diameter = diameter
        
        self.set_electromagnetic_wave(wavelength, frequency)
        self.wavenumber = 2.0*np.pi/self.wavelength
        self.size_parameter = scatt_utils.size_parameter(0.5*diameter,
                                                         self.wavelength)
        
        self.set_dielectric_properties(refractive_index,
                                       dielectric_permittivity)
        
        self.set_scattering_geometry([theta_inc, theta_sca, phi_inc, phi_sca,
                                      alpha, beta])
       
        
    def set_electromagnetic_wave(self, wavelength, frequency):
        """ Convenience setter of the properties of the incoming electromagnetic
        wave

        This setter function resolve the ambiguity of specifing either 
        wavelength or the frequency which should not be set indipendently in 
        orderto avoid internal inconsistencies
        
        Parameters
        ----------
        wavelength : scalar real wavelength [meters] of the incoming 
            electromagnetic wave

        frequency : scalar real frequency [Hz] of the incoming electromagnetic
            wave
        
        """
        if (wavelength is None):
            if (frequency is None):
                raise AttributeError('Either frequency or wavelength' +
                    'must be set')
            else:
                self.frequency = frequency
                self.wavelength = light_speed/frequency
        elif (frequency is None):
            self.wavelength = wavelength
            self.frequency = light_speed/wavelength
        else:
            raise AttributeError('Both frequency and wavelength have been'
                + 'defined')
                
        
    def set_dielectric_properties(self, refractive_index, 
                                  dielectric_permittivity):
        """ Convenience setter of the dielectric properties of the scatterer
        instance
        
        This setter resolve the ambiguity of specifing either dielectric 
        permittivity or the refractive index of the medium which should not be
        set independently in order to avoid internal inconsistencies
        
        Parameters
        ----------
        refractive_index : scalar complex refractive index [dimensionless] of 
            the scattering target

        dielectric_permittivity : scalar complex relative dielectric 
            permittivity [dimensionless] of the scattering target
        
        """
        
        if (refractive_index is None):
            if (dielectric_permittivity is None):
                raise AttributeError('Dielectric permittivity or refractive' +
                    ' index should be defined')
            else:
                self.dielectric_permittivity = np.array(dielectric_permittivity)
                self.refractive_index = ref_utils.eps2n(self.dielectric_permittivity)
                self.K2 = ref_utils.K2(self.dielectric_permittivity)
        elif (dielectric_permittivity is None):
            self.refractive_index = np.array(refractive_index)
            self.dielectric_permittivity = ref_utils.n2eps(self.refractive_index)
            self.K2 = ref_utils.K2(self.dielectric_permittivity)
        else:
            raise AttributeError('Both dielectric permittivity and refractive'
                +' index have been defined')
            
    def set_scattering_geometry(self,geometry):
        """ Convenience setter of the scattering geometry that takes as input
        a 4-element array containing all 4 incident and scattering angles
        
        Any update of the geometry should call the scattering_angle function in
        order to correctly update this value

        Parameters
        ----------
        geometry : A tuple of 6 elements containing (theta_inc, theta_sca, 
            phi_inc, phi_sca, alpha, beta) [rad]
            
        """
        
        (self.theta_inc, self.theta_sca, self.phi_inc, self.phi_sca,
         self.alpha, self.beta) = geometry
        self.scatt_angle = scatt_utils.scattering_angle(self.theta_inc,
                                                        self.theta_sca,
                                                        self.phi_inc,
                                                        self.phi_sca)



class T_Matrix(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a T_Matrix instance')
        raise NotImplementedError('T_Matrix is not implemented yet')

class Liu_DB(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a Liu_DB instance')
        raise NotImplementedError('Liu_DB is not implemented yet')

class Hong_DB(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a Hong_DB instance')
        raise NotImplementedError('Hong_DB is not implemented yet')

class Leinonen_DB(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a Leinonen_DB instance')
        raise NotImplementedError('Leinonen_DB is not implemented yet')

class Aydin_DB(Scatterer):
    def __init__(self):
        Scatterer.__init(self)
        print('I am a Aydin_DB instance')
        raise NotImplementedError('Aydin_DB is not implemented yet')
