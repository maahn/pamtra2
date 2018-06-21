# -*- coding: utf-8 -*-
""" scattering.Rayleigh.py

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

Rayleigh spherical scatterer object and member functions

"""

import numpy as np
import sys

from .scatterer import Scatterer

from pamtra2.libs.refractiveIndex import utilities as ref_utils

class RayleighScatt(Scatterer):
    """
    This is class implement the Rayleigh model of scattering for a sphere
    No check is performed to the actual validity of the Rayleigh approximation
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
        
        Scatterer.__init__(self,
                           diameter = diameter,
                           frequency = frequency,
                           wavelength = wavelength,
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
        self.S3 = 0.0
        self.S4 = 0.0
        
        
