# -*- coding: utf-8 -*-
""" singleScattering.rayleigh.py

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

import sys

import numpy as np
from pamtra2.libs.refractiveIndex import utilities as ref_utils

from .scatterer import Scatterer
from .scattering_utilities import transformation_matrices


class RayleighScatt(Scatterer):
    """
    This is class implement the Rayleigh model of scattering for a sphere
    No check is performed to the actual validity of the Rayleigh approximation
    Inherits from Scatterer, no additional argument with respect to the
    Scatterer class
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
        
        self.Cabs = 4.*self.size_parameter*self.K.imag*self.geometric_cross_section
        self.Csca = 8.*self.size_parameter**4*self.K2*self.geometric_cross_section/3.0
        self.Cext = self.Cabs + self.Csca
        self.Cbck = 4.*self.size_parameter**4*self.K2*self.geometric_cross_section
        
        S1 = self.wavenumber**2*self.K*(self.diameter*0.5)**3
        S2 = S1*np.cos(self.scatt_angle)
        S34 = 0.0 + 0.0j
        Ra, Rb = transformation_matrices(self.rot_alpha, self.rot_beta, self.phi_inc, self.phi_sca)
        self.S = Rb@np.array([[S2, S34], [S34, S1]])@Ra.T
