# -*- coding: utf-8 -*-
""" scattering.Mie.py

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

Mie spherical scatterer object and member functions

"""
import sys

from .scatterer import Scatterer
from . import cMie

import numpy as np

try:
    import pamtra2.libs.refractiveIndex.utilities as ref_utils
except:
    sys.path.append('../../refractiveIndex/')
    from refractiveIndex import utilities as ref_utils

class MieScatt(Scatterer):
    """
    This is class implement the Mie model of scattering for a sphere
    """
    def __init__(self,
                 diameter = 1.0e-3,
                 frequency = None,
                 wavelength = None,
                 refractive_index = None,
                 dielectric_permittivity = None,
                 theta_inc = 0.0,
                 phi_inc = 0.0,
                 theta_sca = 0.0,
                 phi_sca = 0.0):
        
        Scatterer.__init__(self,
                           diameter = diameter,
                           frequency = frequency,
                           wavelength = wavelength,
                           refractive_index = refractive_index,
                           dielectric_permittivity = dielectric_permittivity,
                           theta_inc = theta_inc,
                           phi_inc = phi_inc,
                           theta_sca = theta_sca,
                           phi_sca = phi_sca)

        #print('I am a Mie instance')
        self.geometric_cross_section = np.pi*self.diameter*self.diameter*0.25
        self.K = ref_utils.K(self.dielectric_permittivity)
        
        Q = cMie.mie(self.wavelength,self.diameter,self.refractive_index)
        
        self.Cext = Q[0]*self.geometric_cross_section
        self.Csca = Q[1]*self.geometric_cross_section
        self.Cabs = Q[2]*self.geometric_cross_section
        self.Cbck = Q[3]*self.geometric_cross_section
        
        #self.S1 = -1.5j*self.size_parameter**3.0
        #self.S2 = self.S1*np.cos(self.scatt_angle)
        #self.S3 = 0.0
        #self.S4 = 0.0

