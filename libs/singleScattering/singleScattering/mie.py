# -*- coding: utf-8 -*-
""" singleScattering.mie.py

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

import numpy as np

from pamtra2.libs.refractiveIndex import utilities as ref_utils

from . import cMie
from .scatterer import Scatterer


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
        
        Q, theta, vecS1, vecS2 = cMie.mie(self.wavelength,
                                          self.diameter,
                                          self.refractive_index)
        # Here I apply the dimension and convention conversion factor (-j/k)
        # in order to compare to what Mishenko T-Matrix is giving
        # TODO It might be beneficial if I document the convention somewhere
        # I REALLY DO NOT KNOW WHY: Here I have to switch S1 and S2 and take
        # out the (-) sign. In another program I found (-j/k)... still
        # missing something or messing everything
        self.S1 = 1.j*np.interp(self.scatt_angle, theta, vecS2)/self.wavenumber
        self.S2 = 1.j*np.interp(self.scatt_angle, theta, vecS1)/self.wavenumber
        self.S3 = 0.0
        self.S4 = 0.0

        self.Cext = Q[0]*self.geometric_cross_section
        self.Csca = Q[1]*self.geometric_cross_section
        self.Cabs = Q[2]*self.geometric_cross_section
        self.Cbck = Q[3]*self.geometric_cross_section

