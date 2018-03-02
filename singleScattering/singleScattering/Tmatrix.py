# -*- coding: utf-8 -*-
""" scattering.Tmatrix.py

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

Tmatrix spheroid scatterer object and member functions

"""

import numpy as np
import sys

from .scatterer import Scatterer

try:
    from refractiveIndex import utilities as ref_utils
except:
    sys.path.append('../../refractiveIndex/')
    from refractiveIndex import utilities as ref_utils

class T_Matrix(Scatterer):
    """
    This is class implement the Tmatrix model of scattering for a spheroid
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
                 phi_sca = 0.0
                 aspect_ratio=0.0):
        
        Scatterer.__init__(self,
                           diameter = diameter,
                           frequency = frequency,
                           refractive_index=refractive_index,
                           dielectric_permittivity=dielectric_permittivity,
                           theta_inc = theta_inc,
                           phi_inc = phi_inc,
                           theta_sca = theta_sca,
                           phi_sca = phi_sca)
        print('I am a T_Matrix instance')
        raise NotImplementedError('T_Matrix is not implemented yet')
