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

import sys

import numpy as np
from pamtra2.libs.refractiveIndex import utilities as ref_utils

from .scatterer import Scatterer
from pamtra2.libs.singleScattering import fTMat

from scipy.integrate import dblquad

RADIUS_EQUAL_VOLUME = 1.0
RADIUS_EQUAL_AREA = 0.0
RADIUS_MAXIMUM = 2.0

SHAPE_SPHEROID = -1
SHAPE_CYLINDER = -2
SHAPE_CHEBYSHEV = 1

class TmatrixScatt(Scatterer):
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
                 phi_sca = 0.0,
                 alpha = 0.0, # we introduce alpha and beta euler angles
                 beta = 0.0,  # for orientation
                 aspect_ratio=1.0):
        
        Scatterer.__init__(self,
                           diameter = diameter,
                           frequency = frequency,
                           refractive_index=refractive_index,
                           dielectric_permittivity=dielectric_permittivity,
                           theta_inc = theta_inc,
                           phi_inc = phi_inc,
                           theta_sca = theta_sca,
                           phi_sca = phi_sca)
                                
        self.geometric_cross_section = np.pi*self.diameter*self.diameter*0.25
        self.K = ref_utils.K(self.dielectric_permittivity)
        self.radius = 0.5*self.diameter
        self.radius_type = RADIUS_EQUAL_VOLUME
        self.aspect_ratio = aspect_ratio
        self.shape = SHAPE_SPHEROID
        self.ddelt = 1e-3
        self.ndgs = 2
        self.alpha = alpha
        self.beta = beta
        
        self._init_tmatrix()
        self.S, self.Z = self.get_SZ()
        print(self.S,self.Z)
        self.Csca = self.scattering_xsect()
        self.Cext = 2.*self.wavelength*self.S[1,1].imag  # horizontal polarization
        self.Cbck = 2.*np.pi*(self.Z[0,0]-self.Z[0,1]-self.Z[1,0]+self.Z[1,1])  # horizontal polarization
        self.Cabs = self.Cext-self.Csca 

    def _init_tmatrix(self):
        """Initialize the T-matrix.
        """

        if self.radius_type == RADIUS_MAXIMUM:
            # Maximum radius is not directly supported in the original
            # so we convert it to equal volume radius
            radius_type = Scatterer.RADIUS_EQUAL_VOLUME
            radius = self.equal_volume_from_maximum()
        else:
            radius_type = self.radius_type
            radius = self.radius
        
        m = self.refractive_index
        self.nmax = fTMat.calctmat(radius, radius_type, self.wavelength,
                                   self.refractive_index.real,
                                   self.refractive_index.imag,
                                   self.aspect_ratio, self.shape, self.ddelt,
                                   self.ndgs)

    def get_SZ(self, alpha=None, beta=None):
        """Get the S and Z matrices for a single orientation.
        """
        if alpha == None:
            alpha = self.alpha
        if beta == None:
            beta = self.beta
            
        (self._S_single, self._Z_single) = fTMat.calcampl(self.nmax,
                                                          self.wavelength,
                                                          self.theta_inc,
                                                          self.theta_sca,
                                                          self.phi_inc,
                                                          self.phi_sca,
                                                          alpha, beta)
        return (self._S_single, self._Z_single)

    def scattering_xsect(self):
        """Calculates the scattering cross section by integrating over the all
           the whole 4pi solid scattering angle
        """
        old_theta = self.theta_sca
        old_phi = self.phi_sca
        
        def diff_xsect(theta, phi):
            """Differential scattering cross section multiplied by sin(theta) to be
               integrated over 4pi in order to get the scattering cross section
            """
            self.theta_sca = theta
            self.phi_sca = phi        
            S, Z = self.get_SZ()        
            I = Z[0, 0] - Z[0, 1] # horizontal polarization
            return I * np.sin(theta)

        xsect = dblquad(diff_xsect, 0.0, 2*np.pi,
                        lambda x: 0.0, lambda x: np.pi)[0]
        self.theta_sca = old_theta
        self.phi_sca = old_phi
        return xsect

