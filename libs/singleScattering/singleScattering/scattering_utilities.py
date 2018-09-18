""" singleScattering.scattering_utilities.py

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

This submodule implements useful functions that helps in the definition of the
scatterer geometry and the conversion among diverse scattering frameworks
"""

import numpy as np


class amplitude_matrix(object):
    """ Complex 2x2 amplitude matrix as defined in Bohren and Huffman
    | S2   S3 |
    |         |
    | S4   S1 |


    """

    def __init__(self, S):
        self.matrix = S + 0.0j  # ensure it is complex valued

    @property
    def S1(self):
        return self.matrix[1, 1]

    @property
    def S2(self):
        return self.matrix[0, 0]

    @property
    def S3(self):
        return self.matrix[0, 1]

    @property
    def S4(self):
        return self.matrix[1, 0]

    def to_mueller(self):
        return amplitude2mueller(self)


class scattering_matrix(object):
    """ Real 4x4 scattering matrix as defined in Bohren and Huffman

    """

    def __init__(self, Z):
        # ensure it is real (TODO maybe it is better to just have a check)
        self.matrix = Z.real


def amplitude2mueller(ampl):
    """ This function implement the conversion between complex 2x2 amplitude matrix
    to the real 4x4 scattering Mueller matrix according to Bohren Huffman pp 65
    actually according to Mishchenko (2000) pp 51-52
    It should not be dependent on the reference fo S, so it is invariant by rotations
    around the propagation vector
    However, better to check, I see a sign problem in Z13 and Z14, also Z23 and Z24
    and probably more ...

    """

    mueller = np.ndarray((4, 4))
    S1_2 = (ampl.S1 * ampl.S1.conjugate()).real
    S2_2 = (ampl.S2 * ampl.S2.conjugate()).real
    S3_2 = (ampl.S3 * ampl.S3.conjugate()).real
    S4_2 = (ampl.S4 * ampl.S4.conjugate()).real
    mueller[0, 0] = 0.5 * (S2_2 + S1_2 + S4_2 + S3_2)
    mueller[0, 1] = 0.5 * (S2_2 - S1_2 + S4_2 - S3_2)
    mueller[0, 2] = -(ampl.S2 * ampl.S3.conjugate() + ampl.S1 * ampl.S4.conjugate()).real
    mueller[0, 3] = -(ampl.S2 * ampl.S3.conjugate() - ampl.S1 * ampl.S4.conjugate()).imag

    mueller[1, 0] = 0.5 * (S2_2 - S1_2 - S4_2 + S3_2)
    mueller[1, 1] = 0.5 * (S2_2 + S1_2 - S4_2 - S3_2)
    mueller[1, 2] = -(ampl.S2 * ampl.S3.conjugate() - ampl.S1 * ampl.S4.conjugate()).real
    mueller[1, 3] = -(ampl.S2 * ampl.S3.conjugate() + ampl.S1 * ampl.S4.conjugate()).imag

    mueller[2, 0] = -(ampl.S2 * ampl.S4.conjugate() + ampl.S1 * ampl.S3.conjugate()).real
    mueller[2, 1] = -(ampl.S2 * ampl.S4.conjugate() - ampl.S1 * ampl.S3.conjugate()).real
    mueller[2, 2] = (ampl.S2 * ampl.S1.conjugate() + ampl.S3 * ampl.S4.conjugate()).real
    mueller[2, 3] = (ampl.S2 * ampl.S1.conjugate() + ampl.S4 * ampl.S3.conjugate()).imag

    mueller[3, 0] = -(ampl.S4 * ampl.S2.conjugate() + ampl.S1 * ampl.S3.conjugate()).imag
    mueller[3, 1] = -(ampl.S4 * ampl.S2.conjugate() - ampl.S1 * ampl.S3.conjugate()).imag
    mueller[3, 2] = (ampl.S1 * ampl.S2.conjugate() - ampl.S3 * ampl.S4.conjugate()).imag
    mueller[3, 3] = (ampl.S1 * ampl.S2.conjugate() - ampl.S3 * ampl.S4.conjugate()).real

    return scattering_matrix(mueller)

size_parameter = lambda radius, wavelength: 2.0 * np.pi * radius / wavelength


class spheroid(object):

    def __init__(self):
        raise NotImplementedError

    def mass(self):
        raise NotImplementedError

    def volume(self):
        raise NotImplementedError

    def density(self):
        raise NotImplementedError

    def effective_volume_diameter(self):
        raise NotImplementedError


def scattering_angle(theta_inc, theta_sca, phi_inc, phi_sca):
    """ Calculates the scattering angle in radians given the full set of four
    angles that defines the scattering geometry from the incident and the
    scattered wave directions

    The reference frame is assumed to be set with the polar angle theta
    measured from the vertical axis z, and the azimuth angle phi measured from
    the axis x???

    It also returns the rotation angle alpha which is the angle between the
    scattering plane and the plane defined by the scattering direction and z.
    This angle is needed to transform S matrix components from the Bohren
    and Huffman convention (parallel and perpendicular) to the Mishchenko
    convention (theta and phi) and viceversa

    """
    sin_inc = np.sin(theta_inc)
    sin_sca = np.sin(theta_sca)
    phidiff = abs(phi_sca - phi_inc)
    if (phidiff > np.pi):
        phidiff = 2.*np.pi-phidiff # we need phidiff [0,pi]
    phidiff = abs(phi_sca - phi_inc)
    acos_th = sin_inc * sin_sca * np.cos(phidiff) \
              + np.cos(theta_inc) * np.cos(theta_sca)
    th = np.arccos(acos_th)
    factor = np.sin(phidiff)/np.sin(th)
    alpha = np.arcsin(sin_sca*factor)
    beta = np.arcsin(sin_inc*factor)
    if np.isnan(alpha):
        alpha = 0.0
    if np.isnan(beta):
        beta = 0.0
    return th, alpha, beta


def rotation_matrix(alpha, beta):
    """ This function returns  two 2x2 rotation matrix specified by the
    arguments alpha and beta, it should comes handy when you have to transform
    amplitude matrices from the Bohren and Huffman (parallel, perpendicular)
    basis to the Mishchenko (theta, phi) convention
    Rb is actually an improper rotation (rotoreflection)

    """
    sina = np.sin(alpha)
    cosa = np.cos(alpha)
    sinb = np.sin(beta)
    cosb = np.cos(beta)

    return np.array([[cosa, -sina], [-sina, -cosa]]), np.array([[-cosb, sinb], [-sinb, -cosb]])