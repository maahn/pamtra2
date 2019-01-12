# -*- mode: python -*-
# -*- coding: utf-8 -*-
"""
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
"""

# cython: boundscheck=False
# cython: wraparound=False
# Comments above are special. Please do not remove.
cimport numpy as np  # needed for function arguments
import numpy as np  # needed for np.empty_like

# Import the namespace c_Mie described in c_Mie.pxd
cimport c_Mie

ctypedef np.float32_t float_t
ctypedef np.float64_t double_t
ctypedef np.int32_t int_t
ctypedef np.complex128_t complex_t


def mie(wavelength, size, m,  nangles=180):


    nIterations = len(wavelength)

    Q = np.zeros((nIterations, 4,), dtype=np.float64)
    S1 = np.zeros((nIterations, nangles,), dtype=np.complex128)
    S2 = np.zeros((nIterations, nangles,), dtype=np.complex128)

    for ii in range(nIterations):
        Q[ii], theta, S1[ii], S2[ii] = mie_one(
            wavelength[ii],  size[ii],  m[ii], nangles=nangles
            )

    return Q, theta, S1, S2



def mie_one(double_t wavelength, double_t size, complex_t m, nangles=180):
    """ This is a python high level interface to the C version Mie included in
    c_Mie external module.

        Parameters
        ----------
        wavelength : scalar-double
            The wavelength of the incoming electromagnetic radiation.
            Can be in any measuring unit but we suggest SI [meters]

        size : scalar-double (TODO become array)
            Diameter of the scattering sphere. For consistency it must have the
            same measuring unit of the wavelength.

        m : scalar-complex (TODO become array)
            The complex refractive index of the scattering sphere

        nangles : scalar-integer
            Number of angles to partition the 0-pi range for the calculation of
            the elements of the amplitude matrix. By default it is set to 180,
            so S is computed every 1 deg, but can be increased for accuracy in
            postprocessing interpolation

        Returns
        -------
        Q : array(4)-double (to become nd-double)
            Array of 4 efficiencies Qext, Qsca, Qabs, Qbck. To be multiplied by the
            geometric corss-section to get Cext, Csca, Cabs and Cbck

        theta : array-double [rad]
            Array of angles in radians at which the element of the complex amplitude
            matrix S1 and S2 are computed
        
        S1 : array(nangles)-complex (TODO become nd-complex)
            S1 elements of the amplitude matrix (S22 element in Bohren and Huffman)
            one for each scattering angle

        S2 : array(nangles)-complex (TODO become nd-complex)
            S2 elements of the amplitude matrix (S11 element in Bohren and Huffman)
            one for each scattering angle

    """
    cdef double_t x
    x = np.pi * size / wavelength

    cdef int_t nt
    nt = nangles
    cdef np.ndarray[dtype = double_t, ndim = 1, mode = "c"] theta
    theta = np.linspace(0.0, np.pi, nt, dtype='d')

    cdef np.ndarray[dtype = complex_t, ndim = 1, mode = "c"] S1
    S1 = np.zeros(nt, dtype=np.complex128)
    cdef np.ndarray[dtype = complex_t, ndim = 1, mode = "c"] S2
    S2 = np.zeros(nt, dtype=np.complex128)

    cdef np.ndarray[dtype = double_t, ndim = 1, mode = "c"] Q
    Q = np.zeros(4, dtype='d')

    c_Mie.Mie(x, m, nt, < double * > theta.data,
              < double complex * > S1.data, < double complex * > S2.data,
              < double * > Q.data)

    return Q, theta, S1, S2


def mie_coated(np.ndarray[dtype=double_t, ndim=1, mode='c'] x,
               np.ndarray[dtype=complex_t, ndim=1, mode='c'] m):
    raise NotImplementedError(
        'mie_coated does nothing at the moment except receiving arguments')
    print('I am python3 c cmie and refractive index is ', m)
    c_Mie.Mie(m)


def mie_Nlayers(np.ndarray[dtype=double_t, ndim=1, mode='c'] x,
                np.ndarray[dtype=complex_t, ndim=1, mode='c'] m):
    raise NotImplementedError(
        'mie_Nlayers does nothing at the moment except receiving arguments')
    print('I am python3 c cmie and refractive index is ', m)
    c_Mie.Mie(m)
