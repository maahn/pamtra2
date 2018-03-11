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

#cython: boundscheck=False
#cython: wraparound=False
# Comments above are special. Please do not remove.
cimport numpy as np  # needed for function arguments
import numpy as np # needed for np.empty_like

# Import the namespace c_Mie described in c_Mie.pxd
cimport c_Mie

ctypedef np.float32_t float_t
ctypedef np.float64_t double_t
ctypedef np.int32_t int_t
ctypedef np.complex128_t complex_t

#def mie(double_t x, complex_t m):#double_t wl,double_t size,complex_t m, th=None):
def mie(double_t wavelength, double_t size, complex_t m):
    print('I am python3 c cmie and refractive index is ',m,' and wl is ',wavelength,' and size is ', size)
    cdef double_t x
    x=np.pi*size/wavelength
    print('I am python3 c cmie and refractive index is ',m,' and x is ',x)
    cdef int_t nt
    nt = 180
    cdef np.ndarray[dtype=double_t, ndim=1, mode="c"] theta
    theta = np.linspace(0.0, np.pi, nt, dtype='d')
    cdef np.ndarray[dtype=complex_t, ndim=1, mode="c"] S1
    S1 = np.zeros(nt, dtype=np.complex128)
    cdef np.ndarray[dtype=complex_t, ndim=1, mode="c"] S2
    S2 = np.zeros(nt, dtype=np.complex128)
    c_Mie.Mie(x, m, nt, <double*>theta.data, <double complex*> S1.data, <double complex*> S2.data)
    print(S1[0],S2[0])
    print(S1[-1],S2[-1])

def mie_coated(np.ndarray[dtype= double_t,ndim=1,mode='c'] x,
               np.ndarray[dtype=complex_t,ndim=1,mode='c'] m):
    raise NotImplementedError('mie_coated does nothing at the moment except receiving arguments')
    print('I am python3 c cmie and refractive index is ',m)
    c_Mie.Mie(m)

def mie_Nlayers(np.ndarray[dtype= double_t,ndim=1,mode='c'] x,
                np.ndarray[dtype=complex_t,ndim=1,mode='c'] m):
    raise NotImplementedError('mie_Nlayers does nothing at the moment except receiving arguments')
    print('I am python3 c cmie and refractive index is ',m)
    c_Mie.Mie(m)
