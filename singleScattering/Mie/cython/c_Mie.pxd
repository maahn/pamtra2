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

# This could be part of the  cMie.pyx file, but it could be useful to have a specific namespace (import c_Mie)
# Here we make a wrapper to the functions we want to call from cython

cdef extern from "../src/cMie.h":
    int Mie(double x, double complex m, int nt, double theta[], double complex S1[], double complex S2[]);
    