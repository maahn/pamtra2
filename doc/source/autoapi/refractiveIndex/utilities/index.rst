:mod:`refractiveIndex.utilities`
================================

.. py:module:: refractiveIndex.utilities

.. autoapi-nested-parse::

   refractive.utilities module

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

   This module provides a short list of utilities and complementary functions
   for the refractive index module.

   Basic conversion from refractive index to dielectric permittivity
   (and viceversa) is implemented.
   The module also provides a conversion function from dielectric permittivity to
   radar dielectric factor K2 which is of great importance in radar applications



Module Contents
---------------


.. data:: speed_of_light
   :annotation: = 299792458.0

   

.. function:: eps2n(eps)


.. function:: n2eps(n)


.. function:: wavenumber(frequency=None, wavelength=None)


.. function:: K(eps)

   Rayleigh complex dielectric factor
   This is basically the K complex factor that defines the Radar dielectric
   factor |K|**2. It is useful in Rayleigh theory to define absorption cross
   section from its imaginary part

   :param eps: nd array of complex relative dielectric constants
   :type eps: complex

   :returns: Rayleigh complex dielectric factor K
   :rtype: nd - float


.. function:: K2(eps)

   Radar dielectric factor |K|**2

   :param eps: nd array of complex relative dielectric constants
   :type eps: complex

   :returns: Radar dielectric factor |K|**2 real
   :rtype: nd - float


