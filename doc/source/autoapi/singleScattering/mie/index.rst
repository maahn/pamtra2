:mod:`singleScattering.mie`
===========================

.. py:module:: singleScattering.mie

.. autoapi-nested-parse::

   singleScattering.mie.py

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



Module Contents
---------------


.. py:class:: MieScatt(diameter=0.001, frequency=None, wavelength=None, refractive_index=None, dielectric_permittivity=None, theta_inc=0.0, phi_inc=0.0, theta_sca=0.0, phi_sca=0.0)

   Bases: :class:`singleScattering.scatterer.Scatterer`

   This is class implement the Mie model of scattering for a sphere


