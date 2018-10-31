# -*- coding: utf-8 -*-
""" single scattering module core submodule

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

This module provides a list of scattering technique suitable (and commonly used)
for the computation of the single scattering properties of hydrometeors at
microwave frequencies; future upgrades might consider inclusion of further
frequency bands.
This core file loads submodules and provide handy functions to
consistently call the various scattering models

Example
-------
    $ python
    >>> import singleScattering
    >>> singleScattering.scattering(frequency, size, refractive_index, **kwargs)
	and returns the most important scattering quantities Cext, Csca, Cabs, Cbk

Notes
-----
    At the moment it is not possible to call the functions implemented in this 
    module using nd-arrays, but it is on due for the next development steps.
    This feature will come of particular handy when integration over PSD will
    come to play

All of the argument quatities must be provided in SI units
"""

# Generic
import numpy as np

# Complementary library
# Implemented with basic functionality
from . import mie, rayleigh, tmatrix, ssrg, scatterer, scattering_utilities

models_list = ['Rayleigh (Ray)', 'Mie', 'Tmatrix (TMM)', 
               'Self-Similar Rayleigh-Gans (SSRG)', 'LiuDB', 'LeinonenDB', 
               'AydinDB', 'HongDB', 'ChalmersDB', 'OpenSSP (KwoDB)']

def scattering(diameters,
               frequencies=None, wavelengths=None,
               refractive_indices=None, dielectric_permittivities=None,
               orientation=None, # placeholder, look at what is implemented in Tmatrix
               model=None,
               **kwargs):
    """Scattering properties according to the passed scatterer parameters

    Parameters
    ----------
    diameters : float
        sizes of the scattering particle [meters]
    frequencies : float
        frequencies of the incoming electromagnetic wave [Hz]
        might be substituted by wavelengths
    wavelengths : float
        wavelengths of the incoming electromagnetic wave [Hz]
        might be substituted by frequencies [meters]
    refractive_indices : complex 
        effective refractive index of the scattering material (with respect to 
        ambient refractive index) can be substituted by dielectric_permittivities
    dielectric_permittivities : complex
        effective relative dielectric permittivity of the scattering material 
        (with respect to ambient dielectric properties) can be substituted by 
        refractive_indices
    orientation : placeholder for orientation of the scatterer
    model : one of the model_list (might be substituted by the shortname in the parenthesis)
    
    **kwargs : additional arguments to be passed to the requested model

    Returns
    -------
    nd - float
        array of scattering properties (Cext, Csca, Cabs, Cbck) where
        Cext : Total extinction cross section [meters**2] in the direction 
               of propagation
        Csca : Total extinction cross section [meters**2] in the direction 
               of propagation
        Cabs : Total extinction cross section [meters**2] in the direction 
               of propagation
        Cbck : Radar backscattering cross section [meters**2]
               => 4*pi dCsca(pi)/dOmega

    Raises
    ------
    AttributeError
        If an uncorrect list of arguments is passed

    """

    if (orientation is not None):
        raise NotImplementedError("At the moment the library only computes non oriented properties")

    if ((model == 'Rayleigh') or (model == 'Ray')):
        scatt = rayleigh.RayleighScatt(diameters, frequencies, wavelengths,
                                       refractive_indices,
                                       dielectric_permittivities, **kwargs)
    elif (model == 'Mie') or (model == 'MIE'):
        scatt = mie.MieScatt(diameters, frequencies, wavelengths,
                             refractive_indices, dielectric_permittivities,
                             **kwargs)
    elif ((model == 'Tmatrix') or (model == 'TMM')):
        scatt = tmatrix.TmatrixScatt(diameters, frequencies, wavelengths,
                                     refractive_indices,
                                     dielectric_permittivities, **kwargs)
    elif ( (model == 'Self-Similar Rayleigh-Gans') or (model=='SSRG') ):
        scatt = ssrg.SsrgScatt(diameters, frequencies, wavelengths,
                               refractive_indices, dielectric_permittivities,
                               **kwargs)
    elif ((model == 'LiuDB') ):
        raise NotImplementedError(model+" capabilities will be soon added")
    elif ((model == 'LeinonenDB') ):
        raise NotImplementedError(model+" capabilities will be soon added")
    elif ((model == 'AydinDB') ):
        raise NotImplementedError(model+" capabilities will be soon added")
    elif ((model == 'HongDB') ):
        raise NotImplementedError(model+" capabilities will be soon added")
    elif ((model == 'ChalmersDB') ):
        raise NotImplementedError(model+" capabilities will be soon added")
    elif ((model == 'OpenSSP') or (model == 'KwoDB')):
        raise NotImplementedError(model+" capabilities will be soon added")
    else:
        raise AttributeError("I do not recognize the %s as a " 
                             "valid substance I can only compute"
                             " dielectric properties of %s" % (
                                model, models_list))
    return scatt.Cext, scatt.Csca, scatt.Cabs, scatt.Cbck, scatt.S
