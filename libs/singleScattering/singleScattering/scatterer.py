""" scattering.scatterer.py

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


This module implements the scatterer class which is the parent class for
all of the specific scattering models.
The parent scatterer class is intended to provide an abstraction layer or a
common interface to the various scattering models.

"""

import sys

import numpy as np
from pamtra2.libs.refractiveIndex import utilities as ref_utils

from . import scattering_utilities as scatt_utils

light_speed = 299792458.


class Scatterer(object):
    """ Parent Scatterer class from which every scattering method inherits

    Attributes:
        diameter: Equivalent size of the target in meters. The user should know
            how to properly define it.
        frequency, wavelength: The frequency [Hz] or the wavelength [meters] of
            incident light, see set_electromagnetic_wave function for details.
        refractive_index, dielectric_permittivity: The refractive index or the
            relative dielectric permittivity of the scattering target, see
            set_dielectric_properties for further details
        alpha, beta, gamma: Euler angles of the particle orientation... 
                            ... still do not know if I need all of them
        theta_inc, theta_sca:
        phi_inc, phi_sca:

    """

    def __init__(self,
                 diameter=1.,
                 frequency=None,
                 wavelength=None,
                 refractive_index=None,
                 dielectric_permittivity=None,
                 alpha=0.0,
                 beta=0.0,
                 theta_inc=0.0,
                 theta_sca=0.0,
                 phi_inc=0.0,
                 phi_sca=0.0,
                 aspect_ratio=1,
                 volume=-99,
                 ):

        # first, convert inputs to arrays
        diameter, diameter_scalar = self.makeArray(diameter)
        refractive_index, refractive_index_scalar = self.makeArray(
            refractive_index)
        dielectric_permittivity, dielectric_permittivity_scalar = self.makeArray(
            dielectric_permittivity)
        frequency, frequency_scalar = self.makeArray(frequency)
        wavelength, wavelength_scalar = self.makeArray(wavelength)
        aspect_ratio, aspect_ratio_scalar = self.makeArray(aspect_ratio)
        volume, volume_scalar = self.makeArray(volume)

        if (
            diameter_scalar and
            refractive_index_scalar and
            dielectric_permittivity_scalar and
            frequency_scalar and
            wavelength_scalar and
            aspect_ratio_scalar and
            volume_scalar
        ):
            self.scalar_input = True
        else:
            self.scalar_input = False

        self.diameter = diameter
        self.aspect_ratio = aspect_ratio
        self.volume = volume
        
        self.set_electromagnetic_wave(wavelength, frequency)

        self.set_dielectric_properties(refractive_index,
                                       dielectric_permittivity)

        self.ravel_input()

        self.wavenumber = 2.0*np.pi/self.wavelength
        self.size_parameter = scatt_utils.size_parameter(0.5*self.diameter,
                                                         self.wavelength)
        self.set_scattering_geometry([theta_inc, theta_sca, phi_inc, phi_sca,
                                      alpha, beta])

    def makeArray(self, var):
        # to make sure it is an array
        if var is None:
            scalar_input = True
        else:
            var = np.asarray(var)
            scalar_input = False
            if var.ndim == 0:
                var = var[np.newaxis]  # Makes x 1D
                scalar_input = True
        return var, scalar_input

    def ravel_input(self):

        # make sure input is flatt

        out = np.broadcast_arrays(
            self.diameter,
            self.refractive_index,
            self.K2,
            self.dielectric_permittivity,
            self.frequency,
            self.wavelength,
            self.aspect_ratio,
            self.volume,
        )
        (
            self.diameter,
            self.refractive_index,
            self.K2,
            self.dielectric_permittivity,
            self.frequency,
            self.wavelength,
            self.aspect_ratio,
            self.volume,
        ) = out

        self.shapeIn = self.diameter.shape

        self.diameter = self.diameter.ravel()
        self.refractive_index = self.refractive_index.ravel()
        self.K2 = self.K2.ravel()
        self.dielectric_permittivity = self.dielectric_permittivity.ravel()
        self.frequency = self.frequency.ravel()
        self.wavelength = self.wavelength.ravel()
        self.aspect_ratio = self.aspect_ratio.ravel()
        self.volume = self.volume.ravel()



    def unravel_output(self):

        # black to original shape!

        if self.scalar_input:
            # if scalars were initialy provided, make arrays scalar again

            self.S = np.squeeze(self.S)
            self.Cabs = np.squeeze(self.Cabs)
            self.Csca = np.squeeze(self.Csca)
            self.Cext = np.squeeze(self.Cext)
            self.Cbck = np.squeeze(self.Cbck)
        else:
            self.S = self.S.reshape(self.shapeIn+(2, 2,))
            self.Cabs = self.Cabs.reshape(self.shapeIn)
            self.Csca = self.Csca.reshape(self.shapeIn)
            self.Cext = self.Cext.reshape(self.shapeIn)
            self.Cbck = self.Cbck.reshape(self.shapeIn)

    def set_electromagnetic_wave(self, wavelength, frequency):
        """ Convenience setter of the properties of the incoming electromagnetic
        wave

        This setter function resolve the ambiguity of specifing either 
        wavelength or the frequency which should not be set indipendently in 
        orderto avoid internal inconsistencies

        Parameters
        ----------
        wavelength : scalar real wavelength [meters] of the incoming 
            electromagnetic wave

        frequency : scalar real frequency [Hz] of the incoming electromagnetic
            wave

        """
        if (wavelength is None):
            if (frequency is None):
                raise AttributeError('Either frequency or wavelength' +
                                     'must be set')
            else:
                self.frequency = frequency
                self.wavelength = light_speed/frequency
        elif (frequency is None):
            self.wavelength = wavelength
            self.frequency = light_speed/wavelength
        else:
            raise AttributeError('Both frequency and wavelength have been'
                                 + 'defined')

    def set_dielectric_properties(self, refractive_index,
                                  dielectric_permittivity):
        """ Convenience setter of the dielectric properties of the scatterer
        instance

        This setter resolve the ambiguity of specifing either dielectric 
        permittivity or the refractive index of the medium which should not be
        set independently in order to avoid internal inconsistencies

        Parameters
        ----------
        refractive_index : scalar complex refractive index [dimensionless] of 
            the scattering target

        dielectric_permittivity : scalar complex relative dielectric 
            permittivity [dimensionless] of the scattering target

        """

        if (refractive_index is None):
            if (dielectric_permittivity is None):
               # raise AttributeError('Dielectric permittivity or refractive' +
               #    ' index should be defined')
               # need to allow this for DB based scattering ???
                self.refractive_index = None
                self.dielectric_permittivity = None
            else:
                self.dielectric_permittivity = np.array(
                    dielectric_permittivity)
                self.refractive_index = ref_utils.eps2n(
                    self.dielectric_permittivity)
                self.K2 = ref_utils.K2(self.dielectric_permittivity)
        elif (dielectric_permittivity is None):
            self.refractive_index = np.array(refractive_index)
            self.dielectric_permittivity = ref_utils.n2eps(
                self.refractive_index)
            self.K2 = ref_utils.K2(self.dielectric_permittivity)
        else:
            raise AttributeError('Both dielectric permittivity and refractive'
                                 + ' index have been defined')

    def set_scattering_geometry(self, geometry):
        """ Convenience setter of the scattering geometry that takes as input
        a 4-element array containing all 4 incident and scattering angles

        Any update of the geometry should call the scattering_angle function in
        order to correctly update this value

        Parameters
        ----------
        geometry : A tuple of 6 elements containing (theta_inc, theta_sca, 
            phi_inc, phi_sca, alpha, beta) [rad]

        """

        (self.theta_inc, self.theta_sca, self.phi_inc, self.phi_sca,
         self.alpha, self.beta) = geometry
        angles = scatt_utils.scattering_angle(self.theta_inc, self.theta_sca,
                                              self.phi_inc, self.phi_sca)
        self.scatt_angle, self.rot_alpha, self.rot_beta = angles

    def estimate_amplitude_matrix(self, S1, S2, S34, Ra, Rb):

        S1234 = np.zeros(S1.shape + (2, 2,)) + 0.0j
        S1234[..., 0, 0] = S2
        S1234[..., 0, 1] = S34
        S1234[..., 1, 0] = S34
        S1234[..., 1, 1] = S1

        # @ operator applies only to two last dimesnions.
        self.S = Rb@S1234@Ra.T


class Liu_DB(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a Liu_DB instance')
        raise NotImplementedError('Liu_DB is not implemented yet')


class Hong_DB(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a Hong_DB instance')
        raise NotImplementedError('Hong_DB is not implemented yet')


class Leinonen_DB(Scatterer):
    def __init__(self):
        Scatterer.__init__(self)
        print('I am a Leinonen_DB instance')
        raise NotImplementedError('Leinonen_DB is not implemented yet')


class Aydin_DB(Scatterer):
    def __init__(self):
        Scatterer.__init(self)
        print('I am a Aydin_DB instance')
        raise NotImplementedError('Aydin_DB is not implemented yet')
