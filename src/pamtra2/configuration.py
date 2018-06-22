# -*- coding: utf-8 -*-
# from __future__ import division, absolute_import, print_function
# from builtins import super

import json
from copy import deepcopy

import numpy as np

from . import constants


class SettingsSection(dict):

    """
    Like a dict, but with .freeze() method and .__frozen attribute. If frozen,
    no more keys can be added.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__frozen = False

    def freeze(self):
        self.__frozen = True

    def unfreeze(self):
        self.__frozen = False

    def __setitem__(self, key, value):
        if self.__frozen and key not in self:
            raise ValueError('Frozen dictionary does not contain key %s' % key)
        super().__setitem__(key, value)


DEFAULT_RADAR_SIMULATOR = SettingsSection({
    'randomSeed': 0,
    "aliasingNyquistInterv": 1,
})

DEFAULT_RADAR_PROPERTIES = SettingsSection({
    'integrationTime':  1.4,  # MMCR Barrow during ISDAC
    'fwhrBeamwidthDeg': 0.31,  # full width half radiation beam width MMCR
    # Barrow
    # dielectric constant |K|² (always for liquid water by convention) for
    # the radar equation
    'k2': 0.92,
    # : number of average spectra for noise variance reduction, typical range
    # [1 150]
    "nAve": 150,
    "maxV": 7.885,  # : MinimumNyquistVelocity in m/sec
    "minV": -7.885,  # : MaximumNyquistVelocity in m/sec
    "pNoise1000": -32.23,  # dBz in 1 km height, mean value for BArrow MMCR
    # during ISDAC
    "smoothSpectrum": True,  # Smooth spectrum before estimating moments
})


DEFAULT_HYDROMETEOR_PROPERTIES_BY_TYPE = {
    'liquid': SettingsSection({
        'refractiveIndexModel': "Ellison",
        # 'refractiveIndexMix_snow':None,
        # 'massSizeA_snow' : None,
        # 'massSizeB_snow' : None,
        # 'areaSizeA_snow' : None,
        # 'areaSizeB_snow' : None,
        # 'maxDensity_snow'  :  None,
    }),
    'ice': SettingsSection({
        'refractiveIndexModel': "Matzler_2006",
        # 'refractiveIndexMix_snow':None,
        # 'massSizeA_snow' : None,
        # 'massSizeB_snow' : None,
        # 'areaSizeA_snow' : None,
        # 'areaSizeB_snow' : None,
        # 'maxDensity_snow'  :  None,
    }),
    'snow': SettingsSection({
        'refractiveIndexModel': "Matzler_2006",
        'refractiveIndexMix_snow': "Bruggeman",
        'massSizeA_snow': 0.0121,
        'massSizeB_snow': 1.9,
        'areaSizeA_snow': np.pi,
        'areaSizeB_snow': 2.,
        'maxDensity_snow':  constants.rhoIce,
    }),
}


DEFAULT_GENERAL = SettingsSection({
    'verbosity': 0,
})

DEFAULT_SETTINGS = SettingsSection({
    'general': DEFAULT_GENERAL,
    # 'heymsfield10_particles', #move
    'fallVelocityRelation': 'khvorostyanov01_spheres',
    'radarProperties': SettingsSection({}),
    'hydrometeorProperties': SettingsSection({}),
    'radarSimulator': DEFAULT_RADAR_SIMULATOR,
})


class Settings(SettingsSection):
    def __init__(self, *args):
        """
        Create pamtra2 Settings

        if type(args[0]) is str: load settings from JSON file
        if type(args[0]) is dict: get settings from arguments
        if list of frequencies,list of hydrometeors: populate
        settings with default values including radar defualt
        for each frequency and hydrometeor default for each hydrometeor.
        """
        if isinstance(args[0], dict):
            super().__init__(args[0])
        elif isinstance(args[0], str):
            with open(args[0], 'r') as jsonfile:
                settings = json.load(jsonfile, object_hook=SettingsSection)
            super().__init__(settings)
        else:
            super().__init__(deepcopy(DEFAULT_SETTINGS))
            freqs, hydrometeors = args
            for freq in freqs:
                self['radarProperties'][repr(freq)] = deepcopy(
                    DEFAULT_RADAR_PROPERTIES)
                self['radarProperties'][repr(freq)].freeze()
            for hydrometeor in hydrometeors:
                self['hydrometeorProperties'][hydrometeor] = {}
                # self['hydrometeorProperties'][hydrometeor].freeze()

        for k in self.keys():
            try:
                self[k].freeze()
            except AttributeError:
                print(k)
                pass

    def __setitem__(self, key, value):
        raise ValueError('Only children of this dictionary can be modified')

    def to_json(self, fname):
        with open(fname, 'w') as jsonfile:
            json.dump(self, jsonfile)
