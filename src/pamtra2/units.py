# -*- coding: utf-8 -*-

import collections
import warnings

from . import constants

"""
Central location for all units
"""

units = collections.defaultdict(
    lambda: 'n/a',
    {
        'aspectRatio': '-',
        'crossSectionArea': 'm^2',
        'density': 'kg/m^3',
        'frequency': 'Hz',
        'height': 'm',
        'horizontalWind': 'm/s',
        'mass': 'kg',
        'noiseMean': 'dB',
        'sizeDistribution': 'm^(-4)',
        'pathIntegratedAttBottomUp': 'dB',
        'pathIntegratedAttTopDown': 'dB',
        'pressure': 'Pa',
        'relativeHumidity': '%',
        'sizeBounds': 'm',
        'sizeBoundsWidth': 'm',
        'sizeCenter': 'm',
        'specificAttenuation': 'dB/m',
        'temperature': 'K',
        'verticalWind': 'm/s',
        'wavelength': 'm',
        'hydrometeorContent': 'kg/m^3',
    }
)


def kelvin2Celsius(kelvin):
    return kelvin - constants.tFreezing


def celsius2Kelvin(celsius):
    return celsius + constants.tFreezing
