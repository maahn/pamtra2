# -*- coding: utf-8 -*-

import warnings
import collections

from . import constants

"""
Central location for all units
"""

units = collections.defaultdict(
    lambda: 'n/a',
    {
        'aspectRatio': '-',
        'attenuation': 'dB/m',
        'crossSectionArea': 'm^2',
        'density': 'kg/m^3',
        'frequency': 'Hz',
        'height': 'm',
        'horizontalWind': 'm/s',
        'mass': 'kg',
        'noiseMean': 'dB',
        'particleSizeDistribution': 'm^(-4)',
        'pathIntegratedAttBottomUp': 'dB',
        'pathIntegratedAttTopDown': 'dB',
        'pressure': 'Pa',
        'relativeHumidity': '%',
        'sizeBounds': 'm',
        'sizeBoundsWidth': 'm',
        'sizeCenter': 'm',
        'temperature': 'K',
        'verticalWind': 'm/s',
        'wavelenght': 'm',

    }
)


def kelvin2Celsius(kelvin):
    return kelvin - constants.tFreezing
