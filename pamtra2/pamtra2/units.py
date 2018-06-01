# -*- coding: utf-8 -*-

import warnings

"""
Central location for all units
"""


class unitDict(dict):
    '''
    dict class which returns a default value instead of the item if the
    key does not exist
    '''

    def __init__(self, *args, **kwargs):
        self.default = 'n/a'
        return super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            warnings.warn("no units found for %s." % key)
            return self.default


units = unitDict(
    {
        'aspectRatio': '-',
        'crossSectionArea': 'm^2',
        'density': 'kg/m^3',
        'height': 'm',
        'horizontalWind': 'm/s',
        'mass': 'kg',
        'particleSizeDistribution': 'm^(-4)',
        'pressure': 'Pa',
        'relativeHumidity': '%',
        'sizeCenter': 'm',
        'temperature': 'K',
        'verticalWind': 'm/s',
    }
)
