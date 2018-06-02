# -*- coding: utf-8 -*-

import warnings
import collections

"""
Central location for all units
"""


units = collections.defaultdict(
    lambda: 'n/a',
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
