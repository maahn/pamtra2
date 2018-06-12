# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr

from .. import units
from .. import helpers


class instrument(object):
    def __init__(
        self,
        parent,
        frequencies='all',
        **settings,
    ):
        if frequencies == 'all':
            frequencies = parent.profile.frequency
        elif not hasattr(frequencies, '__iter__'):
            frequencies = [frequencies]
        self.frequencies = frequencies
        self.settings = settings
        self.parent = parent

        self.profile = parent.profile.sel(frequency=frequencies)
        self.hydrometeorProfiles = helpers.AttrDict()
        for hydro in parent.hydrometeors.keys():
            self.hydrometeorProfiles[hydro] = parent.hydrometeors[
                hydro].profile.sel(frequency=frequencies)

        self.results = xr.Dataset()

        return self.results
