# -*- coding: utf-8 -*-
'''Routines to create Pamtra2 objects from radiosondes.

'''
import numpy as np
import xarray as xr

from .. import units
from .. import core


def arm_interpolatedsonde(
    interpolatedsonde,
    maxHeight,
    frequencies,
    hydrometeors=[],
):

    '''Import ARM's interpolatedsonde product

    Parameters
    ----------

    interpolatedsonde : xr.Dataset
        Dataset with radiosonde data. E.g. result of
        xr.open_mfdataset('/path/*nc')
    maxHeight : float
        maximum height to consider in m
    frequencies : list of float, optional
        list of frequencies to consider
    hydrometeors : list of hydrometeor names for preparing the data structure
        default []
    Returns
    -------

    Pamtra2 object :
      Pamtra2 object with radiosonde data.
    '''

    profile = interpolatedsonde.isel(height=np.where(
        interpolatedsonde.height < maxHeight/1000)[0])

    for var in ['temp', 'rh', 'bar_pres', 'wspd']:

        profile[var] = profile[var].where(profile['qc_%s' % var] == 0)

    renameDict = {
        'temp': 'temperature',
        'rh': 'relativeHumidity',
        'bar_pres': 'pressure',
        'wspd': 'horizontalWind',
        'height': 'layer',
    }
    profile = profile.rename(
        renameDict
    )
    profile = profile[list(renameDict.values())]

    profile['temperature'] = units.celsius2Kelvin(profile.temperature)
    profile['pressure'] = profile.pressure * 1000
    profile['layer'] = profile.layer * 1000

    profile['height'] = profile.layer
    profile['layer'].values = np.arange(len(profile.layer))

    profile['frequency'] = xr.DataArray(
        frequencies, dims=['frequency'], coords=[frequencies])

    pam2 = core.pamtra2(
        profile=profile,
        frequencies=frequencies,
        nLayer=len(profile.height),
        additionalDims={'time': profile.time},
        hydrometeors=hydrometeors,
    )

    return pam2
