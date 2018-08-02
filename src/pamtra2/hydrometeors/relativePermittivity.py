# -*- coding: utf-8 -*-
import warnings

import numpy as np
import xarray as xr

from ..libs import refractiveIndex
from .. import constants

'''Wrapper module to handle the refractiveIndex library from within Pamtra2.
'''

# Some work out of the box with Pamtra2:
water_turner_kneifel_cadeddu = refractiveIndex.water.turner_kneifel_cadeddu
water_ellison = refractiveIndex.water.ellison
ice_warren_brandt_2008 = refractiveIndex.ice.warren_brandt_2008


# These routines need to be wrapped to handle xarray objects.
def ice_matzler_2006(
        temperature,
        frequency,
        checkTemperatureForRelativePermittivity=False
):

    if (temperature < 240).any():
        warnings.warn('ice_matzler_2006 defined only above > 240K')

    return refractiveIndex.ice.matzler_2006(
        temperature,
        frequency,
        checkTemperature=checkTemperatureForRelativePermittivity
    )


def ice_iwabuchi_yang_2011(temperature, frequency):
    args = xr.broadcast(temperature, frequency)
    relativePermittivity = xr.apply_ufunc(
        refractiveIndex.ice.iwabuchi_yang_2011,
        *args,
        kwargs={},
        output_dtypes=[np.complex],
        dask='parallelized',
        vectorize=True,
    )
    return relativePermittivity


# Copy doc strings
ice_matzler_2006.__doc__ = refractiveIndex.ice.matzler_2006.__doc__
ice_iwabuchi_yang_2011.__doc__ = refractiveIndex.ice.iwabuchi_yang_2011.__doc__


# Even though the mixing formulas work out of the box, we modify them to
# account for density as input
def mixing_maxwell_garnett(relativePermittivityIce, density):
    mix = density/constants.rhoIce
    assert (np.asarray(mix) <= 1).all()
    assert (np.asarray(mix) > 0).all()
    return refractiveIndex.mixing.maxwell_garnett(relativePermittivityIce, mix)


def mixing_bruggeman(relativePermittivityIce, density):
    mix = density/constants.rhoIce
    assert (np.asarray(mix) <= 1).all()
    assert (np.asarray(mix) > 0).all()
    return refractiveIndex.mixing.bruggeman(relativePermittivityIce, mix)


def mixing_sihvola(relativePermittivityIce, density):
    mix = density/constants.rhoIce
    assert (np.asarray(mix) <= 1).all()
    assert (np.asarray(mix) > 0).all()
    return refractiveIndex.mixing.sihvola(relativePermittivityIce, mix)


# Copy doc strings
mixing_maxwell_garnett.__doc__ = refractiveIndex.mixing.maxwell_garnett.\
    __doc__.replace('eps', 'relativePermittivityIce').replace('mix', 'density')
mixing_bruggeman.__doc__ = refractiveIndex.mixing.bruggeman.__doc__.replace(
    'eps', 'relativePermittivityIce').replace('mix', 'density')
mixing_sihvola.__doc__ = refractiveIndex.mixing.sihvola.__doc__.replace(
    'eps', 'relativePermittivityIce').replace('mix', 'density')
