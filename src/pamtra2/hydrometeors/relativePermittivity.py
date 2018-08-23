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


# mixing wrapper
def _mixing_wrapper(eps1, density1, func=None):
    mix1 = density1/constants.rhoIce
    mix2 = 1.0-mix1
    assert (np.asarray(mix1) <= 1).all()
    assert (np.asarray(mix1) > 0).all()

    eps2 = complex(1.0, 0.0)+(0.0*eps1)

    mix = (mix1, mix2)
    eps = (eps1, eps2)
    return func(eps, mix)


def mixing_sihvola(relativePermittivityIce, density):
    args = xr.broadcast(relativePermittivityIce, density)
    relativePermittivity = xr.apply_ufunc(
        _mixing_wrapper,
        *args,
        kwargs={'func': refractiveIndex.mixing.sihvola},
        output_dtypes=[np.complex],
        dask='parallelized',
        vectorize=True,
    )
    return relativePermittivity


def mixing_bruggeman(relativePermittivityIce, density):
    args = xr.broadcast(relativePermittivityIce, density)
    relativePermittivity = xr.apply_ufunc(
        _mixing_wrapper,
        *args,
        kwargs={'func': refractiveIndex.mixing.bruggeman},
        output_dtypes=[np.complex],
        dask='parallelized',
        vectorize=True,
    )
    return relativePermittivity


def mixing_maxwell_garnett(relativePermittivityIce, density):
    args = xr.broadcast(relativePermittivityIce, density)
    relativePermittivity = xr.apply_ufunc(
        _mixing_wrapper,
        *args,
        kwargs={'func': refractiveIndex.mixing.maxwell_garnett},
        output_dtypes=[np.complex],
        dask='parallelized',
        vectorize=True,
    )
    return relativePermittivity


# Copy doc strings
mixing_maxwell_garnett.__doc__ = refractiveIndex.mixing.maxwell_garnett.\
    __doc__.replace('eps', 'relativePermittivityIce').replace('mix', 'density')
mixing_bruggeman.__doc__ = refractiveIndex.mixing.bruggeman.__doc__.replace(
    'eps', 'relativePermittivityIce').replace('mix', 'density')
mixing_sihvola.__doc__ = refractiveIndex.mixing.sihvola.__doc__.replace(
    'eps', 'relativePermittivityIce').replace('mix', 'density')
