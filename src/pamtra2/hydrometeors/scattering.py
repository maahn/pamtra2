# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr

from ..libs import singleScattering

# required because apply_ufunc is picky about args and kwargs...


def _MieRayleighWrapper(diameter,
                        wavelength,
                        relativePermittivity,
                        model='Rayleigh',
                        ):

    if ((model == 'Rayleigh') or (model == 'Ray')):
        scatt = singleScattering.rayleigh.RayleighScatt(
            diameter,
            wavelength=wavelength,
            dielectric_permittivity=relativePermittivity,
        )
    elif (model == 'Mie'):
        scatt = singleScattering.mie.MieScatt(
            diameter,
            wavelength=wavelength,
            dielectric_permittivity=relativePermittivity,
        )
    return np.stack([scatt.Cext, scatt.Csca, scatt.Cabs, scatt.Cbck], axis=-1)


def Mie(
    sizeCenter,
    wavelength,
    relativePermittivity,
):
    """Simple Wrapper for singleScattering.Mie.MieScatt to
    make sure it works with xr.DataArrays.
    """

    kwargs = dict(model='Mie')
    scatteringProperty = xr.apply_ufunc(
        _MieRayleighWrapper,
        sizeCenter,
        wavelength,
        relativePermittivity,
        kwargs=kwargs,
        output_core_dims=[['scatteringProperty']],
        output_dtypes=[sizeCenter.dtype],
        output_sizes={'scatteringProperty': 4},
        dask='parallelized',
    )

    return scatteringProperty


def Rayleigh(
    sizeCenter,
    wavelength,
    relativePermittivity,
):
    """Simple Wrapper for singleScattering.Rayleigh.RayleighScatt to
    make sure it works with xr.DataArrays.
    """

    kwargs = dict(model='Rayleigh')
    scatteringProperty = xr.apply_ufunc(
        _MieRayleighWrapper,
        sizeCenter,
        wavelength,
        relativePermittivity,
        kwargs=kwargs,
        output_core_dims=[['scatteringProperty']],
        output_dtypes=[sizeCenter.dtype],
        output_sizes={'scatteringProperty': 4},
        dask='parallelized',
    )

    return scatteringProperty
