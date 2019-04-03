# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr

from ..libs import singleScattering
from .. import constants
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


def _SSRGWrapper(diameter,
                 ssrg_volume,
                 aspect_ratio,
                 wavelength,
                 relativePermittivity,
                 ssrg_parameters='HW14'
                 ):

    scatt = singleScattering.ssrg.SsrgScatt(
        diameter,
        wavelength=wavelength,
        dielectric_permittivity=relativePermittivity,
        volume=ssrg_volume,
        aspect_ratio=aspect_ratio,
        ssrg_parameters=ssrg_parameters
    )

    return np.stack([scatt.Cext, scatt.Csca, scatt.Cabs, scatt.Cbck], axis=-1)


def _TMatrixWrapper(diameter,
                    aspect_ratio,
                    wavelength,
                    relativePermittivity,
                    ):

    scatt = singleScattering.tmatrix.TmatrixScatt(
        diameter,
        wavelength=wavelength,
        dielectric_permittivity=relativePermittivity,
        aspect_ratio=aspect_ratio,
        alpha=0.0,  # needs to be exposed to Pamtra2!
        beta=0.0,  # needs to be exposed to Pamtra2!

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


def SSRG(
    sizeCenter,
    mass,
    aspectRatio,
    wavelength,
    relativePermittivity,
    ssrgParameters='HW14',
):
    """Simple Wrapper for singleScattering.SSRG.SSRGScatt to
    make sure it works with xr.DataArrays.
    """

    kwargs = dict(ssrg_parameters=ssrgParameters)

    volume_ssrg = mass/constants.rhoIce

    scatteringProperty = xr.apply_ufunc(
        _SSRGWrapper,
        sizeCenter,
        volume_ssrg,
        aspectRatio,
        wavelength,
        relativePermittivity,
        kwargs=kwargs,
        output_core_dims=[['scatteringProperty']],
        output_dtypes=[sizeCenter.dtype],
        output_sizes={'scatteringProperty': 4},
        dask='parallelized',
        vectorize=True,
    )

    return scatteringProperty


def TMatrix(
    sizeCenter,
    aspectRatio,
    wavelength,
    relativePermittivity,
):
    """Simple Wrapper for singleScattering.SSRG.SSRGScatt to
    make sure it works with xr.DataArrays.
    """

    kwargs = dict()

    scatteringProperty = xr.apply_ufunc(
        _TMatrixWrapper,
        sizeCenter,
        aspectRatio,
        wavelength,
        relativePermittivity,
        kwargs=kwargs,
        output_core_dims=[['scatteringProperty']],
        output_dtypes=[sizeCenter.dtype],
        output_sizes={'scatteringProperty': 4},
        dask='parallelized',
        vectorize=True,
    )

    return scatteringProperty
