# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr

import singleScattering

# required because apply_ufunc is picky about args and kwargs...
def _scatteringWrapper(diameter,
                       wavelength,
                       refractiveIndex,
                       model='Rayleigh',
                       ):

    if ((model == 'Rayleigh') or (model == 'Ray')):
        scatt = singleScattering.Rayleigh.RayleighScatt(
                  diameter,
                  wavelength=wavelength,
                  refractive_index=refractiveIndex,
                  )
    elif (model == 'Mie'):
        scatt = singleScattering.Mie.MieScatt(
                  diameter,
                  wavelength=wavelength,
                  refractive_index=refractiveIndex,
                  )
    return np.stack([scatt.Cext, scatt.Csca, scatt.Cabs, scatt.Cbck], axis=-1)


def _dimensionToVariables(darray, dimension, variables):
    profile = {}
    for ii, var in enumerate(variables):
        profile[var] = darray.isel(**{dimension: ii})
    profile = xr.Dataset(profile)
    return profile


def Mie(
  sizeCenter,
  wavelength,
  refractiveIndex,
  ):

    kwargs = dict(model='Mie')
    scatteringProperty = xr.apply_ufunc(
        _scatteringWrapper,
        sizeCenter,
        wavelength,
        refractiveIndex,
        kwargs=kwargs,
        output_core_dims=[['scatteringProperty']],
        output_dtypes=[sizeCenter.dtype],
        output_sizes={'scatteringProperty': 4},
        dask='parallelized',
        vectorize=True,  # mie code needs to accept vectors!
    )

    return scatteringProperty


def Rayleigh(
  sizeCenter,
  wavelength,
  refractiveIndex,
  ):

    kwargs = dict(model='Rayleigh')
    scatteringProperty = xr.apply_ufunc(
        _scatteringWrapper,
        sizeCenter,
        wavelength,
        refractiveIndex,
        kwargs=kwargs,
        output_core_dims=[['scatteringProperty']],
        output_dtypes=[sizeCenter.dtype],
        output_sizes={'scatteringProperty': 4},
        dask='parallelized',
    )

    # variables = ['extinctionCrossSection', 'scatterCrossSection',
    #              'absorptionCrossSection', 'backscatterCrossSection']
    # dimension = 'scatteringProperty'

    return scatteringProperty

