# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

import numpy as np
import xarray as xr
import scipy.special
import numba

# all functions have to accept hydrometeor size as first parameter

from . import mass


def modifiedGamma(sizeCenter, N0, lambd, mu, gamma):
    """
    classical modifed gamma distribution

    Parameters
    ----------
    sizeCenter : array_like
      particle size at center of size bin
    N0 : array_like
      N0 prefactor (default None)
    lambd : float or array_like
      lambda parameter (default array)
    mu : float or array_like
      mu parameter (default array)
    gamma : float or array_like
      gamma parameter (default array)

    Returns
    -------

    N : array
      particle size distribution with shape = N0.shape + sizeCenter.shape
    """

    N = N0 * sizeCenter**mu * np.exp(-lambd * sizeCenter**gamma)

    return N


def gamma(sizeCenter, N0, lambd, mu):
    """
    classical gamma distribution
    """
    gamma = 1.
    return modifiedGamma(
      sizeCenter=sizeCenter, N0=N0, lambd=lambd, mu=mu, gamma=gamma
      )


def exponential(sizeCenter, N0, lambd):
    """
    classical exponential distribution
    """
    mu = 0.
    gamma = 1.
    return modifiedGamma(
      sizeCenter=sizeCenter, N0=N0, lambd=lambd, mu=mu, gamma=gamma
      )


def exponentialField(sizeCenter, temperature, lambd):
    """
    classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T)
    """
    N0 = _exponentialField(temperature)
    return exponential(sizeCenter, N0, lambd)


def exponentialFieldWC(
  sizeCenter, temperature, waterContent, massSizeA, massSizeB
  ):
    """
    classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T). Lambda is
    estimated from the WC which requires a power law mass-size
    relation.
    """

    N0 = _exponentialField(temperature)
    lambd = _exponentialWC2Lambda(N0, waterContent, massSizeA, massSizeB)
    return exponential(sizeCenter, N0, lambd)


def exponentialFieldReff(sizeCenter, temperature, effectiveRadius):
    """
    classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T). Lambda is
    estimated from the effective radius.
    """
    N0 = _exponentialField(temperature)
    lambd = _exponentialReff2Lambda(effectiveRadius)
    return exponential(sizeCenter, N0, lambd)


def exponentialN0WC(
    sizeCenter,
    N0,
    waterContent,
    massSizeA=mass.powerLawLiquidPrefactor,
    massSizeB=mass.powerLawLiquidExponent,
):
    """
    classical exponential distribution constrained with N0 and LWC.
    """

    lambd = _exponentialWC2Lambda(N0, waterContent, massSizeA, massSizeB)
    return exponential(sizeCenter, N0, lambd)


########################################################################
# Helper functions
########################################################################


def _exponentialField(temperature):
    """
    N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T)
    """
    if not np.all(temperature > 0):
        raise ValueError('temperature has to be larger 0K')
    if np.any(np.isnan(temperature)):
        raise ValueError('Found NAN in temperature')

    N0 = 7.628e6 * np.exp(0.107 * (273.15 - temperature))

    if np.any(np.isnan(N0)):
        raise ValueError('Found NAN in N0')

    return N0


def _exponentialWC2Lambda(N0, WC, massSizeA, massSizeB):
    """
    Estimate lambda of exponential distribution from WC and Ntot
    """

    if np.any(np.isnan(WC)):
        raise ValueError('Found NAN in WC')
    if np.any(np.isnan(massSizeA)):
        raise ValueError('Found NAN in massSizeA')
    if np.any(np.isnan(massSizeB)):
        raise ValueError('Found NAN in massSizeB')

    lambd = (massSizeA * N0 * scipy.special.gamma(massSizeB+1.) /
             WC)**(1. / (massSizeB+1.))

    if np.any(np.isnan(lambd)):
        raise ValueError('Found NAN in lambd')

    return lambd


_exponentialNtot2Lambda = _exponentialWC2Lambda


def _exponentialReff2Lambda(effectiveRadius):
    """
    Estimate lambda of exponential distribution from effective radius
    """
    if np.any(np.isnan(effectiveRadius)):
        raise ValueError('Found NAN in effectiveRadius')

    lambd = 3. / effectiveRadius

    if np.any(np.isnan(lambd)):
        raise ValueError('Found NAN in lambd')

    return lambd
