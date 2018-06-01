# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

import numpy as np
import scipy.special

# all functions have to accept hydrometeor size as first parameter

from . import mass


def modifiedGamma(sizeCenter, N0, lambd, mu, gamma):
    """classical modified gamma distribution

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    N0 : array_like
        N0 pre-factor
    lambd : float or array_like
        lambda parameter
    mu : float or array_like
        mu parameter
    gamma : float or array_like
        gamma parameter

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """

    N = N0 * sizeCenter**mu * np.exp(-lambd * sizeCenter**gamma)

    return N


def gamma(sizeCenter, N0, lambd, mu):
    """classical gamma distribution

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    N0 : array_like
        N0 pre-factor
    lambd : float or array_like
        lambda parameter
    mu : float or array_like
        mu parameter

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """
    gamma = 1.
    return modifiedGamma(
      sizeCenter=sizeCenter, N0=N0, lambd=lambd, mu=mu, gamma=gamma
      )


def exponential(sizeCenter, N0, lambd):
    """classical exponential distribution

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    N0 : array_like
        N0 pre-factor
    lambd : float or array_like
        lambda parameter

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """
    mu = 0.
    gamma = 1.
    return modifiedGamma(
      sizeCenter=sizeCenter, N0=N0, lambd=lambd, mu=mu, gamma=gamma
      )


def exponentialField(sizeCenter, temperature, lambd):
    """classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T)

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    temperature : array_like
        ambient temperature [K]
    lambd : float or array_like
        lambda parameter

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """

    N0 = _exponentialField(temperature)
    return exponential(sizeCenter, N0, lambd)


def exponentialFieldWC(
  sizeCenter, temperature, waterContent, massSizeA, massSizeB
  ):
    """classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T). Lambda is
    estimated from the WC which requires a power law mass-size
    relation.

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    temperature : array_like
        ambient temperature [K]
    waterContent : array_like
        hydrometeor water content
    massSizeA : float
        pre-factor mass-size power law
    massSizeB : float
        exponent mass-size power law

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """

    N0 = _exponentialField(temperature)
    lambd = _exponentialWC2Lambda(N0, waterContent, massSizeA, massSizeB)
    return exponential(sizeCenter, N0, lambd)


def exponentialFieldReff(sizeCenter, temperature, effectiveRadius):
    """classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T). Lambda is
    estimated from the effective radius.

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    temperature : array_like
        ambient temperature [K]
    effectiveRadius : array_like
        hydrometeor effective radius

    Returns
    -------
    size distribution : array_like
        calculated size distribution
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
    """classical exponential distribution constrained with N0 and LWC.

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    N0 : array_like
        N0 pre-factor
    waterContent : array_like
        hydrometeor water content
    massSizeA : float
        pre-factor mass-size power law
        (Default value = mass.powerLawLiquidPrefactor)
    massSizeB : float
        exponent mass-size power law
        (Default value = mass.powerLawLiquidExponent)

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """

    lambd = _exponentialWC2Lambda(N0, waterContent, massSizeA, massSizeB)
    return exponential(sizeCenter, N0, lambd)


########################################################################
# Helper functions
########################################################################


def _exponentialField(temperature):
    """N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T)

    Parameters
    ----------
    temperature : array_like
        ambient temperature [K]

    Returns
    -------
    size distribution : array_like
        calculated size distribution
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
    """Estimate lambda of exponential distribution from WC and Ntot

    Parameters
    ----------
    N0 : array_like
        N0 pre-factor
    WC : array_like
        hydrometeor water content
    massSizeA : float
        pre-factor mass-size power law
    massSizeB : float
        exponent mass-size power law

    Returns
    -------
    size distribution : array_like
        calculated size distribution
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
    """Estimate lambda of exponential distribution from effective radius

    Parameters
    ----------
    effectiveRadius : array_like
        hydrometeor effective radius

    Returns
    -------
    size distribution : array_like
        calculated size distribution
    """
    if np.any(np.isnan(effectiveRadius)):
        raise ValueError('Found NAN in effectiveRadius')

    lambd = 3. / effectiveRadius

    if np.any(np.isnan(lambd)):
        raise ValueError('Found NAN in lambd')

    return lambd
