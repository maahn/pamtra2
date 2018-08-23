# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import warnings

import numpy as np
import scipy.special
import xarray as xr

from . import mass
from .. import units

# input names are not arbritrary and have to follow Pamtra2 defaults!

# if this is too slow think about implementing @vectorize
# https://numba.pydata.org/numba-doc/dev/user/vectorize.html


def fromSizeDistribution(sizeDistribution, sizeBoundsWidth):
    '''
    convert size distribution [1/m4] to number concentration [1/m3]

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    '''
    return sizeDistribution * sizeBoundsWidth


def monoDisperse(sizeCenter, Ntot, nBins):
    """constant distribution

    Parameters
    ----------
    Ntot : array_like
        total number of particles for the whole size spectruum
    nBins : int
        number of size bins

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """
    if isinstance(sizeCenter, xr.DataArray):
        N = xr.ones_like(sizeCenter)
    else:
        N = np.ones_like(sizeCenter)
    N = (Ntot/nBins) * N
    return N


def monoDisperseWC(hydrometeorContent, mass):
    """constant distribution for a given water content

    Parameters
    ----------
    hydrometeorContent : array_like
        hydrometeor water content [kg/m^3]
    mass : array_like
        mass of hydrometeors [kg]
    nBins : int
        number of size bins
    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """
    if isinstance(mass, xr.DataArray):
        mass1Bin = mass.sum('sizeBin')
        N = xr.ones_like(mass)
    else:
        mass1Bin = mass.sum(-1)
        N = np.ones_like(mass)

    N = (hydrometeorContent / mass1Bin) * N

    return N


def modifiedGamma(sizeCenter, sizeBoundsWidth, N0, lambd, mu, gamma):
    """classical modified gamma distribution

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
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
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """

    N = N0 * sizeCenter**mu * np.exp(-lambd * sizeCenter**gamma)

    return N * sizeBoundsWidth


def gamma(sizeCenter, sizeBoundsWidth, N0, lambd, mu):
    """classical gamma distribution

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
    N0 : array_like
        N0 pre-factor
    lambd : float or array_like
        lambda parameter
    mu : float or array_like
        mu parameter

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """
    gamma = 1.
    return modifiedGamma(
        sizeCenter=sizeCenter, sizeBoundsWidth=sizeBoundsWidth,
        N0=N0, lambd=lambd, mu=mu, gamma=gamma
    )


def exponential(sizeCenter, sizeBoundsWidth, N0, lambd):
    """classical exponential distribution

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
    N0 : array_like
        N0 pre-factor
    lambd : float or array_like
        lambda parameter

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """
    mu = 0.
    gamma = 1.
    return modifiedGamma(
        sizeCenter=sizeCenter, sizeBoundsWidth=sizeBoundsWidth,
        N0=N0, lambd=lambd, mu=mu, gamma=gamma
    )


def exponentialMarshallPalmer(sizeCenter, sizeBoundsWidth, rainRate):
    """classical exponential distribution for rain following Marshall
    Palmer, 1948.

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
    rainRate : array_like
        rain rate in mm/hour!

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """

    N0 = 0.08 * 100**4
    lambd = 41*rainRate**(-0.21)*100
    return exponential(sizeCenter, sizeBoundsWidth, N0, lambd)


def exponentialField(sizeCenter, sizeBoundsWidth, temperature, lambd):
    """classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T)

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
    temperature : array_like
        ambient temperature [K]
    lambd : float or array_like
        lambda parameter

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """

    N0 = _exponentialField(temperature)
    return exponential(sizeCenter, sizeBoundsWidth, N0, lambd)


def exponentialFieldWC(
    sizeCenter, sizeBoundsWidth, temperature, hydrometeorContent, massSizeA,
    massSizeB
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
    sizeBoundsWidth : array_like
        particle size bin width
    temperature : array_like
        ambient temperature [K]
    hydrometeorContent : array_like
        hydrometeor water content [kg/m^3]
    massSizeA : float
        pre-factor mass-size power law
    massSizeB : float
        exponent mass-size power law

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """

    N0 = _exponentialField(temperature)
    lambd = _exponentialWC2Lambda(N0, hydrometeorContent, massSizeA, massSizeB)
    return exponential(sizeCenter, sizeBoundsWidth, N0, lambd)


def exponentialFieldReff(sizeCenter, sizeBoundsWidth, temperature, effectiveRadius):
    """classical exponential distribution. N0 is estimated using
    Field et al. (2005 QJRM, end of page 2008 + end of page 2009
    for the relation between N_0 and N_0,23) n_0 = n_0(T). Lambda is
    estimated from the effective radius.

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
    temperature : array_like
        ambient temperature [K]
    effectiveRadius : array_like
        hydrometeor effective radius

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """
    N0 = _exponentialField(temperature)
    lambd = _exponentialReff2Lambda(effectiveRadius)
    return exponential(sizeCenter, sizeBoundsWidth, N0, lambd)


def exponentialN0WC(
    sizeCenter,
    sizeBoundsWidth,
    N0,
    hydrometeorContent,
    massSizeA=mass.powerLawLiquidPrefactor,
    massSizeB=mass.powerLawLiquidExponent,
):
    """classical exponential distribution constrained with N0 and LWC.

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    sizeBoundsWidth : array_like
        particle size bin width
    N0 : array_like
        N0 pre-factor
    hydrometeorContent : array_like
        hydrometeor water content [kg/m^3]
    massSizeA : float
        pre-factor mass-size power law
        (Default value = mass.powerLawLiquidPrefactor)
    massSizeB : float
        exponent mass-size power law
        (Default value = mass.powerLawLiquidExponent)

    Returns
    -------
    number concentration : array_like
        calculated number concentration (NOT normalized, i.e. unit is 1/m3
        instead of 1/m4)
    """

    lambd = _exponentialWC2Lambda(N0, hydrometeorContent, massSizeA, massSizeB)
    return exponential(sizeCenter, sizeBoundsWidth, N0, lambd)


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
    N0
    """

    if np.any(np.isnan(temperature)):
        raise ValueError('Found NAN in temperature')

    N0 = 7.628e6 * np.exp(-0.107 * units.kelvin2Celsius(temperature))

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
        hydrometeor water content [kg/m3]
    massSizeA : float
        pre-factor mass-size power law
    massSizeB : float
        exponent mass-size power law

    Returns
    -------
    lambda
    """

    warnings.warn('Truncation effect on the PSD are not considered. '
                  'I.e., typically mass is lost!')

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
    lambda
    """

    raise NotImplementedError('See tests, they fail by a factor of 2?')

    if np.any(np.isnan(effectiveRadius)):
        raise ValueError('Found NAN in effectiveRadius')

    lambd = 3. / effectiveRadius

    if np.any(np.isnan(lambd)):
        raise ValueError('Found NAN in lambd')

    return lambd
