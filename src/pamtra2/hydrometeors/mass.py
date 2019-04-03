# -*- coding: utf-8 -*-
import numpy as np

from .. import constants

# input names are not arbritrary and have to follow Pamtra2 defaults!

# if this is too slow think about implementing @vectorize
# https://numba.pydata.org/numba-doc/dev/user/vectorize.html

powerLawLiquidPrefactor = np.pi/6 * constants.rhoWater
powerLawIcePrefactor = np.pi/6 * constants.rhoIce

powerLawLiquidExponent = 3.
powerLawIceExponent = 3.


def waterSphere(sizeCenter):
    """mass for water spheres

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin

    Returns
    -------
    mass : array_like
        particle mass
    """
    return powerLaw(
        sizeCenter, powerLawLiquidPrefactor, powerLawLiquidExponent)


def iceSphere(sizeCenter):
    """mass for ice spheres

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin

    Returns
    -------
    mass : array_like
        particle mass
    """
    return powerLaw(
        sizeCenter, powerLawIcePrefactor, powerLawIceExponent)


def powerLaw(sizeCenter, massSizeA, massSizeB):
    """classical mass size relation as power law

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    massSizeA : array_like
        mass size pre factor
    massSizeB : float or array_like
        mass size exponent

    Returns
    -------
    mass : array_like
        particle mass
    """

    m = massSizeA*sizeCenter**massSizeB

    return m


def ellipsoid(sizeCenter, aspectRatio, density):
    """mass of a fixed-density ellipsoid

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    density : array_like
        fixed particle density
    aspectRatio : array_like
        particle aspect ratio

    Returns
    -------
    mass : array_like
        particle mass
    """

    if np.all(aspectRatio <= 1):
        mass = oblateEllipsoid(
            sizeCenter, aspectRatio, density)
    elif np.all(aspectRatio > 1):
        mass = prolateEllipsoid(
            sizeCenter, aspectRatio, density)
    else:
        mass = oblateEllipsoid(
            sizeCenter, aspectRatio, density)
        mass[aspectRatio > 1] = prolateEllipsoid(
            sizeCenter, aspectRatio, density)[aspectRatio > 1]
    return mass


def oblateEllipsoid(sizeCenter, aspectRatio, density):
    """mass of a fixed-density oblate ellipsoid

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    density : array_like
        fixed particle density
    aspectRatio : array_like
        particle aspect ratio <1

    Returns
    -------
    mass : array_like
        particle mass
    """

    if np.any(aspectRatio > 1):
        raise ValueError(
            'oblate ellipsoids only')

    A = B = sizeCenter
    C = sizeCenter * aspectRatio
    volume = np.pi/6 * A * B * C

    return volume * density


def prolateEllipsoid(sizeCenter, aspectRatio, density):
    """mass of a fixed-density prolate ellipsoid

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    density : array_like
        fixed particle density
    aspectRatio : array_like
        particle aspect ratio > 1

    Returns
    -------
    mass : array_like
        particle mass
    """

    if np.any(aspectRatio < 1):
        raise ValueError(
            'prolate ellipsoids only')

    A = B = sizeCenter / aspectRatio
    C = sizeCenter
    volume = np.pi/6 * A * B * C

    return volume * density


