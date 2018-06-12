# -*- coding: utf-8 -*-
import numpy as np
from .. import constants

# input names are not arbritrary and have to follow Pamtra2 defaults!


def softEllipsoid(sizeCenter, aspectRatio, mass, minDensity=100,
                  maxDensity=constants.rhoIce):
    """oblate (AspectRatio <1) or prolate (asectRatio >0) soft spheres

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    aspectRatio : array_like
        particle aspect ratio
    mass : array_like
        particle mass
    minDensity : float
         minimal density (Default value = 100)
    maxDensity : float
         maximum density (Default value = constants.rhoIce)

    Returns
    -------
    density : array_like
        density of hydrometeor
    """

    if np.all(aspectRatio <= 1):
        density = softOblateEllipsoid(
            sizeCenter, aspectRatio, mass,
            minDensity=minDensity,
            maxDensity=maxDensity
            )
    elif np.all(aspectRatio > 1):
        density = softProlateEllipsoid(
            sizeCenter, mass, minDensity=minDensity,
            maxDensity=maxDensity
            )
    else:
        density = softOblateEllipsoid(
            sizeCenter, aspectRatio, mass, minDensity=minDensity,
            maxDensity=maxDensity
            )
        density[aspectRatio > 1] = softProlateEllipsoid(
            sizeCenter, mass,
            minDensity=minDensity, maxDensity=maxDensity
            )[aspectRatio > 1]

    return density


def softOblateEllipsoid(sizeCenter, aspectRatio, mass, minDensity=100,
                        maxDensity=constants.rhoIce):
    """oblate (AspectRatio <1)  soft spheres

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    aspectRatio : array_like
        particle aspect ratio
    mass : array_like
        particle mass
    minDensity : float
         minimal density (Default value = 100)
    maxDensity : float
         maximum density (Default value = constants.rhoIce)

    Returns
    -------
    density : array_like
        density of hydrometeor
    """
    density = (6. * mass) / (np.pi * sizeCenter**3. * aspectRatio)

    density[density < minDensity] = minDensity
    density[density > maxDensity] = maxDensity

    return density


def softProlateEllipsoid(sizeCenter, aspectRatio, mass, minDensity=100,
                         maxDensity=constants.rhoIce):
    """prolate (asectRatio >0) soft spheres

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin
    aspectRatio : array_like
        particle aspect ratio
    mass : array_like
        particle mass
    minDensity : float
         minimal density (Default value = 100)
    maxDensity : float
         maximum density (Default value = constants.rhoIce)

    Returns
    -------
    density : array_like
        density of hydrometeor
    """
    density = (6. * mass * aspectRatio ** 2.) / (np.pi * sizeCenter**3.)

    density[density < minDensity] = minDensity
    density[density > maxDensity] = maxDensity

    return density


water = constants.rhoWater
ice = constants.rhoIce
