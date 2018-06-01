# -*- coding: utf-8 -*-
import numpy as np
from .. import constants

# All functions must accept sizeCenter,aspectRatio,mass as input.


def softEllipsoid(sizeCenter, aspectRatio, mass, minDensity=100,
                  maxDensity=constants.rhoIce):
    """
    oblate (AspectRatio <1) or prolate (asectRatio >0) soft spheres
    """

    if np.all(aspectRatio <= 1):
        density = softOblateEllipsoid(sizeCenter, aspectRatio, mass, minDensity=minDensity,
                                      maxDensity=maxDensity)
    elif np.all(aspectRatio > 1):
        density = softProlateEllipsoid(sizeCenter, mass, minDensity=minDensity,
                                       maxDensity=maxDensity)
    else:
        density = softOblateEllipsoid(sizeCenter, aspectRatio, mass, minDensity=minDensity,
                                      maxDensity=maxDensity)
        density[aspectRatio > 1] = softProlateEllipsoid(sizeCenter, mass,
            minDensity=minDensity, maxDensity=maxDensity)[aspectRatio > 1]

    return density


def softOblateEllipsoid(sizeCenter, aspectRatio, mass, minDensity=100,
                        maxDensity=constants.rhoIce):
    """
    oblate (AspectRatio <1)  soft spheres
    """
    density = (6. * mass) / (np.pi * sizeCenter**3. * aspectRatio)

    density[density < minDensity] = minDensity
    density[density > maxDensity] = maxDensity

    return density

def softProlateEllipsoid(sizeCenter, aspectRatio, mass, minDensity=100,
                         maxDensity=constants.rhoIce):
    """
    prolate (asectRatio >0) soft spheres
    """
    density = (6. * mass * aspectRatio** 2.) / (np.pi * sizeCenter**3.)

    density[density < minDensity] = minDensity
    density[density > maxDensity] = maxDensity

    return density



water = constants.rhoWater
ice = constants.rhoIce

# # -*- coding: utf-8 -*-
# import numpy as np
# from .. import constants

# def softSphere(sizeCenter,aspectRatio,mass=None,minDensity=100,maxDensity=constants.rhoIce):
#   """
#   oblate (AspectRatio <1) or prolate (asectRatio >0) soft spheres
#   """
  
#   massArr = np.zeros(sizeCenter.shape) * np.nan
#   if callable(mass[0]):
#       print('callable')
#       func, kwargs = mass
#       args = [sizeCenter]
#       massArr[:] = func(*args, **kwargs)
#   else:
#       print('not callable',mass)
#       massArr[:] = mass

#   densityArr = np.zeros(sizeCenter.shape) * np.nan
#   densityArr[aspectRatio<=1] = (6. * massArr[aspectRatio<=1]) / (np.pi *  sizeCenter[aspectRatio<=1]**3. * aspectRatio[aspectRatio<=1])
#   densityArr[aspectRatio>1] = (6. * massArr[aspectRatio>1] * aspectRatio[aspectRatio>1]**2.) / (np.pi *  sizeCenter[aspectRatio>1]**3.)

#   densityArr[densityArr<minDensity] = minDensity
#   densityArr[densityArr>maxDensity] = maxDensity

#   return densityArr, massArr


# def solidSphere(sizeCenter,aspectRatio, density=None):

#   if not np.any(aspectRatio!=1):
#     raise NotImplementedError('aspectRatio != 1 not implemented yet!')

#   densityArr = np.zeros(sizeCenter.shape) * np.nan
#   densityArr[:] = density

#   massSizeA = np.pi/6 * densityArr
#   massSizeB = 3.

#   massArr =  massSizePowerLaw(sizeCenter,massSizeA,massSizeB)

#   return densityArr, massArr


# def massSizePowerLaw(sizeCenter,massSizeA=None,massSizeB=None):
#   """
#   classical mass size relation as power law

#   Parameters
#   ----------
#   sizeCenter : array_like
#     particle size at center of size bin
#   massSizeA : array_like
#     mass size pre factor
#   massSizeB : float or array_like
#     mass size exponent

#   Returns
#   -------

#   m : array
#     particle mass
#   """

#   m = massSizeA*sizeCenter**massSizeB

#   return m

