# -*- coding: utf-8 -*-
import numpy as np

from .. import constants

#all functions have to accept hydrometeor size as first parameter. 
# NOT Soft sphere function must also accept  aspectRatio anddensity)

powerLawLiquidPrefactor = np.pi/6 * constants.rhoWater
powerLawIcePrefactor = np.pi/6 * constants.rhoIce

powerLawLiquidExponent= 3.
powerLawIceExponent= 3.



def powerLaw(sizeCenter,massSizeA,massSizeB):
  """
  classical mass size relation as power law

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

  m : array
    particle mass
  """

  m = massSizeA*sizeCenter**massSizeB

  return m



def ellipsoid(sizeCenter,aspectRatio,density):
  """
  classical mass size relation as power law

  Parameters
  ----------
  sizeCenter : array_like
    particle size at center of size bin
  density : array_like
    fixed particle density

  Returns
  -------

  m : array
    particle mass
  """

  if np.any(aspectRatio!=1):
    raise NotImplementedError('aspectRatio!=1 not implemented yet. Patch this function.')

  massSizeA = np.pi/6 * density
  massSizeB = 3.

  return powerLaw(sizeCenter,massSizeA,massSizeB)


