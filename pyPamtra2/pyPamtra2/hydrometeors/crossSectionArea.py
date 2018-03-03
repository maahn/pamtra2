# -*- coding: utf-8 -*-
import numpy as np


#all functions have to accept hydrometeor size as first parameter  



def powerLaw(sizeCenter,areaSizeA,areaSizeB):
  """
  classical cross section area size relation as power law

  Parameters
  ----------
  sizeCenter : array_like
    particle size at center of size bin
  areaSizeA : array_like
    area size pre factor
  areaSizeB : float or array_like
    area size exponent

  Returns
  -------

  a : array
    particle area
  """

  area = areaSizeA*sizeCenter**areaSizeB

  return area

def sphere(sizeCenter):
  """
  cross section of a spehere is equal to circle area

  Parameters
  ----------
  sizeCenter : array_like
    particle size at center of size bin

  Returns
  -------

  a : array
    particle area
  """

  areaSizeA = np.pi/4.*sizeCenter**2
  areaSizeB = 2

  return powerLaw(sizeCenter,areaSizeA,areaSizeB)

