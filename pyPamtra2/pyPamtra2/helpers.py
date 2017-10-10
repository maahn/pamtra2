# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function
from builtins import super 

from collections import OrderedDict
import numpy as np
import xarray as xr

from . import decorators

@decorators.NDto2DtoND(referenceIn=0,noOfInDimsToKeep=1, convertInputs=[0,1],convertOutputs=[0],verbosity=10)
def rayleigh(diameter,K2,frequency):

  """
  To do: move this routine into separte module together with t-matrix, mie etc code.
  """
  C = 299792458.

  K2 = np.asarray(K2)
  diameter = np.asarray(diameter)

  wavelength = C / (frequency*1e9)  
  prefactor = np.pi**5 * K2 / wavelength**4
  back_spec =  prefactor[:,np.newaxis] * diameter**6

  return back_spec


