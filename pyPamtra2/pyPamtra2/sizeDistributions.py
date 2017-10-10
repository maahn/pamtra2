# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

import numpy as np
import xarray as xr


#all functions have to accept the pyPamtra2 object and hydroIndex as first and second parameter  

# to do: add more options


def exponential(pam,hydroIndex,N0,lambd):
  """
  classical exponential distribution
  """
  diameters = pam.data.hydroSize.values[...,hydroIndex,:]

  return N0 * np.exp(lambd * diameters)


