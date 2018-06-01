# -*- coding: utf-8 -*-
import numpy as np
from copy import deepcopy

from . import decorators


@decorators.NDto2DtoND(
  referenceIn=0, noOfInDimsToKeep=1, convertInputs=[0, 1],
  convertOutputs=[0], verbosity=10
  )
def rayleigh(diameter, K2, frequency):
    """
    To do: move this routine into separte module together with t-matrix,
    mie etc code.
    """
    C = 299792458.

    K2 = np.asarray(K2)
    diameter = np.asarray(diameter)

    wavelength = C / (frequency*1e9)
    prefactor = np.pi**5 * K2 / wavelength**4
    back_spec = prefactor[:, np.newaxis] * diameter**6

    return back_spec


def concatDicts(*dicts):
    """
    concatenate (ordered) dicts
    """
    dMerged = deepcopy(dicts[0])
    for dd in range(1, len(dicts)):
        dMerged.update(dicts[dd])
    return dMerged
