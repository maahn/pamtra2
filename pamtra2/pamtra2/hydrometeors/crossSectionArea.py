# -*- coding: utf-8 -*-
import numpy as np

# if this is too slow think about implementing @vectorize
# https://numba.pydata.org/numba-doc/dev/user/vectorize.html


def powerLaw(sizeCenter, areaSizeA, areaSizeB):
    """classical cross section area size relation as power law

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
    area : array_like
        cross section area
    """

    area = areaSizeA*sizeCenter**areaSizeB

    return area


def sphere(sizeCenter):
    """cross section of a spehere is equal to circle area

    Parameters
    ----------
    sizeCenter : array_like
        particle size at center of size bin

    Returns
    -------
    area : array_like
        cross section area
    """

    areaSizeA = np.pi/4.
    areaSizeB = 2.

    return powerLaw(sizeCenter, areaSizeA, areaSizeB)
