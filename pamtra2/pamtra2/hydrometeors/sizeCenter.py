# -*- coding: utf-8 -*-
import numpy as np

# input names are not arbritrary and have to follow Pamtra2 defaults!

# if this is too slow think about implementing @vectorize
# https://numba.pydata.org/numba-doc/dev/user/vectorize.html


def linspace(nBins, Dmin=None, Dmax=None):
    """
    Distribute size bins linearly.

    Parameters
    ----------
    nBins : int
        number of bins
    Dmin :
         minimum size (Default value = None)
    Dmax :
         maximum size (Default value = None)

    Returns
    -------
    size : array_like
        particle size at center of bin.
    """

    return np.linspace(Dmin, Dmax, nBins)


def logspace(nBins, Dmin=None, Dmax=None):
    """
    Distribute size bins in log scale.

    Parameters
    ----------
    nBins : int
        number of bins
    Dmin :
         minimum size (Default value = None)
    Dmax :
         maximum size (Default value = None)

    Returns
    -------
    size : array_like
        particle size at center of bin.
    """

    return np.logspace(np.log10(Dmin), np.log10(Dmax), nBins)
