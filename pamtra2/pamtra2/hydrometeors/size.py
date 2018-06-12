# -*- coding: utf-8 -*-
import numpy as np

# input names are not arbritrary and have to follow Pamtra2 defaults!

# if this is too slow think about implementing @vectorize
# https://numba.pydata.org/numba-doc/dev/user/vectorize.html


def linspaceBounds(nBins=None, Dmin=None, Dmax=None):
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
        particle size at bounds of bin. length is nBins+1
    """

    return np.linspace(Dmin, Dmax, nBins+1)


def logspaceBounds(nBins=None, Dmin=None, Dmax=None):
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
        particle size at bounds of bin. length is nBins+1
    """

    return np.logspace(np.log10(Dmin), np.log10(Dmax), nBins+1)


def boundsToMid(sizeBounds):
    """
    Return mid points of size bounds

    Parameters
    ----------
    sizeBounds : array_like (len nBins+1)
        size boundaries

    Returns
    -------
    sizeCenter : array_like
        particle size at center bin. length is nBins
    """
    return sizeBounds[..., :-1] + 0.5*np.diff(sizeBounds)


def boundsWidth(sizeBounds):
    """
    Return width of size bounds

    Parameters
    ----------
    sizeBounds : array_like (len nBins+1)
        size boundaries

    Returns
    -------
    boundsWidth : array_like
        width of size bounds
    """
    return np.diff(sizeBounds)

