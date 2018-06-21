# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np

# from .constants import *

'''
Functions to deal with wind observations.

'''


def circular_mean(angles):
    """
    Compute the arithmetic circular mean, not ignoring NaNs.

    Parameters
    ----------
    angles : list or array
       The angles for averaging in radians.

    Returns
    -------

    mean : float
        The circular mean in radians.

    """
    if np.any(np.isnan(angles)):
        return np.nan
    else:
        return nan_circular_mean(angles)


def nan_circular_mean(angles):
    """
    Compute the arithmetic circular mean, ignoring NaNs.

    Parameters
    ----------
    angles : list or array
       The angles for averaging in radians.

    Returns
    -------

    mean : float
        The circular mean in radians.

    """
    x = np.nansum(np.cos(angles))
    y = np.nansum(np.sin(angles))
    mean = np.arctan2(y, x)
    if mean < 0:
        mean = mean + (np.pi*2)
    return mean


def circular_mean_deg(angles):
    """
    Compute the arithmetic circular mean, not ignoring NaNs.

    Parameters
    ----------
    angles : list or array
       The angles for averaging in degrees.

    Returns
    -------

    mean : float
        The circular mean in degrees.

    """
    if np.any(np.isnan(angles)):
        return np.nan
    else:
        return nan_circular_mean_deg(angles)


def nan_circular_mean_deg(angles):
    """
    Compute the arithmetic circular mean, ignoring NaNs.

    Parameters
    ----------
    angles : list or array
       The angles for averaging in degrees.

    Returns
    -------

    mean : float
        The circular mean in degrees.

    """
    return np.rad2deg(nan_circular_mean(np.deg2rad(angles)))
