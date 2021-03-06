# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot)
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)

from __future__ import absolute_import, division, print_function

import numpy as np

from . import constants, humidity

# from .due import due, Doi

# from .temperature import *


__all__ = ["moist_rho_rh", "moist_rho_q"]


def moist_rho_rh(p, T, rh, qm=0):
    """
    Compute the arithmetic circular mean, ignoring NaNs.

    Parameters
    ----------
    p :
        Pressure in Pa
    T:
        Temperature in K
    rh:
        Relative humidity in Pa/Pa
    qm: optional
        sum of mixing ratios in kg/kg of other species which contribute to
        the air mass! (ice, snow, cloud etc.)

    Returns
    -------

    float :
        Density of moist air [kg/m^3]

    Example
    -------
    moist_rho_rh(p,T,rh,q_ice,q_snow,q_rain,q_cloud,q_graupel,q_hail)

    """
    q = humidity.rh2q(rh, T, p)

    return moist_rho_q(p, T, q, qm)


def moist_rho_q(p, T, q, qm=0):
    """
    Compute the arithmetic circular mean, ignoring NaNs.

    Parameters
    ----------
    p :
        Pressure in Pa
    T:
        Temperature in K
    q:
        specific humidity in kg/kg
    qm: optional
        sum of mixing ratios in kg/kg of other species which contribute to
        the air mass! (ice, snow, cloud etc.)

    Returns
    -------

    float :
        Density of moist air [kg/m^3]

    Example
    -------
    moist_rho_q(p,T,q,q_ice,q_snow,q_rain,q_cloud,q_graupel,q_hail)

    """

    moist_rho_q = p / (constants.Rair * T *
                       (1 + (constants.Rvapor /
                        constants.Rair-1) * q - qm))

    return moist_rho_q
