# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot)
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)

from __future__ import absolute_import, division, print_function
import numpy as np

# from .due import due, Doi

import meteo_si.constants
import meteo_si.density
import meteo_si.temperature


__all__ = ['a2e', 'e2a', "e2q", "q2e", "rh2q", "rh2a", "rh_to_iwv",
           "e_sat_gg_ice", "e_sat_gg_water", "q2rh", "a2rh"]


def a2e(a, T):
    """
    Calculate water vapor pressure from the absolute humidity and air
    temperature.

    Parameters
    ----------
    a:
        absolute humidity [kg / m3]
    T:
        Temperature in K

    Returns
    -------

    float :
        vapor pressure [Pa]

    """
    e = a * T * meteo_si.constants.Rvapor

    return e


def e2a(e, T):
    """
    Calculate the absolute humidity from water vapor pressure and air
    temperature.

    Parameters
    ----------
    e:
        vapor pressure [Pa]
    T:
        Temperature in K

    Returns
    -------

    float :
        absolute humidity [kg / m3]

    """

    a = e / (T * meteo_si.constants.Rvapor)

    return a


def e_sat_gg_water(T):
    """
    Calculates the saturation pressure over water after "Guide to
    Meteorological Instruments and Methods of Observation" (CIMO Guide)
    (WMO, 2008).

    Parameters
    ----------
    T:
        Temperature in K

    Returns
    -------

    float :
        saturation pressure [Pa]

    """
    T = meteo_si.temperature.kelvin_2_celsius(T)
    e_sat_gg_water = 100 * 6.112 * np.exp(17.62 * T / (243.12 + T))
    return e_sat_gg_water


def e_sat_gg_ice(T):
    """
    Calculates the saturation pressure over water after "Guide to
    Meteorological Instruments and Methods of Observation" (CIMO Guide)
    (WMO, 2008).

    Parameters
    ----------
    T:
        Temperature in K

    Returns
    -------

    float :
        saturation pressure [Pa]

    """
    T = meteo_si.temperature.kelvin_2_celsius(T)
    e_sat_gg_ice = 100 * 6.112 * np.exp(22.46 * T / (272.62 + T))
    return e_sat_gg_ice


def e2q(e, p):
    """
    Calculate the specific humidity from vapor pressure and air
    pressure.

    Parameters
    ----------
    e:
        vapor pressure [Pa]
    p:
        pressure [Pa]

    Returns
    -------

    float :
        specific humidity [kg / kg]

    """
    q = meteo_si.constants.Mwml * e / (p - (1 - meteo_si.constants.Mwml) * e)
    return q


def q2e(q, p):
    """
    Calculate water vapor pressure from the specific humidity and air
    pressure.

    Parameters
    ----------
    q:
        specific humidity [kg / kg]
    p:
        pressure [Pa]

    Returns
    -------

    float :
        vapor pressure [Pa]

    """

    e = p / ((meteo_si.constants.Mwml / q)+1-meteo_si.constants.Mwml)
    return e


def rh2q(rh, T, p, e_sat_func=e_sat_gg_water):
    """
    Calculate the specific humidity from relative humidity, air temperature,
    and pressure.

    Parameters
    ----------
    rh:
        Relative humidity in Pa / Pa
    T:
        Temperature in K
    p:
        pressure [Pa]
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        specific humidity [kg / kg]

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        if np.any(rh > 5):
            raise TypeError("rh must not be in %")

    eStar = e_sat_func(T)
    e = rh*eStar
    q = e2q(e, p)
    del e, eStar
    return q


def rh2a(rh, T, e_sat_func=e_sat_gg_water):
    """
    Calculate the absolute humidity from relative humidity, air temperature,
    and pressure.

    Parameters
    ----------
    rh:
        Relative humidity in Pa / Pa
    T:
        Temperature in K
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        absolute humidity [kg / m3]

    """

    with np.errstate(divide='ignore', invalid='ignore'):
        if np.any(rh > 5):
            raise TypeError("rh must not be in %")

    e = rh*e_sat_func(T)
    a = e / (meteo_si.constants.Rvapor*T)
    return a


def a2rh(a, T, e_sat_func=e_sat_gg_water):
    """
    Calculate the absolute humidity from relative humidity, air temperature,
    and pressure. Source: Kraus, 'Die Atmosphäre der Erde', Chapter 8.1.2

    Parameters
    ----------
    a:
        absolute humidity [kg / m3]
    T:
        Temperature in K
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        relative humidity [kg / kg]

    """

    e = a*(meteo_si.constants.Rvapor * T)
    rh = e / e_sat_func(T)
    return rh


def q2rh(q, T, p, e_sat_func=e_sat_gg_water):
    """
    Calculate relative humidity from specific humidity. Source: Kraus, 'Die
    Atmosphäre der Erde', Chapter 8.1.2

    Parameters
    ----------
    q:
        specific humidity [kg / kg]
    T:
        Temperature in K
    p:
        pressure [Pa]
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        relative humidity [kg / kg]

    """

    e = p / (meteo_si.constants.Mwml *
             ((1 / q)+(1 / (meteo_si.constants.Mwml) - 1)))

    eStar = e_sat_func(T)
    rh = e / eStar
    return rh


def rh_to_iwv(relhum_lev, temp_lev, press_lev, hgt_lev):
    """
    Integrate relative humidity to obtain the integrated water vapor (IWV)
    column.

    Parameters
    ----------
    relhum_lev:
        relative humidity at levels humidity [Pa / Pa]
    temp_lev:
        Temperature at levels [K]
    press_lev:
        pressure at levels [Pa]
    hgt_levels:
        altitude of levels [m]

    Returns
    -------

    float :
        IWV [kg / m^2]

    """

    dz = np.diff(hgt_lev, axis=-1)
    relhum = (relhum_lev[..., 0:-1] + relhum_lev[..., 1:]) / 2.
    temp = (temp_lev[..., 0:-1] + temp_lev[..., 1:]) / 2.

    xp = -1.*np.log(press_lev[..., 1:] / press_lev[..., 0:-1]) / dz
    press = -1.*press_lev[..., 0:-1] / xp*(np.exp(-xp*dz)-1.) / dz

    q = rh2q(relhum, temp, press)
    rho_moist = meteo_si.density.moist_rho_q(press, temp, q)

    return np.sum(q*rho_moist*dz)
