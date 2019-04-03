# -*- coding: utf-8 -
# (c) M. Maahn, 2018

import numpy as np
from pamtra2.libs import meteo_si

from . import pamgasabs_lib

__version__ = '0.1'


def calculate_gas_absorption_rosenkranz98(
  frequency, temperature, waterVaporPressure, pressure, sumResults=False,
  verbosity=0):
    """
    Microwave Gas absorption accoding to Rosenkranz 1998.

    Parameters
    ----------
    frequency : array_like
        frequency [Hz]
    tempk : array_like
        temperature [K]
    waterVaporPressure : array_like
        water vapor density [Pa]
    pres : array_like
        pressure [Pa]
    sumResults : bool, optional
        Sum absAir and absWv when returning (default False)
    verbosity : int, optional
        Verbose level of Fortran module (default 0)

    Returns
    -------
    absAir : array_like
        extinction by dry air [Np/m]
    absWv : array_like
        extinction by water vapor [Np/m]

    """

    pamgasabs_lib.report_module.verbose = verbosity

    # to GHz
    frequencyGHz = frequency/1e9

    absoluteHumidity = meteo_si.humidity.e2a(waterVaporPressure, temperature)

    error, absair, abswv = pamgasabs_lib.rosen98_gasabs(
        frequencyGHz, temperature, absoluteHumidity, pressure)

    # to m
    absair = absair / 1000.
    abswv = abswv / 1000.

    if error > 0:
        raise RuntimeError('Error in Fortran routine rosen98_gasabs')

    if sumResults:
        return absair + abswv
    else:
        return absair, abswv


def calculate_gas_absorption_liebe93(
  frequency, temperature, waterVaporPressure, pressure, verbosity=0):
    """
    Microwave Gas absorption accoding to Rosenkranz 1998.

    Parameters
    ----------
    frequency : array_like
        frequency [Hz]
    tempk : array_like
        temperature [K]
    waterVaporPressure : array_like
        water vapor pressure [Pa]
    pres : array_like
        pressure [Pa]
    verbosity : int, optional
        Verbose level of Fortran module (default 0)

    Returns
    -------
    atmoAbs : array_like
        extinction by moist air [Np/m]
    """

    pamgasabs_lib.report_module.verbose = verbosity

    # to GHz
    frequencyGHz = frequency/1e9
    # to kPa
    pressurekPa = pressure/1000.
    # to Celsius
    temperatureC = _kelvin2Celsius(temperature)
    # to kPa
    waterVaporPressurekPa = waterVaporPressure/1000.
    # we do hydrometeor attenaution somewhere else
    liquidWaterContent = np.zeros_like(temperature)

    error, atmoAbs = pamgasabs_lib.mpm93(
        frequencyGHz,
        pressurekPa,
        waterVaporPressurekPa,
        temperatureC,
        liquidWaterContent,
        )

    # per km to tper m
    atmoAbs = atmoAbs / 1000.

    if error > 0:
        raise RuntimeError('Error in Fortran routine mpm93')

    return atmoAbs


def _kelvin2Celsius(kelvin):
    tFreezing = 273.15
    return kelvin - tFreezing
