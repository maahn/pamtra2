# -*- coding: utf-8 -
# (c) M. Maahn, 2018

from . import pygasabs_lib

__version__ = '0.1'


def calculate_gas_absorption(
  frequency, temperature, absoluteHumidity, pressure):
    """

    Parameters
    ----------
    frequency : array_like
        frequency [Hz]
    tempk : array_like
        temperature [K]
    absoluteHumidity : array_like
        water vapor density [kg/m**3]
    pres : array_like
        pressure [Pa]

    Returns
    -------
   absAir : array_like
        extinction by dry air [Np/m]
   absWv : array_like
        extinction by water vapor [Np/m]

    """
    # to GHz
    frequency = frequency/1e9

    error, absair, abswv = pygasabs_lib.rosen98_gasabs(
        frequency, temperature, absoluteHumidity, pressure)

    # to m
    absair = absair / 1000.
    abswv = abswv / 1000.

    if error > 0:
        raise RuntimeError('Error in Fortran routine calculate_gas_absorption')

    return absair, abswv
