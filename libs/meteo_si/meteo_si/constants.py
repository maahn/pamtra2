# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

'''
Constants relevant for atmospheric sciences.
'''

Grav = 9.80665  # m/s^2 der Wert fuer mittlere Breiten
Rair = 287.04  # J/kg/K
Rvapor = 461.5  # J/kg/K
Cp = 1005.0  # J/kg/K specific heat capacity
Gamma = -Grav/Cp  # =-0.0097..K/m  adiabatic temperature gradient
Lv = 2.5e6  # J/kg  bei 0C Lv heat of vaporization
Mwml = 0.622  # dimlos, Molmassenverhaeltnis
Tnull = -273.15  # degC absolute zero
Kadiab = Rair/Cp  # dimensionless adiabatic exponenet
G = 9.80665  # gravitational acceleration
