# -*- coding: utf-8 -*-
import numpy as np

# from numba import vectorize, float64, float32
from .. import constants


# @vectorize([float32(float32, float32, float32, float32, float32),
#             float64(float64, float64, float64, float64, float64)])
def heymsfield10_particles(
    sizeCenter,
    mass,
    crossSectionArea,
    dynamicViscosity,
    dryAirDensity,
):

    k = 0.5  # defined in the paper
    delta_0 = 8.0
    C_0 = 0.35
    g = constants.gravitation

    area_proj = crossSectionArea/((np.pi/4.)*sizeCenter**2)

    # eq 9
    Xstar = 8.*dryAirDensity*mass*g / \
        (np.pi*area_proj**(1. - k)*dynamicViscosity**2)
    # eq10
    Re = delta_0**2/4. * \
        ((1. + ((4.*np.sqrt(Xstar))/(delta_0**2*np.sqrt(C_0))))**0.5-1)**2

    velSpec = dynamicViscosity*Re/(dryAirDensity*sizeCenter)
    return velSpec


# @vectorize([float32(float32, float32, float32),
#             float64(float64, float64, float64)])
def khvorostyanov01_drops(
    sizeCenter,
    dryAirDensity,
    kinematicViscosity,
):

    # variables to cgs...
    rho_water_cp = constants.rhoWater/1000.  # g/cm³
    g_cp = constants.gravitation*100.  # cm/s
    dryAirDensity_cp = dryAirDensity/1000.  # g/cm³
    sizeCenter_cp = 100.*sizeCenter  # cm
    kinematicViscosity_cp = kinematicViscosity*10000.  # cm² /s
    c1 = 0.0902
    delta0 = 9.06

    lam = 0.47
    # eq. 3.4
    xi = np.exp(-sizeCenter_cp/lam) + (1. - np.exp(-sizeCenter_cp/lam)
                                       )*(1./(1. + (sizeCenter_cp/lam)))

    vB = 4/3.*np.pi*(sizeCenter_cp*0.5)**3*xi  # correct for non-spherity
    A = np.pi*(sizeCenter_cp*0.5)**2  # no correction neccessary for A

    X = 2*vB*(rho_water_cp - dryAirDensity_cp)*g_cp*sizeCenter_cp**2 / \
        (A*dryAirDensity_cp*kinematicViscosity_cp**2)  # eq. 2.7

    bRe = 0.5*c1*X**0.5*((1. + c1*X**0.5)**0.5 - 1.)**(-1) * \
        (1 + c1*X**0.5)**(-0.5)  # eq. 2.12

    aRe = (delta0**2/4.)*((1. + c1*X**0.5)**0.5 - 1.)**2/X**bRe  # eq. 2.13

    velSpec = aRe*kinematicViscosity_cp**(1 - 2*bRe)*(4./3.*g_cp*xi*(
        (rho_water_cp - dryAirDensity_cp)/dryAirDensity_cp))**bRe * \
        sizeCenter_cp**(3*bRe - 1)  # eq. 2.20c

    velSpec = velSpec/100.  # CGS to SI, now m/s

    return velSpec
