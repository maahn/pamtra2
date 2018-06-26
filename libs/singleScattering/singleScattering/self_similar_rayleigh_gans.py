# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 18:41:17 2017

@author: dori
"""

import numpy as np
import pandas as pd
from scipy import interp

leinonen_table = pd.read_csv('/home/dori/develop/pyPamtra2/libs/singleScattering/singleScattering/ssrg_coeffs_jussiagg_simult.dat',
                             delim_whitespace=True)


def leinonen_coeff(D, elwp):
    table = leinonen_table[leinonen_table.ELWP == elwp].set_index('D')
    # print(table.columns)
    beta = interp(D, table.index.values, table.beta_z)
    gamma = interp(D, table.index.values, table.gamma_z)
    kappa = interp(D, table.index.values, table.kappa_z)
    zeta1 = interp(D, table.index.values, table.zeta1_z)
    return beta, gamma, kappa, zeta1

c = 2.99792458e8

# INPUT
diameters = np.linspace(0.001, 0.025, 50)

brandes = lambda D: 7.9e-5*D**2.1
smalles = lambda D: 4.1e-5*D**2.5


def backscattering(frequency, diameters, n, table=None, mass=None):  # SI units
    wavelength = c/frequency
    if mass is None:
        mass = min(brandes(diameters*1.0e3), smalles((diameters*1.0e3)))
    volume = mass/917.

    # CONSTANTS FOR MY PARTICLES
    kappa = 0.190031
    beta = 0.030681461
    gamma = 1.3002167
    zeta1 = 0.29466184

    if table == 'leinonen':
        #kappa = 0.189177
        #beta = 3.06939
        #gamma = 2.53192
        #zeta1 = 0.0709529
        beta, gamma, kappa, zeta1 = leinonen_coeff(diameters, 0.0)
        gamma = -gamma
        # print(kappa, beta, gamma, zeta1)
    

    eps = n*n
    K = (eps-1.0)/(eps+2.0)
    K2 = (K*K.conjugate()).real

    k = 2.0*np.pi/wavelength  # wavenumber
    x = k*diameters           # size parameters

    # from Eq. 1 and Eq. 4 (all elements except PHI_SSRGA(x)) in Hogan (2016)
    prefactor = 9.0*np.pi * k**4. * K2 * volume**2. / 16.0

    term1 = np.cos(x)*((1.0+kappa/3.0)*(1.0/(2.0*x+np.pi)-1.0/(2.0*x-np.pi)) -
                       kappa*(1.0/(2.0*x+3.0*np.pi) - 1.0/(2.0*x-3.0*np.pi)))
    term1 = term1**2.

    # Initialize scattering variables
    c_bsc = 0.0  # backscattering cross section [m2]
    c_abs = 0.0  # absorption cross section [m2]
    c_sca = 0.0  # scattering cross section [m2]

    # define scattering angles for phase function
    theta = (np.arange(301.0)) * (180./300.)
    theta_rad = theta * np.pi/180.

    n_theta = len(theta)
    d_theta_rad = theta_rad[1] - theta_rad[0]

    ph_func = np.ndarray(n_theta)

    # ABSORPTION:
    Kxyz = K  # TODO: here I assume polarizability does not change
    c_abs = 3.*k*volume*Kxyz.imag  # Hogan et al., 2016, Eq. 8
    # print(k,volume,Kxyz,Kxyz.imag)
    # BACKSCATTERING:
    # Initialize the summation in the second term in the braces of Eq. 12
    thesum = 0.

    # Decide how many terms are needed
    jmax = int(5.*x/np.pi + 1.0)

    # Evaluate summation
    for j in range(jmax):
        if j == 0:
            zeta = zeta1
        else:
            zeta = 1.
        jj = j + 1.0
        term_zeta = zeta*(2.0*jj)**(-1.0*gamma)*np.sin(x)**2.
        term_x = ((1.0/(2.0*x+2.0*np.pi*jj)**2.) +
                  (1.0/(2.0*x-2.0*np.pi*jj)**2.))
        increment = term_zeta*term_x
        thesum = thesum + increment
    # Put the terms together
    c_bsc = prefactor*(term1 + beta*thesum)
    # print(mass, c_bsc)

    # SCATTERING PHASE FUNCTION AND SCATTERING CROSS SECTION
    for i_th in range(n_theta):
        new_x = x*np.sin(theta_rad[i_th]*0.5)
        # First term in the braces of Eq. 4 representing the mean structure
        new_term1 = np.cos(new_x)*((1.0+kappa/3.0)*(1.0/(2.0*new_x+np.pi) -
                                                    1.0/(2.0*new_x-np.pi)) -
                                   kappa*(1.0/(2.0*new_x+3.0*np.pi) -
                                          1.0/(2.0*new_x-3.0*np.pi)))
        new_term1 = new_term1**2.

        # Initialize the summation in the second term in the braces of Eq. 12
        new_sum = 0.
        # Decide how many terms are needed
        jmax = int(5.*new_x/np.pi + 1.0)

        # Evaluate summation
        for j in range(jmax):
            if j == 0:
                zeta = zeta1
            else:
                zeta = 1.
            term_a = zeta*(2.0*(j+1.0))**(-1.*gamma)*np.sin(new_x)**2.
            term_b = 1.0/(2.0*new_x + 2.0*np.pi*(j+1.0))**2.
            term_c = 1.0/(2.0*new_x - 2.0*np.pi*(j+1.0))**2.
            increment = term_a*(term_b+term_c)
            new_sum = new_sum + increment

        cos2th = (np.cos(theta_rad[i_th]))**2.
        ph_func[i_th] = prefactor*((1. + cos2th)/2.)*(new_term1 + beta*new_sum)

        c_sca = c_sca + (0.5*ph_func[i_th]*np.sin(theta_rad[i_th])*d_theta_rad)

    # normalize the phase function with c_sca
    ph_func = ph_func/(2.*c_sca)  # norm 2*c_sca is convention used by Janni

    # calculate asymmetry parameter
    asym = 0.
    for i_th in range(n_theta):
        sinth = np.sin(theta_rad[i_th])
        costh = np.cos(theta_rad[i_th])
        asym = asym + ph_func[i_th]*sinth*costh*d_theta_rad

    # check if integral over normalized phase function is indeed one
# ph_func_int = 0.
# for i_th=0, n_theta-1 do begin
#   ph_func_int = ph_func_int +  (ph_func[i_th] * sin(theta_rad[i_th])) * d_theta_rad
# endfor
# print, 'int(ph_func): ', ph_func_int
# if verbose eq 1 then print, 'Scattering phase function(theta): ', ph_func
# if verbose eq 1 then print, 'Backscattering Coeff (m2): ', c_bsc
# if verbose eq 1 then print, 'Scattering Coeff (m2): ', c_sca
# if verbose eq 1 then print, 'Absorption Coeff (m2): ', c_abs
    return [c_bsc, c_abs, c_sca, asym, mass]
