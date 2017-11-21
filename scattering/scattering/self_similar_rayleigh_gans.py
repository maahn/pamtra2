# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 18:41:17 2017

@author: dori
"""

import numpy as np

c = 2.99792458e8

#INPUT
diameters = np.linspace(0.001,0.025,50)
frequency = 94.0e9
wavelength = c/frequency
volume = 1.0e-9

# CONSTANTS FOR MY PARTICLES
kappa = 0.190031
beta  =  0.030681461
gamma = 1.3002167
zeta1 = 0.29466184

# DEFAULT VALUES #
verbose = 0
T       = 268
n       = complex(1.774618673, 0.001197613076)
eps = n*n.conjugate()

N_diameters = len(diameters)

if (N_diameters > 1):
    diameters = diameters[0]
    N_diameters = 1
    print('Currently working on scalars ...  I am just translating IDL')

### Now there is a fantastic IDL routine that check if parameter values are finite

K = (eps-1)/(eps+2)
K2 = K*K.conjugate()

#wavenumber
k = 2*np.pi/wavelength

#size parameter
x = k*diameters

#Factor resulting from Eq. 1 and Eq. 4 (all elements except PHI_SSRGA(x)) in Hogan et al., 2016 
prefactor = 9.0*np.pi * k**4. * K2 * volume**2. / 16.

term1 = cos(x) * ( (1.0+kappa/3.0) * (1.0/(2.0*x+np.pi) - 1.0/(2.0*x-np.pi)) - kappa * (1.0/(2.0*x+3.0*np.pi) - 1.0/(2.0*x-3.0*np.pi)) )
term1 = term1**2.


# Initialize scattering variables
c_bsc = 0d	#backscattering cross section [m2]
c_abs = 0d      #absorption cross section [m2]
c_sca = 0d      #scattering cross section [m2]

#--- define scattering angles for phase function (we could also make it a keyword in the future)
theta = (FINDGEN(200.+1.)) * (180./200.)
theta_rad = theta * !dtor

n_theta = n_elements(theta)
d_theta_rad = theta_rad[1] - theta_rad[0]

ph_func = DBLARR(n_theta)


#**** ABSORPTION:
c_abs = 3. * k * volume * IMAGINARY(-1.*Kxyz)	#Hogan et al., 2016, Eq. 8
  

#**** BACKSCATTERING:
#Initialize the summation in the second term in the braces of Eq. 12
sum = 0.

#-- Decide how many terms are needed
jmax = floor(5.*x/!dpi + 1.0)

#-- Evaluate summation
for j=1, jmax do begin
  if j eq 1 then zeta=zeta1 else zeta=1.
  sum = sum + zeta * (2.0*j)^(-1.*gamma) * sin(x)^2. * ( (1.0/(2.0*x+2.0*!dpi*j)^2.) + (1.0/(2.0*x-2.0*!dpi*j)^2.) )
endfor

#-- Put the terms together
c_bsc = prefactor*(term1 + beta*sum)


#**** SCATTERING PHASE FUNCTION AND SCATTERING CROSS SECTION
for i_th=0, n_theta-1 do begin
  
  new_x = x*sin(theta_rad[i_th]/2.)
  #print, 'new, old x', new_x, x
  
  
  #First term in the braces of Eq. 4 representing the mean structure
  new_term1 = cos(new_x) * ( (1.0+kappa/3.0) * (1.0/(2.0*new_x+!dpi) - 1.0/(2.0*new_x-!dpi)) - $
			kappa * (1.0/(2.0*new_x+3.0*!dpi) - 1.0/(2.0*new_x-3.0*!dpi)) )
  new_term1 = new_term1^2.

  #Initialize the summation in the second term in the braces of Eq. 12
  new_sum = 0.

  #Decide how many terms are needed
  jmax = floor(5.*new_x/!dpi + 1.0)

  #Evaluate summation
  for j=1, jmax do begin
    if j eq 1 then zeta=zeta1 else zeta=1.
    new_sum = new_sum + zeta * (2.0*j)^(-1.*gamma) * sin(new_x)^2. * ( (1.0/(2.0*new_x+2.0*!dpi*j)^2.) + (1.0/(2.0*new_x-2.0*!dpi*j)^2.) )
  endfor
  
  ph_func[i_th] = prefactor * ((1. + (cos(theta_rad[i_th]))^2.)/2.) * (new_term1 + beta*new_sum)
  
  c_sca = c_sca + (0.5 * ph_func[i_th] * sin(theta_rad[i_th]) * d_theta_rad)	#formula 10
  
  
#     print, d_theta
#     print, sin(!dtor*theta[i_th])
#     print,  ph_func[i_th]
#    print, c_sca
#     
#    aa=''
#    read,aa

# print, 'New x: ', new_x
  
endfor


#-- normalize the phase function with c_sca
ph_func[*] = ph_func[*] / REPLICATE(2.*c_sca, n_theta)	#normalization with 2*sigma_sca is convention used by Janni

#-- calculate asymmetry factor
asym = 0.
for i_th=0, n_theta-1 do begin
  asym = asym +  (ph_func[i_th] * sin(theta_rad[i_th]) *cos(theta_rad[i_th]) * d_theta_rad)
endfor

#-- check if integral over normalized phase function is indeed one
ph_func_int = 0.
for i_th=0, n_theta-1 do begin
  ph_func_int = ph_func_int +  (ph_func[i_th] * sin(theta_rad[i_th])) * d_theta_rad
endfor

#print, 'int(ph_func): ', ph_func_int
 

if verbose eq 1 then print, 'Scattering phase function(theta): ', ph_func
if verbose eq 1 then print, 'Backscattering Coeff (m2): ', c_bsc
if verbose eq 1 then print, 'Scattering Coeff (m2): ', c_sca
if verbose eq 1 then print, 'Absorption Coeff (m2): ', c_abs