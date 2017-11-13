# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



from pytmatrix import tmatrix
from pytmatrix import tmatrix_aux
from pytmatrix import refractive
from pytmatrix import scatter
import matplotlib.pyplot as plt

import numpy as np

import sys
sys.path.append('../')

wl = tmatrix_aux.wl_S
m = refractive.mi(tmatrix_aux.wl_S,refractive.ice_density)

pi5 = np.pi**5.0
wl4 = wl**4.0
K = (m**2+1.0)/(m**2-2)
K2 = (K*K.conj()).real
pref = pi5*K2/wl4

scatt = tmatrix.Scatterer(radius=1.0,
                          wavelength=wl,
                          m=m,
                          axis_ratio=1.0)

plt.figure()
ax = plt.gca()
for th in np.linspace(0.0,180):
    scatt.thet = th
    [[S11,S12],[S21,S22]] = scatt.get_S()
    sig = scatter.sca_xsect(scatt)
    ax.scatter(th,1e6*sig,c='m')
    #ax.scatter(th,S11,c='b')
    #ax.scatter(th,S12,c='r')
    #ax.scatter(th,S21,c='g',marker='x')
    #ax.scatter(th,S22,c='k')
#ax.set_ylim([-0.0015,0.0015])
#ax.set_ylim([-0.001355,-0.00135])

plt.figure()
ax = plt.gca()
for s in np.linspace(1,10,50):
    scatt = tmatrix.Scatterer(radius=s,
                              wavelength=wl,
                              m=m,
                              axis_ratio=1.0)
    R = pref * s**6.0
    ax.scatter(2*np.pi*s/wl,scatter.sca_xsect(scatt),c='b')
    ax.scatter(2*np.pi*s/wl,R,c='g',marker='x')