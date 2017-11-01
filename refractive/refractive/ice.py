import numpy as np
"""
Matzler (2006)

"""

B1 = 0.0207
b = 335.
B2 = 1.16e-11
c = 299792458.

def ice(T,f):
    eps1  = 3.1884+(T-273)*9.1e-4
    theta = 300./T-1.
    alpha=(0.00504+0.0062*theta)*np.exp(-22.1*theta)
    deltabeta=np.exp(-9.963+0.0372*(T-273.16))
    betaM=B1*np.exp(b/T)/(T*((np.exp(b/T)-1)*(np.exp(b/T)-1)))+B2*f*f
    beta=betaM+deltabeta
    eps2=alpha/f + beta*f
    N=np.sqrt(eps1*eps1+eps2*eps2)

    r=np.sqrt(N)
    THETA=np.arctan2(eps2,eps1)
    nreal=r*(np.cos(0.5*THETA))
    nimm=r*(np.sin(0.5*THETA))    
    return complex(nreal,nimm)
