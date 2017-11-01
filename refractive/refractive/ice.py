"""
Matzler (2006)

Temperature in Kelvin
Frequency in GHz

"""

import numpy as np

B1 = 0.0207
b = 335.
B2 = 1.16e-11
c = 299792458.

def warren_brandt_2008(T,f):
    print("Not implemented yet falling back to Matzler 2006")
    return matzler_2006(T,f)

def matzler_2006(T,f):
    eps1  = 3.1884+(T-273)*9.1e-4
    theta = 300./T-1.
    alpha=(0.00504+0.0062*theta)*np.exp(-22.1*theta)
    deltabeta=np.exp(-9.963+0.0372*(T-273.16))
    betaM=B1*np.exp(b/T)/(T*((np.exp(b/T)-1)*(np.exp(b/T)-1)))+B2*f*f
    beta=betaM+deltabeta
    eps2=alpha/f + beta*f
    return eps1 + 1j*eps2

#######################################################################################################

def eps(T,f,what="Matzler_2006"):
    if (what == "Matzler_2006"):
        return matzler_2006(T,f)
    else:
        print("I do not recognize the ice refractive index specification, falling back to Matzler 2006")
        return matzler_2006(T,f)

def n(T,f,what="Matzler_2006"):
    return np.sqrt(eps(T,f,what))

#######################################################################################################

if __name__ == "__main__":
    import sys
    n(float(sys.argv[1]),float(sys.argv[2]))
