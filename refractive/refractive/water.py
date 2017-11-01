import numpy as np
"""
Ellison et al. (2005)
IF I REMEMBER CORRECTLY

"""


a0 = 5.7230
a1 = 2.2379e-2
a2 = -7.1237e-4
a3 = 5.0478
a4 = -7.0315e-2
a5 = 6.0059e-4
a6 = 3.6143
a7 = 2.8841e-2
a8 = 1.3652e-1
a9 = 1.4825e-3
a10 = 2.4166e-4

def ellison(Tk,f):
    T = Tk-273.15
    es=(37088.6-82.168*T)/(421.854+T)
    einf=a6+a7*T
    e1=a0+a1*T+a2*T*T
    ni1=(45+T)/(a3+a4*T+a5*T*T)
    ni2=(45+T)/(a8+a9*T+a10*T*T)
    A1=f/ni1
    A2=f/ni2
    eps1=(es-e1)/(1+A1*A1)+(e1-einf)/(1+A2*A2)+einf
    eps2=(es*A1-e1*A1)/(1+A1*A1)+(e1*A2-einf*A2)/(1+A2*A2)
    return eps1 + 1j*eps2

#######################################################################################################

def eps(T,f,what="ellison"):
    if (what == "ellison"):
        return ellison(T,f)
    else:
        print("I do not recognize the ice refractive index specification, falling back to ellison")
        return ellison(T,f)

def n(T,f,what="ellison"):
    return np.sqrt(eps(T,f,what))

#######################################################################################################

if __name__ == "__main__":
    import sys
    n(float(sys.argv[1]),float(sys.argv[2]))
