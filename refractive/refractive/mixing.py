"""
Implement a few effective medium approximation formulas...

Both material dielectric models and EMAs usually work with dielec permittivity
and refractive index is derived as sqrt(eps), hence this modules is implemented
to always take eps as arguments of both the eps and n functions for
computational efficiency

"""

import numpy as np

def maxwell_garnett(eps, mix):
    """Maxwell-Garnett EMA for the refractive index.
       AKA: Rayleigh mixing formula [Sihvola, 1989]
       
       The MG relation for effective medium approximation is derived for the 
       case of small, highly diluted, spatially well separated, optically soft
       and spherical inclusions in a medium. Also the size of the considered
       volume should be small.
       
       The MG relation is asymmetric, the order of medium and inclusion in the
       formula counts. This is a consequence of the aformentioned assumptions.
       
       'J.C.M. Garnett, Philos. Trans. R. Soc. 203, 385 (1904); 205, 237 (1906).'

    Parameters
    ----------
        eps: Tuple of complex
            dielectric permittivity of the media
        mix: Tuple of float
            the volume fractions of the media, len(mix)==len(m)
            (if sum(mix)!=1, these are taken relative to sum(mix))

    Returns:
    ----------
        The Maxwell-Garnett approximation for the complex refractive index of 
        the effective medium

    The first element of the eps and mix tuples is taken as the matrix and the
    second as the inclusion.
    """

    cF = mix[1] / (mix[0]+mix[1]) * (eps[1]-eps[0]) / (eps[1]+2*eps[0])
    return eps[0]*(1.0+2.0*cF) / (1.0-cF)

def bruggeman(eps, mix):
    """Bruggeman EMA for the refractive index.

    For instructions, see mg_refractive in this module, except this routine
    only works for two components.
    
    'D.A.G. Bruggeman, Ann. Phys. 24, 636 (1935)'
    
    Bruggeman model has the advantage with respect to MG of beeing symmetric
    """
    f1 = mix[0]/sum(mix)
    f2 = mix[1]/sum(mix)
    e1 = eps[0] #m[0]**2
    e2 = eps[1] #m[1]**2
    a = -2*(f1+f2)
    b = (2*f1*e1 - f1*e2 + 2*f2*e2 - f2*e1)
    c = (f1+f2)*e1*e2
    return (-b - np.sqrt(b**2-4*a*c))/(2*a)

def sihvola(eps,mix,ni=0.85):
    """Sihvola EMA for the refractive index.

    For instructions, see mg_refractive in this module, except this routine
    only works for two components.

    The original formulation is defaulted to ni=0.85 which has been found to be
    the best for many snow applications. Also, the analitic solution for Sihvola
    modified EMA is way too complicated to be written and computed efficiently:
    a numerically converging solution is applied instead.
    
    A.H. Sihvola 'Self-Consistency Aspects of Dielectric Mixing Theories'
    IEEE Trans. Geos. Rem. Sens. vol 27, n 4, 1989
    """
    print('Sihvola EMA is not yet implemented fallback to Bruggeman')
    return bruggeman(eps,mix)
    
################################################################################

def n(refractive_indices,volume_fractions,model='bruggeman',ni=0.85):
    if model == 'Bruggeman':
        return np.sqrt(bruggeman(refractive_indices**2,volume_fractions))
    elif model == 'Sihvola':
        return np.sqrt(sihvola(refractive_indices**2,volume_fractions,ni=ni))
    elif model == 'Maxwell_Garnett':
        return np.sqrt(maxwell_garnett(refractive_indices**2,volume_fractions))
    else:
        print('Unknown model, fallback to Bruggeman')
        return np.sqrt(bruggeman(refractive_indices**2,volume_fractions))
        
def eps(dielectric_permittivity,volume_fractions,model='bruggeman',ni=0.85):
    if model == 'Bruggeman':
        return bruggeman(dielectric_permittivity,volume_fractions)
    elif model == 'Sihvola':
        return sihvola(dielectric_permittivity,volume_fractions,ni=ni)
    elif model == 'Maxwell_garnett':
        return maxwell_garnett(dielectric_permittivity,volume_fractions)
    else:
        print('Unknown model, fallback to Bruggeman')
        return bruggeman(dielectric_permittivity,volume_fractions)
