"""
Implement a few effective medium approximation formulas...


"""

import numpy as np

def maxwell_garnett(eps, mix):
    """Maxwell-Garnett EMA for the refractive index.

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

    If len(m)==2, the first element is taken as the matrix and the second as 
    the inclusion. If len(m)>2, the media are mixed recursively so that the 
    last element is used as the inclusion and the second to last as the 
    matrix, then this mixture is used as the last element on the next 
    iteration, and so on.
    """

    cF = mix[1] / (mix[0]+mix[1]) * (eps[1]-eps[0]) / (eps[1]+2*eps[0])
    return eps[0]*(1.0+2.0*cF) / (1.0-cF)
#    m = np.sqrt(er)
    
#    if len(m) == 2:
#        cF = float(mix[1]) / (mix[0]+mix[1]) * \
#            (m[1]**2-m[0]**2) / (m[1]**2+2*m[0]**2)
#        er = m[0]**2 * (1.0+2.0*cF) / (1.0-cF)
#        m = np.sqrt(er)
#    else:
#        m_last = maxwell_garnett(m[-2:], mix[-2:])
#        mix_last = mix[-2] + mix[-1]
#        m = maxwell_garnett(m[:-2] + (m_last,), mix[:-2] + (mix_last,))
#    return m


def bruggeman(eps, mix):
    """Bruggeman EMA for the refractive index.

    For instructions, see mg_refractive in this module, except this routine
    only works for two components.
    """
    f1 = mix[0]/sum(mix)
    f2 = mix[1]/sum(mix)
    e1 = eps[0] #m[0]**2
    e2 = eps[1] #m[1]**2
    a = -2*(f1+f2)
    b = (2*f1*e1 - f1*e2 + 2*f2*e2 - f2*e1)
    c = (f1+f2)*e1*e2
    #e_eff = (-b - np.sqrt(b**2-4*a*c))/(2*a)
    return (-b - np.sqrt(b**2-4*a*c))/(2*a)
    #return np.sqrt(e_eff)

def sihvola(eps,mix,ni=0.85):
    """Sihvola EMA for the refractive index.

    For instructions, see mg_refractive in this module, except this routine
    only works for two components.

    The original formulation is defaulted to ni=0.85 which has been found to be
    the best for many snow applications. Also, the analitic solution for Sihvola
    modified EMA is way too complicated to be written and computed efficiently:
    a numerically converging solution is applied instead.
    """
    return bruggeman(eps,mix)
    
################################################################################

def n(refractive_indices,volume_fractions,model='bruggeman',ni=0.85):
    if model == 'bruggeman':
        return np.sqrt(bruggeman(refractive_indices**2,volume_fractions))
    elif model == 'sihvola':
        return np.sqrt(sihvola(refractive_indices**2,volume_fractions,ni=ni))
    elif model == 'maxwell_garnett':
        return np.sqrt(maxwell_garnett(refractive_indices**2,volume_fractions))
    else:
        print('Unknown model, fallback to Bruggeman')
        return np.sqrt(bruggeman(refractive_indices**2,volume_fractions))
        
def eps(dielectric_permittivity,volume_fractions,model='bruggeman',ni=0.85):
    if model == 'bruggeman':
        return bruggeman(dielectric_permittivity,volume_fractions)
    elif model == 'sihvola':
        return sihvola(dielectric_permittivity,volume_fractions,ni=ni)
    elif model == 'maxwell_garnett':
        return maxwell_garnett(dielectric_permittivity,volume_fractions)
    else:
        print('Unknown model, fallback to Bruggeman')
        return bruggeman(dielectric_permittivity,volume_fractions)
