"""
Implement a couple of mixing formulas...


"""

def maxwell_garnett(m, mix):
    """Maxwell-Garnett EMA for the refractive index.

    Args:
       m: Tuple of the complex refractive indices of the media.
       mix: Tuple of the volume fractions of the media, len(mix)==len(m)
            (if sum(mix)!=1, these are taken relative to sum(mix))

    Returns:
       The Maxwell-Garnett approximation for the complex refractive index of 
       the effective medium

    If len(m)==2, the first element is taken as the matrix and the second as 
    the inclusion. If len(m)>2, the media are mixed recursively so that the 
    last element is used as the inclusion and the second to last as the 
    matrix, then this mixture is used as the last element on the next 
    iteration, and so on.
    """

    if len(m) == 2:
        cF = float(mix[1]) / (mix[0]+mix[1]) * \
            (m[1]**2-m[0]**2) / (m[1]**2+2*m[0]**2)
        er = m[0]**2 * (1.0+2.0*cF) / (1.0-cF)
        m = np.sqrt(er)
    else:
        m_last = maxwell_garnett(m[-2:], mix[-2:])
        mix_last = mix[-2] + mix[-1]
        m = maxwell_garnett(m[:-2] + (m_last,), mix[:-2] + (mix_last,))
    return m


def bruggeman(m, mix):
    """Bruggeman EMA for the refractive index.

    For instructions, see mg_refractive in this module, except this routine
    only works for two components.
    """
    f1 = mix[0]/sum(mix)
    f2 = mix[1]/sum(mix)
    e1 = m[0]**2
    e2 = m[1]**2
    a = -2*(f1+f2)
    b = (2*f1*e1 - f1*e2 + 2*f2*e2 - f2*e1)
    c = (f1+f2)*e1*e2
    e_eff = (-b - np.sqrt(b**2-4*a*c))/(2*a)
    return np.sqrt(e_eff)

def sihvola()::
    """Sihvola EMA for the refractive index.

    For instructions, see mg_refractive in this module, except this routine
    only works for two components.

    The original formulation is defaulted to ni=0.85 which has been found to be
    the best for many snow applications
    """
    return bruggeman(m,mix)
