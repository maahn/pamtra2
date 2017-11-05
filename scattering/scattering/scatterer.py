""" scattering.utilities.py

This submodule implements useful functions that helps in the definition of the
scatterer geometry and the conversion among diverse scattering frameworks

"""

import numpy as np

from ...refractive import refractive.utilities as ref_utils

#class amplitude_matrix

#class mueller_matrix

#def amplitude2mueller(ampl):
#    """
#    This function implement the conversion between complex 2x2 amplitude matrix
#    to the real 4x4 scattering Mueller matrix according to Bohren Huffman pp..
#    
#    """
#    mueller =1.0
#    return mueller

class scatterer(object):
    def __init__(self,
                 diameter = 1.,
                 frequency = 1.,
                 refractive_index=None,
                 dielectric_permittivity=None
                 ):
        self.diameter = diameter
        self.frequency = frequency
        self.set_dielectric_properties(refractive_index,dielectric_permittivity)
        
    def set_dielectric_properties(refractive_index,dielectric_permittivity)
        if (refractive_index is None):
            if (dielectric_permittivity is None):
                self.refractive_index = np.array(complex(1.0,0.0))
                self.dielectric_permittivity = ref_utils.n2eps(self.refractive_index)
            else:
                self.dielectric_permittivity = dielectric_permittivity
                self.refractive_index = ref_utils.eps2n(self.dielectric_permittivity)
        elif (dielectric_permittivity is None):
            self.refractive_index = refractive_index
            self.dielectric_permittivity = ref_utils.n2eps(self.refractive_index)
        else:
            raise AttributeError('Both dielectric permittivity and refractive index has been defined')
            
            
        
    

