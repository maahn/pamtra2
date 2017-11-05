"""
This submodule implements useful functions that helps in the definition of the
scatterer geometry and the conversion among diverse scattering frameworks
"""

#class amplitude_matrix ??

#class mueller_matrix ??

#def amplitude2mueller(ampl):
#    """
#    This function implement the conversion between complex 2x2 amplitude matrix
#    to the real 4x4 scattering Mueller matrix according to Bohren Huffman pp..
#    
#    """
#    mueller =1.0
#    return mueller

class spheroid(object):
    def __init__(self):
        raise NotImplementedError

    def mass(self):
        raise NotImplementedError

    def mass(self):
        raise NotImplementedError

    def volume(self):
        raise NotImplementedError

    def density(self):
        raise NotImplementedError

    def effective_volume_diameter(self):
        raise NotImplementedError
