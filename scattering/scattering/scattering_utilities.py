"""
This submodule implements useful functions that helps in the definition of the
scatterer geometry and the conversion among diverse scattering frameworks
"""

import numpy as np

class amplitude_matrix(object):
    """ Complex 2x2 amplitude matrix as defined in Bohren and Huffman
    | S2   S3 |
    |         |
    | S4   S1 |
    
    
    """
    def __init__(self,S):
        self.matrix = S+0.0j # ensure it is complex valued

    @property
    def S1(self):
        return self.matrix[1,1]
        
    @property
    def S2(self):
        return self.matrix[0,0]
    
    @property
    def S3(self):
        return self.matrix[0,1]
    
    @property
    def S4(self):
        return self.matrix[1,0]
        
    def to_mueller(self):
        return amplitude2mueller(self)

class scattering_matrix(object):
    """ Real 4x4 scattering matrix as defined in Bohren and Huffman
    
    """

    def __init__(self,Z):
        self.matrix = Z.real # ensure it is real (TODO maybe it is better to just have a check)


def amplitude2mueller(ampl):
    """ This function implement the conversion between complex 2x2 amplitude matrix
    to the real 4x4 scattering Mueller matrix according to Bohren Huffman pp..
    
    """
    
    mueller = np.ndarray((4,4))
    S1_2 = (ampl.S1*ampl.S1.conj).real
    S2_2 = (ampl.S2*ampl.S2.conj).real
    S3_2 = (ampl.S3*ampl.S3.conj).real
    S4_2 = (ampl.S4*ampl.S4.conj).real
    mueller[0,0] = 0.5*(S2_2+S1_2+S4_2+S3_2)
    mueller[0,1] = 0.5*(S2_2-S1_2+S4_2+S3_2)
    mueller[0,2] = (ampl.S2*ampl.S3.conj+ampl.S1*ampl.S4.conj).real
    mueller[0,3] = (ampl.S2*ampl.S3.conj-ampl.S1*ampl.S4.conj).imag

    mueller[1,0] = 0.5*(S2_2-S2_2-S4_2+S3_2)
    mueller[1,1] = 0.5*(S2_2+S2_2-S4_2-S3_2)
    mueller[1,2] = (ampl.S2*ampl.S3.conj-ampl.S1*ampl.S4.conj).real
    mueller[1,3] = (ampl.S2*ampl.S3.conj+ampl.S1*ampl.S4.conj).imag

    mueller[2,0] = (ampl.S2*ampl.S4.conj+ampl.S1*ampl.S3.conj).real
    mueller[2,1] = (ampl.S2*ampl.S4.conj-ampl.S1*ampl.S3.conj).real
    mueller[2,2] = (ampl.S1*ampl.S2.conj+ampl.S3*ampl.S4.conj).real
    mueller[2,3] = (ampl.S2*ampl.S1.conj+ampl.S4*ampl.S3.conj).imag

    mueller[3,0] = (ampl.S2.conj*ampl.S4+ampl.S3.conj*ampl.S1).imag
    mueller[3,1] = (ampl.S2.conj*ampl.S4-ampl.S3.conj*ampl.S1).imag
    mueller[3,2] = (ampl.S1*ampl.S2.conj-ampl.S3*ampl.S4.conj).imag
    mueller[3,3] = (ampl.S1*ampl.S2.conj-ampl.S3*ampl.S4.conj).real

    return scattering_matrix(mueller)

size_parameter = lambda radius, wavelength: 2.0*np.pi*radius/wavelength

class spheroid(object):
    def __init__(self):
        raise NotImplementedError

    def mass(self):
        raise NotImplementedError

    def volume(self):
        raise NotImplementedError

    def density(self):
        raise NotImplementedError

    def effective_volume_diameter(self):
        raise NotImplementedError

def scattering_angle(theta_inc,theta_sca,phi_inc,phi_sca):
    """ Calculates the scattering angle in radians given the full set of four
    angles that defines the scattering geometry from the incident and the
    scattered wave directions
    
    The reference frame is assumed to be set with the polar angle theta
    measured from the vertical axis z, and the azimuth angle phi measured from
    the axis x???
    
    """
    acos_th = np.sin(theta_inc)*np.sin(theta_sca)*np.cos(phi_inc-phi_sca) + np.cos(theta_inc)*np.cos(theta_sca)
    return np.arccos(acos_th)