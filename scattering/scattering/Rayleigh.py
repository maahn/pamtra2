import numpy as np

def rayleigh(diameter,K2,frequency):

  """
  This is Rayleigh for spheres
  """
  C = 299792458.

  K2 = np.asarray(K2)
  diameter = np.asarray(diameter)

  wavelength = C / (frequency*1e9)  
  prefactor = np.pi**5 * K2 / wavelength**4
  back_spec =  prefactor[:,np.newaxis] * diameter**6

return back_spec
