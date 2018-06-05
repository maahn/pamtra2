  # -*- coding: utf-8 -
# (c) M. Maahn, 2017

from __future__ import division, absolute_import, print_function

import numpy as np

from pamtra2 import decorators, constants
from . import pyPamtraRadarSimulatorLib



__version__ = '0.1'


def calcSpectralBroadening(
    edr,
    horizontalWind,
    height,
    wavelength,
    beamwidthDeg=0.2,
    integrationTime=60,
    kolmogorov = 0.5,
    verbosity = 0
    ):

  """
  Estimate the spectral broadening due to turbulence and horizontal wind. 

  Parameters
  ----------

  edr : array_like or float
    Eddy dissipation rate in m2/s3
  horizontalWind : array_like or float
    horizontal wind in m/s
  height : array_like or float
    radar range in m
  beamwidthDeg :float
    Full width half power Beamwidth in degrees
  integrationTime : float
    Integration time in seconds
  wavelength : float
    wavelength in m
  kolmogorov : float, optional
    Kolmogorov constant (default 0.5)
  verbosity : integer, optional
    Define verbosity level (default 0)


  Returns
  -------

  specbroad : array_like or float
    Spectral Broadening in m/s (?)
 

  Reference: Shupe, M. D., P. Kollias, M. Poellot, and E. Eloranta, 2008: On Deriving Vertical Air Motions from Cloud Radar Doppler Spectra. J. Atmos. Oceanic Technol., 25, 547–557, doi:10.1175/2007JTECHA1007.1.


  """


  pyPamtraRadarSimulatorLib.report_module.verbose = verbosity

  assert len(edr.shape) == 1, 'edr has to be 1D'
  assert len(horizontalWind.shape) == 1, 'horizontalWind has to be 1D'
  assert len(height.shape) == 1, 'height has to be 1D'

  error,specbroad = pyPamtraRadarSimulatorLib.radar_spectral_broadening(
    edr,
    horizontalWind,
    height,
    beamwidthDeg,
    integrationTime,
    wavelength,
    kolmogorov
    )

  if error>0:
    raise RuntimeError('Error in Fortran routine estimate_spectralbroadening')

  return specbroad


calcSpectralBroadening_ND = decorators.NDto2DtoND(
  referenceIn=0,
  noOfInDimsToKeep=0,
  convertInputs=[0,1,2],
  convertOutputs=[0],
  verbosity=10)(calcSpectralBroadening)

def simulateRadarSpectrum(
  diameterSpec,
  backSpec,
  massSpec,
  rhoSpec,
  areaSpec,
  pathIntegratedAtenuattion,
  height,
  temperature,
  pressure,
  verticalWind,
  spectralBroadening,
  wavelength,
  fallVelocityRelation ='heymsfield10_particles',
  radarMaxV =7.885,
  radarMinV =-7.885,
  radarAliasingNyquistInterv = 1,
  radarNFFT = 256,
  radarAirmotion = True,
  radarAirmotionModel = "constant", #"constant","linear","step"
  radarAirmotionVmin = 0,
  radarAirmotionVmax = 0,
  radarAirmotionLinearSteps = 30,
  radarAirmotionStepVmin = 0.5,
  radarPNoise1000 = -30,
  radarK2 =0.93, # dielectric constant |K|² (always for liquid water by convention) for the radar equation
  radarNAve = 150,
  seed  = 0,
  verbosity = 0
  ):

  """
  Convert a spectrum of hydrometeor backscattering (per hydrometeor)
  as a function of size into a merged spectrum as a function of velocity. 
  Note that the optional shapes of all input variables have to be the same.


  Parameters
  ----------
  diameterSpec : array_like
    Hydrometeo Diameter in m. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  backSpec : array_like
    Hydrometeo backscattering in m2(?). Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  temperature : array_lile or float
    Temperature in K
  pressure: array_lile or float
    Pressure in pa
  verticalWind: array_lile or float
    vertical wind in m/s
  wavelength,
    wavelength in m
  rhoSpec: array_like
    Hydrometeo density in kg/m3. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  fallVelocityRelation: str
  massSpec : array_like
    Hydrometeo mass in kg. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  areaSpec: array_like
    Hydrometeo cross section area in m2. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  radarMaxV : float, optional
    maximum radar nyquist velocity in m/s (default 7.885 )
  radarMinV : float, optional
    mimimum radar nyquist velocity in m/s (default -7.885 )
  radarAliasingNyquistInterv : int, optional, >0
    defines how often teh spectrum is folded to consider aliasing (default 1)
  radarNFFT : int, optional
    bins of the radar spectrum (default 256)
  radarAirmotion : bool, optional 
    consider vertical air otion (default True)
  radarAirmotionModel : str, optional
    "constant","linear","step"
  radarAirmotionVmin : float, optional
     -4.0,
  radarAirmotionVmax: float, optional
     +4.0,
  radarAirmotionLinearSteps: int, optional
     30,
  radarAirmotionStepVmin: float, optional
     0.5,
  radarK2, float, optional
    used dielectric constant |K|²  to estimate Ze, always for liquid water by convention (default 0.93), 


  Returns
  =======
  radar_spectrum : array_like
    merged radar Doppler spectrum for all hydrometeors in m2/mm6, i.e. Ze = sum(radar_spectrum)

  """
  pyPamtraRadarSimulatorLib.report_module.verbose = verbosity

  #In case we are using a fixed seed, make sure we start the beginning again. 
  pyPamtraRadarSimulatorLib.random_module.counter = 0


  radarNFFTAliased = radarNFFT *(1+2*radarAliasingNyquistInterv)

  nHeights = diameterSpec.shape[0]
  nHydro = diameterSpec.shape[1]
  particleSpec = np.zeros((nHeights,nHydro,radarNFFTAliased))

  for hh in range(nHydro):
    assert not np.all(np.isnan(backSpec[:,hh,:]))


    #to do: expose vel_spec in case you need nothing else.

    error,particleSpec[:,hh,:], vel_spec = pyPamtraRadarSimulatorLib.radar_spectrum.get_radar_spectrum(
      diameterSpec[:,hh,:],
      backSpec[:,hh,:],
      temperature,
      pressure,
      verticalWind,
      wavelength,
      rhoSpec[:,hh,:],
      fallVelocityRelation,
      massSpec[:,hh,:],
      areaSpec[:,hh,:],
      radarMaxV,
      radarMinV,
      radarAliasingNyquistInterv,
      radarNFFT,
      radarNFFTAliased,
      radarAirmotion,
      radarAirmotionModel,
      radarAirmotionVmin,
      radarAirmotionVmax,
      radarAirmotionLinearSteps,
      radarAirmotionStepVmin,
      radarK2,
      )
    if error>0:
      raise RuntimeError('Error in Fortran routine radar_spectrum')

  #merge all the hydrometeors
  mergedParticleSpec = np.sum(particleSpec,axis=1)

  #estimate noise from value at 1 km:
  radarPNnoise = 10**(0.1*radarPNoise1000) * (height/1000.)**2


  error,radar_spectrum = pyPamtraRadarSimulatorLib.radar_simulator.simulate_radar(wavelength,
    mergedParticleSpec,
    pathIntegratedAtenuattion,
    spectralBroadening,
    radarPNnoise,
    radarMaxV,
    radarMinV,
    radarNFFT,
    radarNAve,
    radarAliasingNyquistInterv,
    radarK2,
    seed,
    )
  if error>0:
    raise RuntimeError('Error in Fortran routine radar_simulator')

  return radar_spectrum

simulateRadarSpectrum_ND = decorators.NDto2DtoND(
  referenceIn=0,
  noOfInDimsToKeep=2,
  convertInputs=list(range(11)),
  convertOutputs=[0],
  verbosity=1
  )(simulateRadarSpectrum)



'''
NEW ROUTINES
'''

def createRadarSpectrum(
  diameterSpec,
  backSpec,
  massSpec,
  rhoSpec,
  areaSpec,
  height,
  temperature,
  pressure,
  verticalWind,
  wavelength,
  fallVelocityRelation ='heymsfield10_particles',
  radarMaxV =7.885,
  radarMinV =-7.885,
  radarAliasingNyquistInterv = 1,
  radarNFFT = 256,
  verbosity = 0,
  radarAirmotion = True,
  radarAirmotionModel = "constant", #"constant","linear","step"
  radarAirmotionVmin = 0,
  radarAirmotionVmax = 0,
  radarAirmotionLinearSteps = 30,
  radarAirmotionStepVmin = 0.5,
  radarK2 =0.93, # dielectric constant |K|² (always for liquid water by convention) for the radar equation
  ):

    pyPamtraRadarSimulatorLib.report_module.verbose = verbosity


    radarNFFTAliased = radarNFFT *(1+2*radarAliasingNyquistInterv)

    nHeights = height.shape[0]
    particleSpec = np.zeros((radarNFFTAliased,nHeights))


    #to do: expose vel_spec in case you need nothing else.

    error, particleSpec, vel_spec = pyPamtraRadarSimulatorLib.radar_spectrum.get_radar_spectrum(
      diameterSpec,
      backSpec,
      temperature,
      pressure,
      verticalWind,
      wavelength,
      rhoSpec,
      fallVelocityRelation,
      massSpec,
      areaSpec,
      radarMaxV,
      radarMinV,
      radarAliasingNyquistInterv,
      radarNFFT,
      radarNFFTAliased,
      radarAirmotion,
      radarAirmotionModel,
      radarAirmotionVmin,
      radarAirmotionVmax,
      radarAirmotionLinearSteps,
      radarAirmotionStepVmin,
      radarK2,
      )
    if error>0:
        raise RuntimeError('Error in Fortran routine radar_spectrum')


    return particleSpec


def simulateRadarSpectrum(
  height,
  mergedParticleSpec,
  pathIntegratedAtenuattion,
  spectralBroadening,
  wavelength,
  fallVelocityRelation='heymsfield10_particles',
  radarMaxV=7.885,
  radarMinV=-7.885,
  radarAliasingNyquistInterv=1,
  radarNFFT=256,
  radarPNoise1000=-30,
  radarK2=0.93,# dielectric constant |K|² (always for liquid water by convention) for the radar equation
  radarNAve=150,
  seed=0,
  verbosity=0
  ):

  """
  Convert a spectrum of hydrometeor backscattering (per hydrometeor)
  as a function of size into a merged spectrum as a function of velocity. 
  Note that the optional shapes of all input variables have to be the same.


  Parameters
  ----------
  diameterSpec : array_like
    Hydrometeo Diameter in m. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  backSpec : array_like
    Hydrometeo backscattering in m2(?). Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  temperature : array_lile or float
    Temperature in K
  pressure: array_lile or float
    Pressure in pa
  verticalWind: array_lile or float
    vertical wind in m/s
  wavelength,
    wavelength in m
  rhoSpec: array_like
    Hydrometeo density in kg/m3. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  fallVelocityRelation: str
  massSpec : array_like
    Hydrometeo mass in kg. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  areaSpec: array_like
    Hydrometeo cross section area in m2. Shape ([optional dimensions], hydrometeor, hydrometeorBin)
  radarMaxV : float, optional
    maximum radar nyquist velocity in m/s (default 7.885 )
  radarMinV : float, optional
    mimimum radar nyquist velocity in m/s (default -7.885 )
  radarAliasingNyquistInterv : int, optional, >0
    defines how often teh spectrum is folded to consider aliasing (default 1)
  radarNFFT : int, optional
    bins of the radar spectrum (default 256)
  radarAirmotion : bool, optional 
    consider vertical air otion (default True)
  radarAirmotionModel : str, optional
    "constant","linear","step"
  radarAirmotionVmin : float, optional
     -4.0,
  radarAirmotionVmax: float, optional
     +4.0,
  radarAirmotionLinearSteps: int, optional
     30,
  radarAirmotionStepVmin: float, optional
     0.5,
  radarK2, float, optional
    used dielectric constant |K|²  to estimate Ze, always for liquid water by convention (default 0.93), 


  Returns
  =======
  radar_spectrum : array_like
    merged radar Doppler spectrum for all hydrometeors in m2/mm6, i.e. Ze = sum(radar_spectrum)

  """
  pyPamtraRadarSimulatorLib.report_module.verbose = verbosity

  #In case we are using a fixed seed, make sure we start the beginning again. 
  pyPamtraRadarSimulatorLib.random_module.counter = 0

  #estimate noise from value at 1 km:
  radarPNnoise = 10**(0.1*radarPNoise1000) * (height/1000.)**2


  error, radar_spectrum = pyPamtraRadarSimulatorLib.radar_simulator.simulate_radar(wavelength,
    mergedParticleSpec,
    pathIntegratedAtenuattion,
    spectralBroadening,
    radarPNnoise,
    radarMaxV,
    radarMinV,
    radarNFFT,
    radarNAve,
    radarAliasingNyquistInterv,
    radarK2,
    seed,
    )
  if error>0:
    raise RuntimeError('Error in Fortran routine radar_simulator')

  return radar_spectrum



