# -*- coding: utf-8 -
# (c) M. Maahn, 2017

from . import pyPamtraRadarMomentsLib
from pamtra2 import decorators

import numpy as np


__version__ = '0.1'


def calc_hildebrandSekhon(spectrum, noAve = 1,verbose=0):

  """
  Calculate the mean and maximum of noise of the linear radar spectrum following Hildebrand and Sekhon 1974.

  Parameters
  ----------

  spectrum : array_like
      linear radar spectrum. Can be 1D or 2D.
  noAve : int, optional
      number of averages (default 1)
  verbose : int, optional
      verbosity level (default 0)

  Returns
  -------

  meanNoise : array_like or float
    spectral mean noise level in linear units, noise power = meanNoise*nFFT
  maxNoise : array_like or float
    spectral maximum noise level in linear units
 
  """

  spectrum = np.asarray(spectrum)
  specShape = np.shape(spectrum)

  assert len(specShape) <=2,  'spectrum must not have more than two dimensions (height, nfft)'

  if len(specShape) == 1:
    spectrum = spectrum.reshape((1,specShape[0]))

  pyPamtraRadarMomentsLib.report_module.verbose = verbose

  error, meanNoise, maxNoise = pyPamtraRadarMomentsLib.hildebrand_sekhon(spectrum,noAve)
  
  if error>0:
    raise RuntimeError('Error in Fortran routine hildebrand_sekhon')

  return meanNoise, maxNoise



@decorators.NDto2DtoND(referenceIn=0,convertInputs=[0],convertOutputs=[0,1,2,3,4])
def calc_radarMoments(spectrum,
    verbose = 0, 
    maxV = 7.885, 
    minV = -7.885, 
    noAve = 150, 
    nPeaks = 3, 
    noiseDistanceFactor = 0, 
    specNoiseMean = None, #linear noise per spectral bin in mm6/m3
    specNoiseMax = None, #linear noise per spectral bin in mm6/m3
    peakMinSnr=-10, 
    peakMinBins=2, 
    smoothSpectrum= True,
    useWiderPeak=False,
    receiverMiscalibration=0,
    ):
  
  """
  Calculates the moments, slopes and edges of the linear radar spectrum. 

  Parameters
  ----------

  spectrum : array_like
      linear radar spectrum [mm⁶/m³]. Can be 1D or 2D.
  verbose : int, optional
      verbosity level (default 0)
  maxV : float, optional
    MVimum Nyquist Velocity (default  7.885)
  minV : float, optional
    Minimum Nyquist Velocity (default  -7.885)
  noAve : int, optional
    No of averages per spectrum (default  150)
  npeaks : int, optional
    No of peaks which should be determined (default  3)
  noiseDistanceFactor : float, optional
    factor between noise and noise max. If 0, noiase max is obtained using calc_hildebrandSekhon (default  0)
  specNoiseMean : array_like or float, optional
    linear mean noise per spectral bin in mm: mm6/m3. If None, it is determined using calc_hildebrandSekhon (default  None)
  specNoiseMax : array_like or float, optional
    linear maximum noise per spectral bin in mm: mm6/m3. If None, it is determined using calc_hildebrandSekhon (default  None)
  peakMinSnr: float, optional
    minimum linear SNR for each peak (default 1.2)
  peakMinBins: int, optional
    minimal number of bins per peak (default 2)
  smooth_spectrum: bool, optional
    smooth spectrum before estiamting moments (default  True)
  useWiderPeak: bool, optional
    include edges into peak (default False)
  receiverMiscalibration, float, optional
    simulate a wrong radar receiver calibration [dB] (default 0)
  Returns
  -------

  spectrumOut : array_like
    radar spectrum with noise removed [mm⁶/m³]
  moments : array_like
    0th - 4th moment [mm⁶/m³, m/s, m/s,-,-]
  slope : array_like
    slope of the peak [dB/(m/s)]
  edge : array_like
    left(0) and right(1) edge the peak [m/s]
  quality : array_like
    quality flag: 1st byte: aliasing; 2nd byte: more peaks present; 7th: no peak found; 8th: principal peak isolated
  noiseMean : float
    mean intgrated noise power in linear units

 """

  spectrum = np.asarray(spectrum)
  specShape = np.shape(spectrum)


  assert len(specShape) <=2, 'spectrum must not have more than two dimensions (height, nfft)'

  if (specNoiseMean is None):
    specNoiseMean, specNoiseMaxHilde = calc_hildebrandSekhon(spectrum, noAve = noAve,verbose=verbose)
  if (specNoiseMax is None):
    if noiseDistanceFactor > 0:
      specNoiseMax = specNoiseMean * noiseDistanceFactor
    else:
      specNoiseMax = specNoiseMaxHilde


  specNoiseMean = np.asarray(specNoiseMean)
  specNoiseMax = np.asarray(specNoiseMax)

  pyPamtraRadarMomentsLib.report_module.verbose = verbose

  #apply a receiver miscalibration:
  if receiverMiscalibration != 0:
    spectrum = spectrum * 10**(0.1*receiverMiscalibration)
    specNoiseMax = specNoiseMax * 10**(0.1*receiverMiscalibration)
    specNoiseMean = specNoiseMean * 10**(0.1*receiverMiscalibration)

  if len(specShape) == 1:
    spectrum = spectrum.reshape((1,specShape[0]))
    specNoiseMean = specNoiseMean.reshape(1)
    specNoiseMax = specNoiseMax.reshape(1)


  noiseMean  = specNoiseMean*spectrum.shape[-1] #nFFT

  assert (len(spectrum.shape) -1) == len(specNoiseMean.shape), 'shape of spectrum and specNoiseMean does not match'
  assert (len(spectrum.shape) -1) == len(specNoiseMax.shape), 'shape of spectrum and specNoiseMax does not match'


  output = pyPamtraRadarMomentsLib.calc_moments.calc_moments_column(
    nPeaks,
    spectrum,
    specNoiseMean,
    specNoiseMax,
    maxV,
    minV,
    smoothSpectrum,
    useWiderPeak,
    peakMinBins,
    peakMinSnr,
    )

  error,spectrumOut,moments,slope,edge,quality = output
  if error>0:
    raise RuntimeError('Error in Fortran routine calc_moments_column')

  return spectrumOut,moments,slope,edge,quality,noiseMean

