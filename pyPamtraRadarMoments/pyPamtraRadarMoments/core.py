# -*- coding: utf-8 -
# (c) M. Maahn, 2017

from . import pyPamtraRadarMomentsLib
# from pamtra2 import decorators

import numpy as np


__version__ = '0.1'


def calc_hildebrandSekhon(spectrum, radarNAve=1, verbosity=0):
    """
    Calculate the mean and maximum of noise of the linear radar spectrum
    following Hildebrand and Sekhon 1974.

    Parameters
    ----------

    spectrum : array_like [mm⁶/m³/(m/s)]
        linear radar spectrum. Can be 1D or 2D.
    radarNAve : int, optional
        number of averages (default 1)
    verbosity : int, optional
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

    assert len(
        specShape) <= 2,  ('spectrum must not have more than two dimensions' +
                           ' (height, nfft)')

    if len(specShape) == 1:
        spectrum = spectrum.reshape((1, specShape[0]))

    pyPamtraRadarMomentsLib.report_module.verbose = verbosity

    error, meanNoise, maxNoise = pyPamtraRadarMomentsLib.hildebrand_sekhon(
        spectrum, radarNAve)

    if error > 0:
        raise RuntimeError('Error in Fortran routine hildebrand_sekhon')

    return meanNoise, maxNoise


# @decorators.NDto2DtoND(referenceIn=0,convertInputs=[0],convertOutputs=[0,1,2,3,4])
def calc_radarMoments(spectrum,
                      verbosity=0,
                      radarMaxV=7.885,
                      radarMinV=-7.885,
                      radarNAve=150,
                      momentsNPeaks=3,
                      momentsNoiseDistanceFactor=0,
                      momentsSpecNoiseMean=None,
                      momentsSpecNoiseMax=None,
                      momentsPeakMinSnr=-10,
                      momentsPeakMinBins=2,
                      momentsSmoothSpectrum=True,
                      momentsUseWiderPeak=False,
                      momentsReceiverMiscalibration=0,
                      ):
    """
    Calculates the moments, slopes and edges of the linear radar spectrum.

    Parameters
    ----------

    spectrum : array_like
        linear radar spectrum [mm⁶/m³/(m/s)]. Can be 1D or 2D.
    verbosity : int, optional
        verbosity level (default 0)
    radarMaxV : float, optional
      MVimum Nyquist Velocity (default  7.885)
    radarMinV : float, optional
      Minimum Nyquist Velocity (default  -7.885)
    radarNAve : int, optional
      No of averages per spectrum (default  150)
    momentsNpeaks : int, optional
      No of peaks which should be determined (default  3)
    momentsNoiseDistanceFactor : float, optional
      factor between noise and noise max. If 0, noiase max is obtained using
      calc_hildebrandSekhon (default  0)
    specNoiseMean : array_like or float, optional
      linear mean noise per spectral bin in mm: mm6/m3. If None, it is
      determined using calc_hildebrandSekhon (default  None)
    momentsSpecNoiseMean : array_like or float, optional
      linear maximum noise per spectral bin in mm: mm6/m3. If None, it is
      determined using calc_hildebrandSekhon (default  None)
    momentsPeakMinSnr: float, optional
      minimum linear SNR for each peak (default 1.2)
    momentsPeakMinBins: int, optional
      minimal number of bins per peak (default 2)
    momentsSmoothSpectrum: bool, optional
      smooth spectrum before estiamting moments (default  True)
    momentsUseWiderPeak: bool, optional
      include edges into peak (default False)
    momentsReceiverMiscalibration, float, optional
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
      quality flag: 1st byte: aliasing; 2nd byte: more peaks present; 7th: no
      peak found; 8th: principal peak isolated
    noiseMean : float
      mean intgrated noise power in linear units

   """

    spectrum = np.asarray(spectrum)
    specShape = np.shape(spectrum)

    assert len(
        specShape) <= 2, ('spectrum must not have more than two dimensions' +
                          ' (height, nfft)')

    if (momentsSpecNoiseMean is None):
        momentsSpecNoiseMean, momentsSpecNoiseMaxHilde = calc_hildebrandSekhon(
            spectrum, radarNAve=radarNAve, verbosity=verbosity)
    if (momentsSpecNoiseMax is None):
        if momentsNoiseDistanceFactor > 0:
            momentsSpecNoiseMax = momentsSpecNoiseMean * momentsNoiseDistanceFactor
        else:
            momentsSpecNoiseMax = momentsSpecNoiseMaxHilde

    momentsSpecNoiseMean = np.asarray(momentsSpecNoiseMean)
    momentsSpecNoiseMax = np.asarray(momentsSpecNoiseMax)

    import pdb

    pyPamtraRadarMomentsLib.report_module.verbose = verbosity

    # apply a receiver miscalibration:
    if momentsReceiverMiscalibration != 0:
        spectrum = spectrum * 10**(0.1*momentsReceiverMiscalibration)
        momentsSpecNoiseMax = momentsSpecNoiseMax * \
            10**(0.1*momentsReceiverMiscalibration)
        momentsSpecNoiseMean = momentsSpecNoiseMean * \
            10**(0.1*momentsReceiverMiscalibration)

    if len(specShape) == 1:
        spectrum = spectrum.reshape((1, specShape[0]))
        momentsSpecNoiseMean = momentsSpecNoiseMean.reshape(1)
        momentsSpecNoiseMax = momentsSpecNoiseMax.reshape(1)

    noiseMean = momentsSpecNoiseMean*spectrum.shape[-1]  # nFFT

    assert (len(spectrum.shape) - 1) == len(
        momentsSpecNoiseMean.shape
    ), 'shape of spectrum and momentsSpecNoiseMean does not match'
    assert (len(spectrum.shape) - 1) == len(
        momentsSpecNoiseMax.shape
    ), 'shape of spectrum and momentsSpecNoiseMax does not match'

    output = pyPamtraRadarMomentsLib.calc_moments.calc_moments_column(
        momentsNPeaks,
        spectrum,
        momentsSpecNoiseMean,
        momentsSpecNoiseMax,
        radarMaxV,
        radarMinV,
        momentsSmoothSpectrum,
        momentsUseWiderPeak,
        momentsPeakMinBins,
        momentsPeakMinSnr,
    )

    error, spectrumOut, moments, slope, edge, quality = output
    if error > 0:
        raise RuntimeError('Error in Fortran routine calc_moments_column')

    return spectrumOut, moments, slope, edge, quality, noiseMean
