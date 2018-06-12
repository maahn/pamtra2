# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr

import pyPamtraRadarSimulator
import pyPamtraRadarMoments

from .. import units
from .. import helpers
from .core import instrument


class _radar(instrument):

    def _calcPIA(self):
        print('PIA = 0!')
        return xr.zeros_like(self.profile.height)


class simpleRadar(_radar):
    def __init__(
        self,
        parent,
        frequencies='all',
        radarK2=0.93,
    ):
        super().__init__(
            parent,
            frequencies=frequencies,
            radarK2=radarK2,
        )

    def solve(self):

        crossSections = self.parent.getIntegratedScatteringCrossSections(
            crossSections=['backscatterCrossSection'],
            frequencies=self.frequencies,
        )

        back = crossSections['backscatterCrossSection']
        wavelength = self.parent.profile.wavelength.sel(
            frequency=self.frequencies)
        K2 = self.settings['radarK2']
        PIA = self.parent.profile.pathIntegratedAtenuattion.sel(
            frequency=self.frequencies)
        Ze_back = 1.e18*(1.e0/(K2*np.pi**5))*back*(wavelength)**4
        Ze_back = Ze_back/10**(0.1*PIA)

        self.results['radarReflectivty'] = 10*np.log10(Ze_back)
        self.results['radarReflectivty'].attrs['unit'] = units.units[
            'radarReflectivty']


        return self.results


class dopplerRadarPamtra(_radar):

    def __init__(
        self,
        parent,
        frequencies='all',
        radarMaxV=7.885,
        radarMinV=-7.885,
        radarAliasingNyquistInterv=4,
        radarNFFT=256,
        verbosity=0,
        radarAirmotion=True,
        radarAirmotionModel="constant",  # "constant","linear","step"
        radarAirmotionVmin=0,
        radarAirmotionVmax=0,
        radarAirmotionLinearSteps=30,
        radarAirmotionStepVmin=0.5,
        # dielectric constant |K|² (always for liquid water by convention) for the radar equation
        radarK2=0.93,
        radarBeamwidthDeg=0.2,
        radarIntegrationTime=60,
        radarPNoise1000=-30,
        radarNAve=150,
        momentsNPeaks=2,
        momentsNoiseDistanceFactor=0,
        momentsSpecNoiseMean=None,
        momentsSpecNoiseMax=None,
        momentsPeakMinSnr=-10,
        momentsPeakMinBins=2,
        momentsSmoothSpectrum=True,
        momentsUseWiderPeak=False,
        momentsReceiverMiscalibration=0,
        seed=0,
    ):

        super().__init__(
            parent,
            frequencies=frequencies,
            radarMaxV=radarMaxV,
            radarMinV=radarMinV,
            radarAliasingNyquistInterv=radarAliasingNyquistInterv,
            radarNFFT=radarNFFT,
            verbosity=verbosity,
            radarAirmotion=radarAirmotion,
            radarAirmotionModel=radarAirmotionModel,
            radarAirmotionVmin=radarAirmotionVmin,
            radarAirmotionVmax=radarAirmotionVmax,
            radarAirmotionLinearSteps=radarAirmotionLinearSteps,
            radarAirmotionStepVmin=radarAirmotionStepVmin,
            radarK2=radarK2,
            radarPNoise1000=radarPNoise1000,
            radarNAve=radarNAve,
            radarBeamwidthDeg=radarBeamwidthDeg,
            radarIntegrationTime=radarIntegrationTime,
            seed=seed,
            momentsNPeaks=momentsNPeaks,
            momentsNoiseDistanceFactor=momentsNoiseDistanceFactor,
            momentsSpecNoiseMean=momentsSpecNoiseMean,
            momentsSpecNoiseMax=momentsSpecNoiseMax,
            momentsPeakMinSnr=momentsPeakMinSnr,
            momentsPeakMinBins=momentsPeakMinBins,
            momentsSmoothSpectrum=momentsSmoothSpectrum,
            momentsUseWiderPeak=momentsUseWiderPeak,
            momentsReceiverMiscalibration=momentsReceiverMiscalibration,
        )

        if len(frequencies) > 1:
            raise NotImplementedError('Sorry, the radar simulator can handle '
                                      'only one frequency at the moment (due '
                                      'to missing support for vectorized '
                                      'settings.)')

    def solve(self):

        for name in self.hydrometeorProfiles.keys():
            if len(self.hydrometeorProfiles[name].sizeBin) <= 1:
                raise ValueError('nbins must be greater than 1 for the'
                                 'spectral radar simulator!')

        self._calcRadarSpectrum()
        self._simulateRadar()
        self._calcMoments()

        for vv in self.results.variables:
            self.results[vv].attrs['unit'] = units.units[vv]

        return self.results

    def _calcRadarSpectrum(self):

        hydroVars = [
            'sizeCenter',
            'sizeBoundsWidth',
            'backscatterCrossSection',
            'fallVelocity',
            'sizeDistribution',
        ]
        profileVars = [
            'verticalWind',
            'wavelength',
        ]
        funcArgVars = [
            'sizeCenter',
            'sizeBoundsWidth',
            'bcsINTEGRATED',
            'fallVelocity',
        ] + profileVars

        radarSpecs = []
        for name in self.hydrometeorProfiles.keys():

            # would be great if this could become an attribute of
            # parent.hydrometeors[name].profile and we could use
            # self.hydrometeorProfiles.
            mergedProfile = self.parent.hydrometeors[
                name
            ].getProfileWithParentAllBroadcasted(
                variables=hydroVars,
                parentVariables=profileVars,
                exclude=['sizeBin'],
            ).sel(frequency=self.frequencies)

            mergedProfile = mergedProfile.stack(
                merged=helpers.concatDicts(
                    self.parent.coords['additional'],
                    self.parent.coords['layer'],
                    self.parent.coords['frequency']
                )
            )

            mergedProfile['bcsINTEGRATED'] = (
                mergedProfile['backscatterCrossSection'] *
                mergedProfile['sizeDistribution']
            )

            argNames, kwargNames = helpers.provideArgKwargNames(
                pyPamtraRadarSimulator.createRadarSpectrum)

            args = []
            for var in funcArgVars:
                args.append(mergedProfile[var])

            assert len(argNames) == len(args)

            input_core_dims = helpers.getInputCoreDims(args, ['sizeBin'])
            kwargs = {}
            for k in kwargNames:
                kwargs[k] = self.settings[k]

            nfft = kwargs['radarNFFT'] * (
                1 + 2*kwargs['radarAliasingNyquistInterv']
            )

            radarSpec = xr.apply_ufunc(
                pyPamtraRadarSimulator.createRadarSpectrum,
                *args,
                kwargs=kwargs,
                input_core_dims=input_core_dims,
                output_core_dims=[('dopplerVelocityAliased',)],
                output_dtypes=[mergedProfile.backscatterCrossSection.dtype],
                output_sizes={'dopplerVelocityAliased': nfft},
                dask='parallelized',
            )
            radarSpecs.append(radarSpec)
        radarSpecs = xr.concat(radarSpecs, dim='hydrometeor')
        radarSpecs = radarSpecs.sum('hydrometeor').unstack('merged')

        self.results['radarIdealizedSpectrum'] = radarSpecs
        self.results['radarIdealizedSpectrum'].attrs['unit'] = units.units[
            'radarIdealizedSpectrum']

        return self.results['radarIdealizedSpectrum']

    def _simulateRadar(self):

        argNames, kwargNames = helpers.provideArgKwargNames(
            pyPamtraRadarSimulator.simulateRadarSpectrum)

        kwargs = {}
        for k in kwargNames:
            kwargs[k] = self.settings[k]

        variables = [
            'height',
            'eddyDissipationRate',
            'horizontalWind',
            'radarIdealizedSpectrum',
            'pathIntegratedAtenuattion',
            'wavelength',
        ]

        mergedProfile = self.parent.profile.copy()
        mergedProfile['radarIdealizedSpectrum'] = self.results[
            'radarIdealizedSpectrum']

        mergedProfile = mergedProfile.stack(
            merged=helpers.concatDicts(
                self.parent.coords['additional'],
                self.parent.coords['layer'],
                self.parent.coords['frequency'])
        )

        args = []
        for var in variables:
            args.append(mergedProfile[var])

        assert len(argNames) == len(args)

        input_core_dims = helpers.getInputCoreDims(
            args, ['dopplerVelocityAliased'])

        radarSpec = xr.apply_ufunc(
            pyPamtraRadarSimulator.simulateRadarSpectrum,
            *args,
            kwargs=kwargs,
            input_core_dims=input_core_dims,
            output_core_dims=[('dopplerVelocity',)],
            output_dtypes=[mergedProfile.radarIdealizedSpectrum.dtype],
            output_sizes={'dopplerVelocity': 256},
            dask='parallelized',
        ).unstack('merged')

        self.results['radarSpectrum'] = radarSpec.assign_coords(
            dopplerVelocity=np.linspace(
                self.settings['radarMinV'],
                self.settings['radarMaxV'],
                self.settings['radarNFFT'],
                endpoint=False
            )
        )

        return self.results['radarSpectrum']

    def _calcMoments(self):

        output_core_dims = [
            ('peak'),
            ('peak'),
            ('peak'),
            ('peak'),
            ('peak'),
            ('peak'),
            ('peak'),
            ('peak'),
            ('peak'),
            tuple(),
            tuple(),
        ]

        output_names = [
            'radarReflectivity',
            'meanDopplerVel',
            'spectrumWidth',
            'skewness',
            'kurtosis',
            'leftSlope',
            'rightSlope',
            'leftEdge',
            'rightEdge',
            'quality',
            'noiseMean',
        ]

        output_sizes = {
            'peak': self.settings['momentsNPeaks'],
        }
        # theseVars
        input_core_dims = [['dopplerVelocity', ]]

        args = [
            self.results.radarSpectrum.stack(merged=helpers.concatDicts(
                self.parent.coords['additional'],
                self.parent.coords['layer'],
                self.parent.coords['frequency']
            )
            )
        ]

        # take care of settings
        argNames, kwargNames = helpers.provideArgKwargNames(
            pyPamtraRadarMoments.calc_radarMoments)
        kwargs = {}
        for k in kwargNames:
            kwargs[k] = self.settings[k]

        assert len(argNames) == len(args)

        moments = helpers.apply_ufunc_extended(
            _calc_radarMoments_wrapper,
            *args,
            kwargs=kwargs,
            output_names=output_names,
            input_core_dims=input_core_dims,
            output_core_dims=output_core_dims,
            output_dtypes=[args[0].dtype],
            output_sizes=output_sizes,
            dask='parallelized',
        )
        moments = moments.unstack('merged')

        moments = moments.assign_coords(
            peak=np.arange(1,self.settings['momentsNPeaks']+1)
            )

        self.results.merge(moments, inplace=True)
        return moments


def _calc_radarMoments_wrapper(*args, **kwargs):
    result = pyPamtraRadarMoments.calc_radarMoments(
        *args, **kwargs)
    spectrumOut, moments, slope, edge, quality, noiseMean = result

    moments[moments == -9999.] = np.nan
    slope[slope == -9999.] = np.nan
    edge[edge == -9999.] = np.nan
    noiseMean[noiseMean == -9999.] = np.nan

    radarReflectivity = 10 * np.log10(moments[:, 0, :])
    meanDopplerVel = moments[:, 1, :]
    spectrumWidth = moments[:, 2, :]
    skewness = moments[:, 3, :]
    kurtosis = moments[:, 4, :]
    leftSlope = slope[:, 0, :]
    rightSlope = slope[:, 1, :]
    leftEdge = edge[:, 0, :]
    rightEdge = edge[:, 1, :]
    quality = quality.astype(moments.dtype)
#     No = noiseMean

    result = (radarReflectivity, meanDopplerVel, spectrumWidth, skewness,
              kurtosis, leftSlope, rightSlope, leftEdge, rightEdge,
              quality, noiseMean)

    return result
