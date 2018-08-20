# -*- coding: utf-8 -*-
import warnings

import numpy as np
import xarray as xr

from .. import helpers, units
from ..libs import pyPamtraRadarMoments, pyPamtraRadarSimulator
from .core import microwaveInstrument


class simpleRadar(microwaveInstrument):
    def __init__(
        self,
        parent,
        frequencies='all',
        radarK2=0.93,
        gaseousAttenuationModel='Rosenkranz98',
        applyAttenuation=None,
        **kwargs
    ):
        super().__init__(
            parent,
            frequencies=frequencies,
            radarK2=radarK2,
            applyAttenuation=applyAttenuation,
            gaseousAttenuationModel=gaseousAttenuationModel,
            **kwargs
        )

    def solve(self):

        self._calcPIA()

        if self.parent.nHydrometeors > 0:

            crossSections = self.parent.getIntegratedScatteringCrossSections(
                crossSections=['backscatterCrossSection'],
                frequencies=self.frequencies,
            )

            back = crossSections['backscatterCrossSection']
            wavelength = self.parent.profile.wavelength.sel(
                frequency=self.frequencies)
            K2 = self.settings['radarK2']

            Ze_back = 1.e18*(1.e0/(K2*np.pi**5))*back*(wavelength)**4

            self.results['radarReflectivity'] = 10*np.log10(Ze_back)

            if self.settings['applyAttenuation'] is None:
                pass
            elif self.settings['applyAttenuation'] == 'bottomUp':
                self.results['radarReflectivity'] -= self.results[
                    'pathIntegratedAttBottomUp']
            elif self.settings['applyAttenuation'] == 'topDown':
                self.results['radarReflectivity'] -= self.results[
                    'pathIntegratedAttTopDown']
            else:
                raise ValueError('Do not understand applyAttenuation: %s. Must'
                                 ' be None, "bottomUp" or "topDown"' %
                                 self.settings['applyAttenuation'])
        else:
            warnings.warn('No hydrometeors found. Calculating only gaseous '
                          'attenuation.')

        for vv in self.results.variables:
            self.results[vv].attrs['unit'] = units.units[vv]

        return self.results

    def _calcPIA(self):
        absHydro = self._calcHydrometeorAbsorption()
        absGas = self._calcGaseousAbsorption()

        # Doviak Zrnic eq 3.12 with 10*np.log10(np.exp(X) = 4.34*X
        specificAttenuation = 10*np.log10(np.exp(absHydro + absGas))

        PIA_bottomup, PIA_topdown = self._attenuation2pia(specificAttenuation)

        self.results['specificAttenuation'] = specificAttenuation
        self.results['pathIntegratedAttBottomUp'] = PIA_bottomup
        self.results['pathIntegratedAttTopDown'] = PIA_topdown

    def _attenuation2pia(self, attenuation, dim='layer'):

        attenuation = attenuation * self.parent.profile.heightBinDepth

        attenuationRev = attenuation.isel(
            layer=attenuation.coords[dim].values[::-1])

        PIA_bottomup = attenuation.cumsum(dim) * 2 - attenuation
        PIA_topdown = attenuationRev.cumsum(dim) * 2 - attenuationRev

        PIA_topdown = PIA_topdown.isel(
            layer=attenuation.coords[dim].values[::-1])
        return PIA_bottomup, PIA_topdown


class dopplerRadarPamtra(simpleRadar):

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
        applyAttenuation=None,
        gaseousAttenuationModel='Rosenkranz98',
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
            applyAttenuation=applyAttenuation,
            gaseousAttenuationModel=gaseousAttenuationModel,
        )

        if len(frequencies) > 1:
            warnings.warn('Note that the radar simulator '
                          'assumes that settings are the same for '
                          'all frequencies (due to '
                          'missing support for vectorized '
                          'settings.)')

    def solve(self):

        for name in self.hydrometeorProfiles.keys():
            if len(self.hydrometeorProfiles[name].sizeBin) <= 1:
                raise ValueError('nbins must be greater than 1 for the'
                                 'spectral radar simulator!')

        self._calcPIA()
        if self.parent.nHydrometeors > 0:

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
                mergedProfile['sizeDistribution'].fillna(0)
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
        radarSpecs = radarSpecs.sum('hydrometeor')
        radarSpecs = xr.Dataset({'radarSpecs': radarSpecs})
        radarSpecs = helpers.xrFastUnstack(radarSpecs, 'merged').radarSpecs

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
            'pathIntegratedAttenuation',
            'wavelength',
        ]

        mergedProfile = self.parent.profile.copy()
        mergedProfile['radarIdealizedSpectrum'] = self.results[
            'radarIdealizedSpectrum']

        if self.settings['applyAttenuation'] is None:
            mergedProfile['pathIntegratedAttenuation'] = xr.zeros_like(
                mergedProfile.height)
        elif self.settings['applyAttenuation'] == 'bottomUp':
            mergedProfile['pathIntegratedAttenuation'] = self.results[
                'pathIntegratedAttBottomUp']
        elif self.settings['applyAttenuation'] == 'topDown':
            mergedProfile['pathIntegratedAttenuation'] = self.results[
                'pathIntegratedAttTopDown']
        else:
            raise ValueError('Do not understand applyAttenuation: %s. Must be'
                             'None, "bottomUp" or "topDown"' %
                             self.settings['applyAttenuation'])

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
        )
        radarSpec = xr.Dataset({'radarSpec': radarSpec})
        radarSpec = helpers.xrFastUnstack(radarSpec, 'merged').radarSpec

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
        moments = helpers.xrFastUnstack(moments, 'merged')

        moments = moments.assign_coords(
            peak=np.arange(1, self.settings['momentsNPeaks']+1)
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
