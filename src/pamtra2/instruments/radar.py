# -*- coding: utf-8 -*-
import warnings

import numpy as np
import xarray as xr

from .. import helpers, units
from .core import microwaveInstrument


class simpleRadar(microwaveInstrument):
    def __init__(
        self,
        name='simpleRadar',
        parent=None,
        frequencies='all',
        radarK2=0.93,
        gaseousAttenuationModel='Rosenkranz98',
        applyAttenuation=None,
        **kwargs
    ):
        super().__init__(
            name=name,
            parent=parent,
            frequencies=frequencies,
            radarK2=radarK2,
            applyAttenuation=applyAttenuation,
            gaseousAttenuationModel=gaseousAttenuationModel,
            **kwargs
        )

    def solve(self):

        self._link_parent()

        self._calcPIA()

        if self.parent.nHydrometeors > 0:

            wavelength = self.parent.profile.wavelength.sel(
                frequency=self.frequencies)
            K2 = self.settings['radarK2']

            Ze_increment = 0.
            MDV_increment = 0.

            for hydro in self.hydrometeorProfiles.keys():
                bsc = self.hydrometeorProfiles[hydro].backscatterCrossSection\
                    .sel(frequency=self.frequencies)
                N = self.hydrometeorProfiles[hydro]\
                    .numberConcentration.fillna(0)
                Ze_bin = 1.e18*(1.e0/(K2*np.pi**5))*bsc*N*(wavelength)**4
                vel = self.hydrometeorProfiles[hydro].fallVelocity
                Ze_increment += Ze_bin.sum(['sizeBin'])
                MDV_increment += (vel*Ze_bin).sum(['sizeBin'])
            # perHydro = [] 
            # for name in self.hydrometeors.keys():
            #     numberConcentration = self.hydrometeors[
            #         name].profile.numberConcentration.fillna(0)
            #     crossSec = self.hydrometeors[name].profile[crossSection]
            #     if frequencies is not None:
            #         crossSec.sel(frequency=frequencies)
            #     thisHydro = crossSec * numberConcentration
            #     perHydro.append(thisHydro)
            # integrated[crossSection] = xr.concat(
            #     perHydro, dim='hydro').sum(['hydro', 'sizeBin'])

            # crossSections = self.parent.getIntegratedScatteringCrossSections(
            #     crossSections=['backscatterCrossSection'],
            #     frequencies=self.frequencies,
            # )

            # back = crossSections['backscatterCrossSection']
            # Ze_back = 1.e18*(1.e0/(K2*np.pi**5))*back*(wavelength)**4

            self.results['radarReflectivity'] = 10*np.log10(Ze_increment)
            self.results['meanDopplerVel'] = MDV_increment/Ze_increment

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

        self.results['specificAttenuation'] = 0
        self.results['pathIntegratedAttBottomUp'] = 0
        self.results['pathIntegratedAttTopDown'] = 0

    def _attenuation2pia(self, attenuation, dim='layer'):

        attenuation = attenuation * self.parent.profile.heightBinDepth

        attenuationRev = attenuation.isel(
            layer=attenuation.coords[dim].values[::-1])

        PIA_bottomup = attenuation.cumsum(dim) * 2 - attenuation
        PIA_topdown = attenuationRev.cumsum(dim) * 2 - attenuationRev

        PIA_topdown = PIA_topdown.isel(
            layer=attenuation.coords[dim].values[::-1])
        return PIA_bottomup, PIA_topdown

