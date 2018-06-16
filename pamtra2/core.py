# -*- coding: utf-8 -*-

from collections import OrderedDict
import numpy as np
import xarray as xr

import meteo_si

from . import helpers
from . import units
from . import dimensions
from . import constants


__version__ = 0.2


class customProfile (xr.Dataset):
    """ """

    def __init__(
        self,
        parent,
        profileVars=[],
        nLayer=None,
        hydrometeors=None,
        frequencies=None,
        additionalDims={},
    ):

        if isinstance(profileVars, xr.Dataset):
            self = profileVars
            return
        else:
            super().__init__(coords=parent.coords['all'])

            for var, coord, dtype in profileVars:
                coord = helpers.concatDicts(
                    *map(lambda x: parent.coords[x], coord)
                )
                thisShape = tuple(map(len, coord.values()))
                self[var] = xr.DataArray(
                    (np.zeros(thisShape)*np.nan).astype(dtype),
                    coords=coord.values(),
                    dims=coord.keys(),
                    attrs={'unit': units.units[var]},
                )

            self['wavelength'] = constants.speedOfLight/self.frequency

            return


class pamtra2(object):
    """ """

    def __init__(
        self,
        profileVars=[
            (
            'height',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'temperature',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'pressure',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'relativeHumidity',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'horizontalWind',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'verticalWind',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'eddyDissipationRate',
                [dimensions.ADDITIONAL, dimensions.LAYER],
                np.float64
            ),
            (
            'waterContent',
                [dimensions.ADDITIONAL, dimensions.LAYER, dimensions.HYDROMETEOR],
                np.float64
            ),
        ],
        nLayer=None,
        hydrometeors=None,
        frequencies=None,
        additionalDims={},
    ):

        self.coords = {}
        self.coords[dimensions.ADDITIONAL] = OrderedDict(additionalDims)
        self.coords[dimensions.LAYER] = OrderedDict(layer=range(nLayer))
        self.coords[dimensions.HYDROMETEOR] = OrderedDict(
            hydrometeor=hydrometeors)
        self.coords[dimensions.FREQUENCY] = OrderedDict(
            frequency=frequencies)
        self.coords['nonCore'] = helpers.concatDicts(
            self.coords[dimensions.ADDITIONAL],
            self.coords[dimensions.LAYER],
        )
        self.coords['all'] = helpers.concatDicts(
            self.coords[dimensions.ADDITIONAL],
            self.coords[dimensions.LAYER],
            self.coords[dimensions.HYDROMETEOR],
            self.coords[dimensions.FREQUENCY],
        )

        self.profile = customProfile(
            self,
            profileVars=profileVars,
            nLayer=nLayer,
            hydrometeors=hydrometeors,
            frequencies=frequencies,
            additionalDims=additionalDims
        )
        self.additionalDims = additionalDims
        self.hydrometeors = helpers.AttrDict()
        for hh in hydrometeors:
            self.hydrometeors[hh] = None

        self.instruments = helpers.AttrDict()

        return

    def getProfileAllBroadcasted(self, variables=None, sel={}):
        if variables is None:
            return xr.broadcast(self.profile.sel(**sel))[0]
        else:
            return xr.broadcast(self.profile.sel(**sel)[variables])[0]

    def getIntegratedScatteringCrossSections(
        self,
        frequencies=None,
        crossSections='all'
    ):

        integrated = {}

        if crossSections == 'all':
            crossSections = [
                'extinctionCrossSection'
                'scatterCrossSection'
                'absorptionCrossSection',
                'backscatterCrossSection',
            ]
        for crossSection in crossSections:
            perHydro = []
            for name in self.hydrometeors.keys():
                sizeDistribution = self.hydrometeors[
                    name].profile.sizeDistribution
                sizeWidth = self.hydrometeors[
                    name].profile.sizeBoundsWidth
                crossSec = self.hydrometeors[name].profile[crossSection]
                if frequencies is not None:
                    crossSec.sel(frequency=frequencies)
                thisHydro = crossSec * sizeWidth * sizeDistribution
                perHydro.append(thisHydro)
            integrated[crossSection] = xr.concat(
                perHydro, dim='hydro').sum(['hydro', 'sizeBin'])
        return integrated

    def addMissingVariables(self):

        self.addHeightBinDepth()
        self.addDryAirDensity()
        self.addAirDensity()
        self.addDynamicViscosity()
        self.addKinematicViscosity()
        self.addSpecificHumidity()
        self.addAbsoluteHumidity()
        self.addWaterVaporPressure()

        return self.profile

    def addHeightBinDepth(self, update=False):

        if (not update) and ('heightBinDepth' in self.profile.keys()):
            return
        else:
            self.profile['heightBinDepth'] = helpers.xrGradient(
                self.parent.profile.height, 'layer')

            return self.profile['heightBinDepth']

    def addAbsoluteHumidity(self):
        '''
        add absolute humidity
        '''

        self.profile['absoluteHumidity'] = meteo_si.humidity.rh2a(
            self.profile.relativeHumidity/100.,
            self.profile.temperature
        )

        return self.profile['absoluteHumidity']

    def addSpecificHumidity(self):
        '''
        add specific humidity
        '''

        self.profile['specificHumidity'] = meteo_si.humidity.rh2q(
            self.profile.relativeHumidity/100.,
            self.profile.temperature,
            self.profile.pressure,
        )

        return self.profile['specificHumidity']

    def addWaterVaporPressure(self):
        '''
        add waterVaporPressure
        '''

        self.profile['waterVaporPressure'] = meteo_si.humidity.q2e(
            self.profile['specificHumidity'],
            self.profile.pressure
        )

        return self.profile['waterVaporPressure']

    def addDryAirDensity(self):
        '''
        add dry air density
        '''

        p = self.profile.pressure
        T = self.profile.temperature
        rh = self.profile.temperature/100.

        self.profile['dryAirDensity'] = meteo_si.density.moist_rho_rh(p, T, rh)

        return self.profile['dryAirDensity']

    def addAirDensity(self):
        '''
        add air density
        '''

        p = self.profile.pressure
        T = self.profile.temperature
        rh = self.profile.temperature/100.
        qm = self.profile.waterContent.sum('hydrometeor')

        self.profile['airDensity'] = meteo_si.density.moist_rho_rh(
            p, T, rh, qm)

        return self.profile['airDensity']

    def addDynamicViscosity(self):
        '''
        dynamic viscosity of dry air
        '''

        self.profile['dynamicViscosity'] = _dynamic_viscosity_air(
            self.profile.temperature)

        return self.profile['dynamicViscosity']

    def addKinematicViscosity(self):
        '''
        kinematic viscosity of dry air
        '''

        if 'dryAirDensity' not in self.profile.keys():
            self.addAirDensity()
        if 'dynamicViscosity' not in self.profile.keys():
            self.addDynamicViscosity()
        self.profile['kinematicViscosity'] = (
            self.profile['dynamicViscosity']/self.profile['dryAirDensity']
        )

        return self.profile['kinematicViscosity']

    @property
    def nHydrometeors(self):
        """ """
        return len(self.hydrometeors)

    @property
    def nInstruments(self):
        """ """
        return len(self.instruments)

    @property
    def nLayer(self):
        """ """
        return len(self.profile.layer)

    def describeHydrometeor(
        self,
        hydrometeorClass,
        solve=True,
        **kwargs,
    ):
        """
        Add hydrometeor properties for one hydrometeor. Hydrometeor is added
        to self.hydrometeors[name]

        Parameters
        ----------
        hydrometeorClass :
            hydrometeor class
        **kwargs : dict
            Arguments handed over to hydrometeorClass to initialize class

        Returns
        -------
        hydrometeorClass :
            Evaluated hydrometeor class
        """

        name = kwargs['name']
        self.hydrometeors[name] = hydrometeorClass(
            self,
            **kwargs
        )

        if solve:
            self.hydrometeors[name].solve()

        return self.hydrometeors[name]

    def addInstrument(
        self,
        instrumentClass,
        name=None,
        frequencies=[],
        solve=True,
        **kwargs,
    ):
        """

        Parameters
        ----------
        name :

        frequencies :
             (Default value = [])

        Returns
        -------

        """

        if not hasattr(frequencies, '__iter__'):
            frequencies = [frequencies]

        self.instruments[name] = instrumentClass(
            self,
            frequencies=frequencies,
            **kwargs
        )

        if solve:
            self.instruments[name].solve()

        return self.instruments[name]


# MOVE TO METEOSI!
def _dynamic_viscosity_air(temperature):
    """
    ! This function returns the dynamic viscosity of dry air in Pa s
    ! Sutherland law
    ! coefficients from F. M. White, Viscous Fluid Flow, 2nd ed., McGraw-Hill,
    ! (1991). Kim et al., arXiv:physics/0410237v1
    """

    mu0 = 1.716e-5  # Pas
    T0 = 273.
    C = 111.  # K

    eta = mu0*((T0 + C)/(temperature + C))*(temperature/T0)**1.5

    return eta


# # MOVE TO METEOSI!
# def _kinematic_viscosity_air(temperature, dryAirDensity):
#     # ! This function returns the kineamtic viscosity_air

#     viscosity = _dynamic_viscosity_air(temperature)
#     return viscosity/dryAirDensity
