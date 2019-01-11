# -*- coding: utf-8 -*-
import warnings
from collections import OrderedDict

import numpy as np
import xarray as xr

from . import constants, dimensions, helpers, units
from .libs import meteo_si

__version__ = 0.2


class customProfile (xr.Dataset):
    """ """

    def __init__(
        self,
        parent,
        profileVars=[],
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

            return


class pamtra2(object):
    """ """

    def __init__(
        self,
        nLayer,
        hydrometeors,
        frequencies,
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
            'hydrometeorContent',
                [dimensions.ADDITIONAL, dimensions.LAYER,
                 dimensions.HYDROMETEOR],
                np.float64
            ),
        ],
        profile=None,
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

        # remove length zero coordinates
        for k1 in self.coords.keys():
            for k2 in list(self.coords[k1].keys()):
                if len(self.coords[k1][k2]) == 0:
                    warnings.warn('Dimension %s has length 0 and was removed'
                                  % k2)
                    del self.coords[k1][k2]

        if profile is None:
            self.profile = customProfile(
                self,
                profileVars=profileVars,
            )
        else:
            self.profile = profile

        self.profile['wavelength'] = constants.speedOfLight /\
            self.profile.frequency

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
                numberConcentration = self.hydrometeors[
                    name].profile.numberConcentration.fillna(0)
                crossSec = self.hydrometeors[name].profile[crossSection]
                if frequencies is not None:
                    crossSec.sel(frequency=frequencies)
                thisHydro = crossSec * numberConcentration
                perHydro.append(thisHydro)
            integrated[crossSection] = xr.concat(
                perHydro, dim='hydro').sum(['hydro', 'sizeBin'])
        return integrated

    def addMissingVariables(self):

        self.addHeightBinDepth()
        self.addSpecificHumidity()
        self.addAbsoluteHumidity()
        self.addDryAirDensity()
        self.addAirDensity()
        self.addDynamicViscosity()
        self.addKinematicViscosity()
        self.addWaterVaporPressure()

        return self.profile

    def addHeightBinDepth(self, update=False):

        if (not update) and ('heightBinDepth' in self.profile.keys()):
            return
        else:
            self.profile['heightBinDepth'] = helpers.xrGradient(
                self.profile.height, 'layer')

            return self.profile['heightBinDepth']

    def addAbsoluteHumidity(self):
        '''
        add absolute humidity
        '''

        args = (
            self.profile.relativeHumidity/100.,
            self.profile.temperature
        )

        self.profile['absoluteHumidity'] = xr.apply_ufunc(
                meteo_si.humidity.rh2a,
                *args,
                kwargs={},
                output_dtypes=[np.float64],
                dask='parallelized',
            )

        return self.profile['absoluteHumidity']

    def addSpecificHumidity(self):
        '''
        add specific humidity
        '''

        args = (
            self.profile.relativeHumidity/100.,
            self.profile.temperature,
            self.profile.pressure,
        )

        self.profile['specificHumidity'] = xr.apply_ufunc(
                meteo_si.humidity.rh2q,
                *args,
                kwargs={},
                output_dtypes=[np.float64],
                dask='parallelized',
            )

        return self.profile['specificHumidity']

    def addWaterVaporPressure(self):
        '''
        add waterVaporPressure
        '''

        args = (
            self.profile['specificHumidity'],
            self.profile.pressure
        )

        self.profile['waterVaporPressure'] = xr.apply_ufunc(
                meteo_si.humidity.q2e,
                *args,
                kwargs={},
                output_dtypes=[np.float64],
                dask='parallelized',
            )

        return self.profile['waterVaporPressure']

    def addDryAirDensity(self):
        '''
        add dry air density
        '''

        p = self.profile.pressure
        T = self.profile.temperature
        q = 0

        args = (p, T, q)
        self.profile['dryAirDensity'] = xr.apply_ufunc(
                meteo_si.density.moist_rho_q,
                *args,
                kwargs={},
                output_dtypes=[np.float64],
                dask='parallelized',
            )
        return self.profile['dryAirDensity']

    def addAirDensity(self):
        '''
        add air density
        '''

        if 'specificHumidity' not in self.profile.keys():
            self.addSpecificHumidity()

        p = self.profile.pressure
        T = self.profile.temperature
        q = self.profile.specificHumidity
        try:
            qm = self.profile.hydrometeorContent.sum('hydrometeor')
        except ValueError:
            warnings.warn('hydrometeor content not considered when calculating'
                          ' air density, because hydrometeor dimension is '
                          'missing. ')
            qm = 0
        except AttributeError:
            warnings.warn('hydrometeor content not considered when calculating'
                          ' air density, because hydrometeorContent variable'
                          ' is missing.')
            qm = 0

        args = (p, T, q, qm)
        self.profile['airDensity'] = xr.apply_ufunc(
                meteo_si.density.moist_rho_q,
                *args,
                kwargs={},
                output_dtypes=[np.float64],
                dask='parallelized',
            )

        return self.profile['airDensity']

    def addDynamicViscosity(self):
        '''
        dynamic viscosity of dry air
        '''

        self.profile['dynamicViscosity'] = xr.apply_ufunc(
                _dynamic_viscosity_air,
                self.profile.temperature,
                kwargs={},
                output_dtypes=[np.float64],
                dask='parallelized',
            )

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

    def addHydrometeor(
        self,
        hydrometeorClass,
        solve=True
    ):
        """
        Add hydrometeor properties for one hydrometeor. Hydrometeor is added
        to self.hydrometeors[name]

        Parameters
        ----------
        hydrometeorClass :
            initialized hydrometeor class
        solve : bool, optional
            Solve hydrometeor description (default true)

        Returns
        -------
        hydrometeorClass :
            linked hydrometeor class
        """

        name = hydrometeorClass.name
        self.hydrometeors[name] = hydrometeorClass
        # link parent object with hydrometeor
        self.hydrometeors[name]._parentFull = self
        if solve:
            self.hydrometeors[name].solve()

        return self.hydrometeors[name]


    def addInstrument(
        self,
        instrumentClass,
        solve=True
    ):
        """

        Parameters
        ----------
        instrumentClass :
            initialized instrument class
        solve : bool, optional
            Solve instrument description (default true)

        Returns
        -------

        """

        name = instrumentClass.name
        self.instruments[name] = instrumentClass

        # link parent object with instruments
        self.instruments[name].parent = self

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
