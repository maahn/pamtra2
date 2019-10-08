# -*- coding: utf-8 -*-

import inspect
import warnings

import numpy as np
import xarray as xr

from . import (aspectRatio, crossSectionArea, density, fallVelocity, mass,
               scattering, size, numberConcentration, relativePermittivity)
from .. import helpers, units
from ..libs import refractiveIndex


def DEFAULT_CALCULATION_ORDER():
    return [
        'sizeBounds',
        'sizeCenter',
        'sizeBoundsWidth',
        'aspectRatio',
        'density',
        'mass',
        'crossSectionArea',
        'numberConcentration',
        'relativePermittivityIce',
        'relativePermittivity',
        'scattering',
        'fallVelocity',
    ]


class hydrometeor(object):
    """generic class to store hydrometeor properties.

        Parameters
        ----------
        parent : pamtra2 object
            content of parent's class
        name : str, optional
            name of the hydrometeor
        discreteProperties : xr.Dataset, optional
            pre-calculated discrete properties
        calculationOrder : optional
            order for estimating the properties
        useFuncArgDefaults : bool, optional
            if parameters are not found in any kwarg, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        **kwargs :
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'numberConcentration'.

        Attributes
        ----------
        name : str, optional
            name of the hydrometeor
        index : int
            hydrometeor index in parent
        kind :  {'liquid', 'ice'}
            liquid or frozen hydrometeor?
        nBins : int
            number of size bins
        discreteProperties : xr.Dataset
            calculated discrete properties
        calculationOrder
            order for estimating the properties
        useFuncArgDefaults : bool
            if parameters are not found in any of funcArgs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        description : dict
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'numberConcentration'.

    """

    def __init__(
        self,
        parent=None,
        name=None,
        discreteProperties=None,
        calculationOrder=None,
        useFuncArgDefaults=True,
        **kwargs
    ):

        self.calculationOrder = calculationOrder
        # self.funcArgs = funcArgs
        self.useFuncArgDefaults = useFuncArgDefaults
        self.description = kwargs
        self.name = name
        # self.index = np.where(parent.profile.hydrometeor.values == name)[0][0]
        self._parentFull = parent
        # self.coords = parent.coords
        # self._parentFull.verbosity = parent.verbosity

        if discreteProperties is None:
            discreteProperties = xr.Dataset(
                coords=dict(
                    sizeBin=range(kwargs['nBins']),
                    sizeBin1=range(kwargs['nBins']+1),
                )
            )
        self.profile = discreteProperties

        return

    def getProfileAllBroadcasted(self, variables=None, sel={}):
        if variables is None:
            return xr.broadcast(self.profile.sel(**sel))[0]
        else:
            return xr.broadcast(self.profile.sel(**sel)[variables])[0]

    def getProfileWithParentAllBroadcasted(
        self,
        variables=None,
        parentVariables=None,
        exclude=[],
    ):
        profile = self.getProfileAllBroadcasted(variables)
        parent = self._parentFull.getProfileAllBroadcasted(
            parentVariables,
            sel={'hydrometeor': self.name},
        )
        profile, parent = xr.broadcast(profile, parent, exclude=exclude)
        merged = xr.merge((profile, parent))
        return merged

    @property
    def _parentProfile(self):
        """Helper function

        Returns
        -------
        _parentProfile
            Limited version of parent.profile. Contains only data belonging to
        the hydrometeor.
        """
        return self._parentFull.profile.sel(hydrometeor=self.name, drop=True)

    def _arrayOrFunc(self, thisKey, thisDesription, **fixedKwargs):
        """Helper function calling functions if required.

        Parameters
        ----------
        thisDesription :
            function or value or xr.DataArray
        **fixedKwargs :
            additional parameters for the function not defined elsewhere

        Returns
        -------
        thisProperty
            value or xr.DataArray
        """

        if callable(thisDesription):
            if self._parentFull.verbosity >= 1:
                print('callable')

            func = thisDesription

            try:
                self._keysToBeUsed.remove(thisKey)
            except ValueError:
                pass

            # inspect function to get the required arguments
            argspec = inspect.getargspec(func)
            funcArgs, funcVarargs, funcKeywords, funcDefaults = argspec

            if funcDefaults is None:
                funcDefaults = {}
            else:
                funcDefaults = dict(
                    zip(funcArgs[-len(funcDefaults):], funcDefaults))

            # where do we find the required data?
            kw4Func = {}
            for k in funcArgs:
                # if k in self.funcArgs.keys():
                #     kw4Func[k] = self.funcArgs[k]
                if k in self.profile.keys():
                    kw4Func[k] = self.profile[k]
                elif k in self._parentProfile.keys():
                    kw4Func[k] = self._parentProfile[k]
                elif k in self.description.keys():
                    kw4Func[k] = self.description[k]
                    try:
                        self._keysToBeUsed.remove(k)
                    except ValueError:
                        pass
                elif k in fixedKwargs.keys():
                    kw4Func[k] = fixedKwargs[k]
                elif self.useFuncArgDefaults and (k in funcDefaults.keys()):
                    kw4Func[k] = funcDefaults[k]
                else:
                    raise KeyError('Did not find %s in provided kwargs or '
                                   'discreteProperties or profile or '
                                   'functions\'s defaultArgs for '
                                   ' %s' % (k, thisKey))

            if self._parentFull.verbosity >= 2:
                print('kw4Func', kw4Func)
            thisProperty = func(**kw4Func)
        else:
            if self._parentFull.verbosity >= 1:
                print('not callable')
            thisProperty = thisDesription

        return thisProperty

    def solve(self):
        """Helper function to estimate all discrete properties of a
         hydrometeor

        Returns
        -------
        discreteProperties
            xr.Dataset with results
        """

        if self._parentFull is None:
            raise AttributeError('set .parent attribute to pamtra2 object')

        if self.calculationOrder is None:
            self.calculationOrder = DEFAULT_CALCULATION_ORDER()

        for key in self.description.keys():
            if key not in self.calculationOrder:
                self.calculationOrder.append(key)
            self.description.keys()

        self._keysToBeUsed = list(self.description.keys())

        for key in self.calculationOrder:

            if key not in self.description.keys():
                print('Did not find information about %s. This might cause'
                      ' trouble later.' % key)
                continue

            value = self.description[key]
            if self._parentFull.verbosity >= 1:
                print(key, value)

            thisProperty = self._arrayOrFunc(key, value)
            if (isinstance(thisProperty, xr.DataArray) and
                    (key in ['sizeCenter', 'sizeBoundsWidth'])):
                # when sizeCenter and sizeBoundsWidth are estimated from
                # sizeBounds with numpy, the dimension name is typically
                # wrong. So we try to rename:
                try:
                    thisProperty = thisProperty.rename({
                        'sizeBin1': 'sizeBin'
                    })
                except ValueError:
                    pass
            if not isinstance(thisProperty, xr.DataArray):
                if (key in ['sizeCenter', 'sizeBoundsWidth']):
                    # import pdb; pdb.set_trace()
                    thisProperty = xr.DataArray(
                        np.asarray(thisProperty.data),
                        coords=[self.profile.sizeBin],
                        attrs={'unit': units.units[key]},
                    )
                elif (key in ['sizeBounds']):
                    thisProperty = xr.DataArray(
                        np.asarray(thisProperty.data),
                        coords=[self.profile.sizeBin1],
                        attrs={'unit': units.units[key]},
                    )
            self.profile[key] = thisProperty

            self.profile[key].attrs.update(
                {'unit': units.units[key]}
            )
        self._postProcessing()

        self._keysToBeUsed = [x for x in self._keysToBeUsed if x not in
                              DEFAULT_CALCULATION_ORDER()]
        if len(self._keysToBeUsed) > 0:
            warnings.warn('The following kwargs were NOT used: '
                          '%s' % self._keysToBeUsed)

        # Apply units
        for k in self.profile.variables:
            self.profile[k].attrs['unit'] = units.units[k]

        # test results
        varsGreaterZero = [
            'sizeBoundsWidth',
            'sizeBounds',
            'sizeCenter',
            'mass',
            'aspectRatio',
            'density',
            'crossSectionArea',
            'relativePermittivity',
            'fallVelocity',
            'nBins',
            'extinctionCrossSection',
            'scatterCrossSection',
            'absorptionCrossSection',
            'backscatterCrossSection',
        ]
        varsGreaterEqualZero = [
            'numberConcentration',
        ]

        for key in varsGreaterZero:
            if key in self.profile.keys():
                assert (self.profile[key] > 0).all()

        for key in varsGreaterEqualZero:
            if key in self.profile.keys():
                assert (self.profile[key] >= 0).all()

        return self.profile

    def _postProcessing(self):
        """
        clean up 'artificial' dimensions. Required becaus xr.apply_ufunc cannot
        return multiple items.
        """
        variables = ['extinctionCrossSection', 'scatterCrossSection',
                     'absorptionCrossSection', 'backscatterCrossSection']
        scatteringProperty = helpers.dimensionToVariables(
            self.profile.scattering,
            'scatteringProperty',
            variables
        )

        for key in scatteringProperty.keys():
            scatteringProperty[key].attrs.update(
                {'unit': units.units[key]}
            )

        self.profile = self.profile.drop('scattering')
        self.profile = self.profile.merge(scatteringProperty)

        return


class softEllipsoidFixedDensity(hydrometeor):
    """hydrometeor class to be used for soft ellipsoids with fixed density."""

    def __init__(self, *args, **kwargs):
        if 'calculationOrder' not in kwargs.keys():
            kwargs['calculationOrder'] = DEFAULT_CALCULATION_ORDER()
        return super().__init__(*args, **kwargs)


class softEllipsoidMassSize(hydrometeor):
    """hydrometeor class to be used for soft ellipsoids with variable
    density depending on the mass-size relation."""

    def __init__(self, *args, **kwargs):
        if 'calculationOrder' not in kwargs.keys():
            kwargs['calculationOrder'] = helpers.swapListItems(
                DEFAULT_CALCULATION_ORDER(),
                'mass',
                'density'
            )

            ii = kwargs['calculationOrder'].index('numberConcentration')
            kwargs['calculationOrder'].insert(ii, 'relativePermittivityIce')

        return super().__init__(*args, **kwargs)


class cloud(softEllipsoidFixedDensity):
    """hydrometeor class for cloud droplets. 

        Parameters
        ----------
        parent : pamtra2 object
            content of parent's class
        name : str, optional
            name of the hydrometeor
        nBins : int, optional
            number of bins. Note taht the spectral radar simulator requires
            >= 2 hydrometeor bins. (default 2)
        numberConcentration : func, optional
            Function to describe the size distribution. default:
            numberConcentration.monoDisperseWC which requires hydrometeorContent
            either as kwarg or in pamtra2's profile
        scattering : func, optional
            Function to describe the single scattering of the hydrometeor.
            default: scattering.Mie
        Dmin : float, optional
            Smallest boundary of size distribution. Note that Dmin != Dmax
            is required also for a monodisperse distribution. (default
            1e-5 - 1e-10)
        Dmax :   float, optional
            Largest boundary of size distribution. Note that Dmin != Dmax
            is required also for a monodisperse distribution. (default
            1e-5 + 1e-10)
        useFuncArgDefaults : bool, optional
            if parameters are not found in any of kwargs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        **kwargs :
            Further properties if required. Also the parameters defined in
            defaultArgs (see source) can be overwritten if provided.

        Attributes
        ----------
        name : str, optional
            name of the hydrometeor
        index : int
            hydrometeor index in parent
        kind :  {'liquid', 'ice'}
            liquid or frozen hydrometeor?
        nBins : int
            number of size bins
        discreteProperties : xr.Dataset
            calculated discrete properties
        calculationOrder
            order for estimating the properties
        useFuncArgDefaults : bool
            if parameters are not found in any of funcArgs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        description : dict
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'numberConcentration'.

    """

    def __init__(
        self,
        *args,
        nBins=2,
        numberConcentration=numberConcentration.monoDisperseWC,
        scattering=scattering.Mie,
        Dmin=1e-5 - 1e-10,
        Dmax=1e-5 + 1e-10,
        **kwargs
    ):

        defaultArgs = {}
        defaultArgs['aspectRatio'] = 1.0
        defaultArgs['mass'] = mass.ellipsoid
        defaultArgs['density'] = density.water
        defaultArgs['crossSectionArea'] = crossSectionArea.sphere
        defaultArgs['sizeCenter'] = size.boundsToMid
        defaultArgs['sizeBounds'] = size.linspaceBounds
        defaultArgs['sizeBoundsWidth'] = size.boundsWidth
        defaultArgs['fallVelocity'] = fallVelocity.khvorostyanov01_drops
        defaultArgs['relativePermittivity'] = relativePermittivity.\
            water_turner_kneifel_cadeddu

        kwargs['nBins'] = nBins
        kwargs['numberConcentration'] = numberConcentration
        kwargs['scattering'] = scattering
        kwargs['Dmin'] = Dmin
        kwargs['Dmax'] = Dmax

        defaultArgs.update(kwargs)

        return super().__init__(*args, **defaultArgs)


class rain(softEllipsoidFixedDensity):
    """hydrometeor class for rain.

        Parameters
        ----------
        parent : pamtra2 object
            content of parent's class
        name : str, optional
            name of the hydrometeor
        nBins : int, optional
            number of bins. Note taht the spectral radar simulator requires
            >= 2 hydrometeor bins. (default 50)
        numberConcentration : func, optional
            Function to describe the size distribution. default:
            numberConcentration.exponentialN0WC which requires hydrometeorContent
            either as kwarg or in pamtra2's profile. In addition, N0 has to be
            provided
        scattering : func, optional
            Function to describe the single scattering of the hydrometeor.
            default: scattering.Mie
        aspectRatio : float, optional
            hydrometeor aspect ratio (default 1.0)
        Dmin : float, optional
            Smallest boundary of size distribution. Note that Dmin != Dmax
            is required also for a monodisperse distribution. (default
            1e-6)
        Dmax :   float, optional
            Largest boundary of size distribution. Note that Dmin != Dmax
            is required also for a monodisperse distribution. (default
            0.008)
        N0 : float, optional
            N0 of exponential distribution. (default 8e6)
        useFuncArgDefaults : bool, optional
            if parameters are not found in any of kwargs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        **kwargs :
            Further properties if required. Also the parameters defined in
            defaultArgs (see source) can be overwritten if provided.

        Attributes
        ----------
        name : str, optional
            name of the hydrometeor
        index : int
            hydrometeor index in parent
        kind :  {'liquid', 'ice'}
            liquid or frozen hydrometeor?
        nBins : int
            number of size bins
        discreteProperties : xr.Dataset
            calculated discrete properties
        calculationOrder
            order for estimating the properties
        useFuncArgDefaults : bool
            if parameters are not found in any of funcArgs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        description : dict
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'numberConcentration'.

    """

    def __init__(
        self,
        *args,
        nBins=50,
        numberConcentration=numberConcentration.exponentialN0WC,
        scattering=scattering.Mie,
        aspectRatio=1.0,  # TODO add non 1 values
        Dmin=1e-6,
        Dmax=0.008,
        N0=8.e6,  # from COSMO 1 moments scheme...
        **kwargs
    ):

        defaultArgs = {}
        defaultArgs['mass'] = mass.ellipsoid
        defaultArgs['density'] = density.water
        defaultArgs['crossSectionArea'] = crossSectionArea.sphere  # TODO add
        # elliptical
        defaultArgs['sizeCenter'] = size.boundsToMid  # TODO add log version
        defaultArgs['sizeBounds'] = size.logspaceBounds
        defaultArgs['sizeBoundsWidth'] = size.boundsWidth
        defaultArgs['fallVelocity'] = fallVelocity.khvorostyanov01_drops
        defaultArgs['relativePermittivity'] = relativePermittivity.\
            water_turner_kneifel_cadeddu

        kwargs['nBins'] = nBins
        kwargs['numberConcentration'] = numberConcentration
        kwargs['scattering'] = scattering
        kwargs['aspectRatio'] = aspectRatio
        kwargs['Dmin'] = Dmin
        kwargs['Dmax'] = Dmax

        defaultArgs.update(kwargs)

        return super().__init__(*args, **defaultArgs)


class ice(softEllipsoidFixedDensity):
    """hydrometeor class for ice.

        Parameters
        ----------
        parent : pamtra2 object
            content of parent's class
        name : str, optional
            name of the hydrometeor
        nBins : int, optional
            number of bins. Note taht the spectral radar simulator requires
            >= 2 hydrometeor bins. (default 2)
        numberConcentration : func, optional
            Function to describe the size distribution. default:
            numberConcentration.monoDisperseWC which requires hydrometeorContent
            either as kwarg or in pamtra2's profile
        scattering : func, optional
            Function to describe the single scattering of the hydrometeor.
            default: scattering.Mie
        Dmin : float, optional
            Smallest boundary of size distribution. Note that Dmin != Dmax
            is required also for a monodisperse distribution. (default
            1e-5 - 1e-10)
        Dmax :   float, optional
            Largest boundary of size distribution. Note that Dmin != Dmax
            is required also for a monodisperse distribution. (default
            1e-5 + 1e-10)
        useFuncArgDefaults : bool, optional
            if parameters are not found in any of kwargs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        **kwargs :
            Further properties if required. Also the parameters defined in
            defaultArgs (see source) can be overwritten if provided.

        Attributes
        ----------
        name : str, optional
            name of the hydrometeor
        index : int
            hydrometeor index in parent
        kind :  {'liquid', 'ice'}
            liquid or frozen hydrometeor?
        nBins : int
            number of size bins
        discreteProperties : xr.Dataset
            calculated discrete properties
        calculationOrder
            order for estimating the properties
        useFuncArgDefaults : bool
            if parameters are not found in any of funcArgs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        description : dict
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'numberConcentration'.

    """

    def __init__(
        self,
        *args,
        nBins=2,
        numberConcentration=numberConcentration.monoDisperseWC,
        scattering=scattering.Mie,
        Dmin=1e-4 - 1e-10,
        Dmax=1e-4 + 1e-10,
        **kwargs
    ):

        defaultArgs = {}
        defaultArgs['aspectRatio'] = 1.0
        defaultArgs['mass'] = mass.ellipsoid
        defaultArgs['density'] = density.ice
        defaultArgs['crossSectionArea'] = crossSectionArea.sphere
        defaultArgs['sizeCenter'] = size.boundsToMid
        defaultArgs['sizeBounds'] = size.linspaceBounds
        defaultArgs['sizeBoundsWidth'] = size.boundsWidth
        defaultArgs['fallVelocity'] = fallVelocity.heymsfield10_particles
        defaultArgs['relativePermittivity'] = relativePermittivity.\
            ice_matzler_2006

        kwargs['nBins'] = nBins
        kwargs['numberConcentration'] = numberConcentration
        kwargs['scattering'] = scattering
        kwargs['Dmin'] = Dmin
        kwargs['Dmax'] = Dmax

        defaultArgs.update(kwargs)

        return super().__init__(*args, **defaultArgs)


class snow(softEllipsoidMassSize):
    """hydrometeor class for snow.  """

    def __init__(
        self,
        *args,
        nBins=2,
        aspectRatio=1.0,  # TODO: implement 0.6
        numberConcentration=numberConcentration.exponentialN0WC,  # TODO:change to WC
        scattering=scattering.Mie,
        Dmin=0.5e-8,
        Dmax=1.e-2,
        N0=1e5,
        massSizeA=0.0121,
        massSizeB=1.9,
        areaSizeA=0.4,
        areaSizeB=1.8,
        minDensity=100,
        **kwargs
    ):

        defaultArgs = {}
        defaultArgs['mass'] = mass.powerLaw
        defaultArgs['crossSectionArea'] = crossSectionArea.powerLaw
        defaultArgs['density'] = density.softEllipsoid,
        defaultArgs['sizeCenter'] = size.boundsToMid
        defaultArgs['sizeBounds'] = size.logspaceBounds
        defaultArgs['sizeBoundsWidth'] = size.boundsWidth
        defaultArgs['fallVelocity'] = fallVelocity.heymsfield10_particles
        defaultArgs['relativePermittivityIce'] = relativePermittivity.\
            ice_matzler_2006
        defaultArgs['relativePermittivity'] = relativePermittivity.\
            mixing_sihvola

        kwargs['nBins'] = nBins
        kwargs['aspectRatio'] = aspectRatio
        kwargs['numberConcentration'] = numberConcentration
        kwargs['scattering'] = scattering
        kwargs['Dmin'] = Dmin
        kwargs['Dmax'] = Dmax
        kwargs['N0'] = N0
        kwargs['massSizeA'] = massSizeA
        kwargs['massSizeB'] = massSizeB
        kwargs['areaSizeA'] = areaSizeA
        kwargs['areaSizeB'] = areaSizeB
        kwargs['minDensity'] = minDensity

        defaultArgs.update(kwargs)

        return super().__init__(*args, **defaultArgs)
