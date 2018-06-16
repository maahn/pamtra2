# -*- coding: utf-8 -*-

import warnings
import numpy as np
import xarray as xr
import inspect

import refractiveIndex

from .. import units
from .. import helpers

from . import aspectRatio
from . import crossSectionArea
from . import density
from . import fallVelocity
from . import mass
from . import scattering
from . import size
from . import sizeDistribution


def DEFAULT_CALCULATION_ORDER():
    return [
        'sizeBounds',
        'sizeCenter',
        'sizeBoundsWidth',
        'aspectRatio',
        'density',
        'mass',
        'crossSectionArea',
        'sizeDistribution',
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
        funcArgs : dict, optional
            additional arguments for the functions describing the hydrometeor
            default {}.
        useFuncArgDefaults : bool, optional
            if parameters are not found in any of funcArgs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        **kwargs :
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'sizeDistribution'.

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
        funcArgs : dict
            additional arguments for the functions describing the hydrometeor
            default {}.
        useFuncArgDefaults : bool
            if parameters are not found in any of funcArgs, discreteProperties,
            parent's profile, then fall back to default values of the function.
            Helpful for debugging. default True.
        description : dict
            All properties of the hydrometeor. Most hydrometeors require at
            least 'sizeCenter', 'aspectRatio', 'mass', 'density',
            'crossSectionArea', and 'sizeDistribution'.

    """

    def __init__(
        self,
        parent,
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
        self.index = np.where(parent.profile.hydrometeor.values == name)[0][0]
        self._parentFull = parent
        self.coords = parent.coords

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

            thisProperty = func(**kw4Func)
        else:
            print('not callable', thisDesription)
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

        if self.calculationOrder is None:
            self.calculationOrder = DEFAULT_CALCULATION_ORDER()

        for key in self.description.keys():
            if key not in self.calculationOrder:
                self.calculationOrder.append(key)
            self.description.keys()

        self._keysToBeUsed = list(self.description.keys())

        for key in self.calculationOrder:

            value = self.description[key]
            print(key, value)

            thisProperty = self._arrayOrFunc(key, value)
            # to do: what if sizeCenter depends on other coords?
            if (key in ['sizeCenter', 'sizeBoundsWidth']):
                thisProperty = xr.DataArray(
                    thisProperty.data,
                    coords=[self.profile.sizeBin],
                    attrs={'unit': units.units[key]},
                )
            elif (key in ['sizeBounds']):
                thisProperty = xr.DataArray(
                    thisProperty.data,
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
        self.profile.merge(scatteringProperty, inplace=True)

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

        return super().__init__(*args, **kwargs)


class cloudDroplet(hydrometeor):
    """hydrometeor class with presets for cloud droplets.  """

    def __init__(
        self,
        *args,
        nBins=2,
        sizeBounds=size.linspaceBounds,
        sizeCenter=size.boundsToMid,
        sizeBoundsWidth=size.boundsWidth,
        sizeDistribution=sizeDistribution.monoDisperse,
        aspectRatio=1.0,
        mass=mass.ellipsoid,
        density=density.water,
        crossSectionArea=crossSectionArea.sphere,
        relativePermittivity=refractiveIndex.water.turner_kneifel_cadeddu,
        scattering=scattering.Rayleigh,
        fallVelocity=fallVelocity.khvorostyanov01_drops,
        Dmin=1e-5 - 1e-10,
        Dmax=1e-5 + 1e-10,
        Ntot=1e9,
        **kwargs,
    ):

        kwargs['nBins'] = nBins
        kwargs['sizeBounds'] = sizeBounds
        kwargs['sizeCenter'] = sizeCenter
        kwargs['sizeBoundsWidth'] = sizeBoundsWidth
        kwargs['sizeDistribution'] = sizeDistribution
        kwargs['aspectRatio'] = aspectRatio
        kwargs['mass'] = mass
        kwargs['density'] = density
        kwargs['crossSectionArea'] = crossSectionArea
        kwargs['relativePermittivity'] = relativePermittivity
        kwargs['scattering'] = scattering
        kwargs['fallVelocity'] = fallVelocity
        kwargs['Dmin'] = Dmin
        kwargs['Dmax'] = Dmax
        kwargs['Ntot'] = Ntot

        return super().__init__(*args, **kwargs)
