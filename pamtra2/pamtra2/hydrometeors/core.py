# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import inspect

from .. import units
from .. import helpers


def DEFAULT_CALCULATION_ORDER():
    return [
                'sizeCenter',
                'aspectRatio',
                'density',
                'mass',
                'crossSectionArea',
                'sizeDistribution',
                'refractiveIndex',
                'scattering',
                ]


class hydrometeor(object):
    """generic class to store hydrometeor properties.

        Parameters
        ----------
        parent : pamtra2 object
            content of parent's class
        name : str, optional
            name of the hydrometeor
        kind :  {'liquid', 'ice'}, optional
            liquid or frozen hydrometeor?
        nBins : int, optional
            number of size bins
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
        nBins=None,
        discreteProperties=None,
        calculationOrder=None,
        useFuncArgDefaults=True,
        **kwargs
    ):

        self._nBins = nBins
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
                    coords=dict(sizeBin=range(self._nBins))
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
                elif k in fixedKwargs.keys():
                    kw4Func[k] = fixedKwargs[k]
                elif self.useFuncArgDefaults and (k in funcDefaults.keys()):
                    kw4Func[k] = funcDefaults[k]
                else:
                    raise KeyError('Did not find %s in provided kwargs or '
                                   'discreteProperties or profile or '
                                   'functions\'s defaultArgs' % k)

            thisProperty = func(**kw4Func)
        else:
            print('not callable', thisDesription)
            thisProperty = thisDesription

        return thisProperty

    def calculateProperties(self):
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

        for key in self.calculationOrder:

            value = self.description[key]
            print(key, value)

            thisProperty = self._arrayOrFunc(key, value, nBins=self._nBins)
            if (key == 'sizeCenter') and 'coords' not in dir(value):
                thisProperty = xr.DataArray(
                    thisProperty,
                    coords=[self.profile.sizeBin]
                    )

            self.profile[key] = thisProperty

            self.profile[key].attrs.update(
                {'unit': units.units[key]}
                )
        self._postProcessing()

        return self.profile

    def _postProcessing(self):
        """
        clean up 'artificial' dimensions required becaus xr.apply_ufunc cannot
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
