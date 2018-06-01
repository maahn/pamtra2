# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import inspect

from .. import units


class hydrometeor(object):
    """
    generic class to store hydrometeor properties.
    """
    def __init__(
        self,
        parent,
        name=None,  # or None, then str(index)
        kind=None,  # liquid, ice
        nBins=None,
        discreteProperties=None,
        calculationOrder=None,
        funcArgs={},
        useFuncArgDefaults=True,
        **kwargs
    ):

        self.name = name
        self.kind = kind
        self.nBins = nBins
        self.calculationOrder = calculationOrder
        self.funcArgs = funcArgs
        self.useFuncArgDefaults = useFuncArgDefaults
        self.description = kwargs

        self.index = np.where(parent.profile.hydrometeor.values == name)[0][0]
        self._parentFull = parent

        if discreteProperties is None:
            discreteProperties = xr.Dataset(
                    coords=dict(sizeBin=range(self.nBins))
                    )
        self.discreteProperties = discreteProperties

        return

    @property
    def _parentProfile(self):
        '''
        Limited version of parent.profile. Contains only data belonging to
        the hydrometeor.
        '''
        return self._parentFull.profile.sel(hydrometeor=self.name, drop=True)

    def _arrayOrFunc(self, thisDesription, **fixedKwargs):

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
                if k in self.funcArgs.keys():
                    kw4Func[k] = self.funcArgs[k]
                elif k in self.discreteProperties.keys():
                    kw4Func[k] = self.discreteProperties[k]
                elif k in self._parentProfile.keys():
                    kw4Func[k] = self._parentProfile[k]
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

        if self.calculationOrder is None:
            self.calculationOrder = self.description.keys()

        for key in self.calculationOrder:
            value = self.description[key]
            print(key, value)

            thisProperty = self._arrayOrFunc(value, nBins=self.nBins)
            if (key == 'sizeCenter') and 'coords' not in dir(value):
                thisProperty = xr.DataArray(
                    thisProperty,
                    coords=[self.discreteProperties.sizeBin]
                    )

            self.discreteProperties[key] = thisProperty

            self.discreteProperties[key].attrs.update(
                {'unit': units.units[key]}
                )

        return self.discreteProperties


class softEllipsoidFixedDensity(hydrometeor):

    def __init__(self, *args, **kwargs):
        if 'calculationOrder' not in kwargs.keys():
            kwargs['calculationOrder'] = [
                'sizeCenter',
                'aspectRatio',
                'density',
                'mass',
                'crossSectionArea',
                'sizeDistribution'
                ]
        return super().__init__(*args, **kwargs)


class softEllipsoidMassSize(hydrometeor):

    def __init__(self, *args, **kwargs):
        if 'calculationOrder' not in kwargs.keys():
            kwargs['calculationOrder'] = [
                'sizeCenter',
                'aspectRatio',
                'mass',
                'density',
                'crossSectionArea',
                'sizeDistribution'
                ]
        return super().__init__(*args, **kwargs)
