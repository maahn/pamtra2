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
            print(key,value)

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


# class properties(object):
#     """
#     class to store hydrometeor properties. 
#     """
#     def __init__(
#         self,
#         parent,
#         name = None, #or None, then str(index)
#         kind = None, #liquid, ice
#         nBins = None,
#         sizeCenter = None,
#         sizeDistribution = None,
#         aspectRatio = None,
#         mass = None,
#         density = None,
#         crossSectionArea = None,
#         discreteProperties = None
#         ):

#         self.name =name
#         self.kind = kind
#         self.nBins = nBins
#         self.sizeCenter =sizeCenter
#         self.sizeDistribution =sizeDistribution
#         self.aspectRatio = aspectRatio
#         self.mass = mass
#         self.density = density
#         self.crossSectionArea =crossSectionArea
        
#         self.frequencies = []

#         self.index = np.where(parent.profile.hydrometeor.values==name)[0][0]
#         self._parent = parent

#         if discreteProperties is None:
#             self.discretePropertiesCoords = OrderedDict(self._parent.profile.coords)
#             self.discretePropertiesCoords.pop('hydrometeor')
#             self.discretePropertiesCoords.update(
#                 OrderedDict(
#                     sizeBin =range(self.nBins)
#                 )
#             )
#             self.discretePropertiesCoordsExt = deepcopy(self.discretePropertiesCoords)
#             self.discretePropertiesCoordsExt.update(
#                 OrderedDict(
#                     frequency =self.frequencies
#                 )
#             )
#             self.discreteProperties = xr.Dataset(coords=self.discretePropertiesCoordsExt)
#         else:
#             self.discreteProperties = discreteProperties
#             self.discretePropertiesCoordsExt = OrderedDict(discreteProperties.coords)
#             self.discretePropertiesCoords = deepcopy(self.discretePropertiesCoordsExt)
#             self.discretePropertiesCoords.pop('frequency')

#         return
        
#     @property
#     def nFrequencies():
#         return len(self.frequencies)

#     def sel(self,**indexers):
#         """
#         Select subset by name
#         see xarray.Dataset.sel() for documentation
#         """
        
#         for k in list(indexers.keys()):
#             if not hasattr(indexers[k], '__iter__'):
#                 indexers[k] = [indexers[k]]
        
#         subSet = hydrometeor(
#         self._parent,
#         name = self.name, #or None, then str(index)
#         kind = self.kind, #liquid, ice
#         nBins = self.nBins,
#         sizeCenter = self.sizeCenter,
#         psd = self.sizeDistribution,
#         aspectRatio = self.aspectRatio,
#         mass = self.mass,
#         density = self.density,
#         crossSectionArea = self.crossSectionArea,
#         discreteProperties=self.discreteProperties.sel(**indexers),
#         )  
        
#         return subSet

#     def isel(self,**indexers):
#         """
#         Select subset by index
#         see xarray.Dataset.sel() for documentation
#         """  
        
#         for k in list(indexers.keys()):
#             if not hasattr(indexers[k], '__iter__'):
#                 indexers[k] = [indexers[k]]
        
#         subSet = hydrometeor(
#         self._parent,
#         name = self.name, #or None, then str(index)
#         kind = self.kind, #liquid, ice
#         nBins = self.nBins,
#         sizeCenter = self.sizeCenter,
#         psd = self.sizeDistribution,
#         aspectRatio = self.aspectRatio,
#         mass = self.mass,
#         density = self.density,
#         crossSectionArea = self.crossSectionArea,
#         discreteProperties=self.discreteProperties.isel(**indexers),
#         )  
        
#         return subSet    
    
#     def _selectHydro(self,arr):
#         '''
#         helper function to get the wanted hydrometeo out of a xarray.
#         '''
#         try:
#             name = arr.name
#         except AttributeError:
#             # print('is not DataArray')
#             pass
#         else:
#             try: 
#                 arr = deepcopy(self._parent.profile[name])
#             except KeyError:
#                 # print('key %s not found in DataArray'%name)
#                 pass
#             else:
#                 try:
#                     arr = arr.sel(hydrometeor=self.name,drop=True)
#                     # print('managed to extract %s'%self.name)
#                 except ValueError:
#                     # print('cannot extract hydrometeor from DataArray %s'%name)
#                     pass
#         return arr
    
#     def _arrayOrFunc(self,thisDesription,**fixedKwargs):
        
#         thisShape = tuple()
#         for val in self.discretePropertiesCoords.values():
#             try:
#                 thisShape = thisShape + (len(val),)
#             except TypeError:
#                 pass
        
#         thisProperty = np.zeros(thisShape) * np.nan
#         if hasattr(thisDesription, '__iter__') and callable(thisDesription[0]):
#             print('callable')
#             try:
#                 func, kwargs = thisDesription
#             except ValueError:
#                 func = thisDesription[0]
#                 kwargs = {}
#             # print(func, kwargs)

#             #inspect function to get the required arguments
#             funcArgs, funcVarargs, funcKeywords, funcDefaults = inspect.getargspec(func)
            
#             if funcDefaults is None:
#                 funcDefaults = {}
#             else:
#                 funcDefaults = dict(zip(funcArgs[-len(funcDefaults):],funcDefaults))

#             #where do we find the required data?
#             kw4Func = {}
#             for k in funcArgs:
#                 if k in kwargs.keys():
#                     kw4Func[k] = kwargs[k]
#                 elif k in self.discreteProperties.keys():
#                     kw4Func[k] = self.discreteProperties[k]
#                 elif k in self._parent.profile.keys():
#                     kw4Func[k] = self._parent.profile[k]
#                 elif k in fixedKwargs.keys():
#                     kw4Func[k] = fixedKwargs[k]
#                 elif k in funcDefaults.keys():
#                     kw4Func[k] = funcDefaults[k]
#                 else:
#                     raise KeyError('Did not find %s in provided kwargs or discreteProperties'
#                         ' or profile or functions\'s defaultArgs'%k )

#             #by transposing we comply with the broadcasting rules and save a loop
#             # use .T instead of np.transpose, becasue of xr
#             # import pdb;pdb.set_trace()
#             kw4Func = toolz.valmap(lambda arr: np.asarray(self._selectHydro(arr)).T, kw4Func)
#             # print(kwargs)
#             thisProperty[:] = func(**kw4Func).T
#         else:
#             print('not callable',thisDesription)
#             thisProperty[:] = thisDesription

#         return thisProperty
    
#     def _calculateHydroMaxDim(self):
#         sizeCenter = self._arrayOrFunc(self.sizeCenter,nBins=self.nBins)
#         self.discreteProperties['sizeCenter'] = xr.DataArray(
#             sizeCenter,
#             coords = self.discretePropertiesCoords.values(),
#             dims = self.discretePropertiesCoords.keys(),
#             attrs = {'unit' : 'm'}
#         )
        
#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
        
#         return self.discreteProperties['sizeCenter']
    
#     def _calculateHydroPsd(self):
        
#         if 'sizeCenter' not in self.discreteProperties.variables:
#             raise RuntimeError('Estimate maximum dimension first!')
        
#         particleSizeDistribution = self._arrayOrFunc(self.sizeDistribution)
#         self.discreteProperties['particleSizeDistribution'] = xr.DataArray(
#             particleSizeDistribution,
#             coords = self.discretePropertiesCoords.values(),
#             dims = self.discretePropertiesCoords.keys(),
#             attrs = {'unit' : 'm^(-4)'}
#         )
        
#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
                
#         return self.discreteProperties['particleSizeDistribution']    
    
#     def _calculateHydroAR(self):
        
#         if 'sizeCenter' not in self.discreteProperties.variables:
#             raise RuntimeError('Estimate maximum dimension first!')
        
#         aspectRatio = self._arrayOrFunc(self.aspectRatio)
#         self.discreteProperties['aspectRatio'] = xr.DataArray(
#             aspectRatio,
#             coords = self.discretePropertiesCoords.values(),
#             dims = self.discretePropertiesCoords.keys(),
#             attrs = {'unit' : 'm^(-4)'}
#         )
        
#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
                
#         return self.discreteProperties['aspectRatio']       
    

#     def _calculateHydroArea(self):
        
#         if 'sizeCenter' not in self.discreteProperties.variables:
#             raise RuntimeError('Estimate maximum dimension first!')
        
#         crossSectionArea = self._arrayOrFunc(self.crossSectionArea)
#         self.discreteProperties['crossSectionArea'] = xr.DataArray(
#             crossSectionArea,
#             coords = self.discretePropertiesCoords.values(),
#             dims = self.discretePropertiesCoords.keys(),
#             attrs = {'unit' : 'm^2'}
#         )
        
#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
                
#         return self.discreteProperties['crossSectionArea']        

#     def _calculateHydroMassDensity(self):
#         if 'sizeCenter' not in self.discreteProperties.variables:
#             raise RuntimeError('Estimate maximum dimension first!')
        
#         #for softEllipsoids, we have to estimate mass FIRST before we can look at density
#         if self.density[0] is softEllipsoid:
#             print ('softsphere!')
#             mass = self._arrayOrFunc(
#                 self.mass,
#                 )
#             density = self._arrayOrFunc(
#                 self.density,
#                 mass =mass,
#                 )
#         #for solid ellispoids, it is density first. 
#         else:
#             print ('NOT softsphere!')
#             density = self._arrayOrFunc(
#                 self.density,
#                 )
#             mass = self._arrayOrFunc(
#                 self.mass,
#                 density = density,
#                 )

#         self.discreteProperties['mass'] = xr.DataArray(
#             mass,
#             coords = self.discretePropertiesCoords.values(),
#             dims = self.discretePropertiesCoords.keys(),
#             attrs = {'unit' : 'kg'}
#         )
#         self.discreteProperties['density'] = xr.DataArray(
#             density,
#             coords = self.discretePropertiesCoords.values(),
#             dims = self.discretePropertiesCoords.keys(),
#             attrs = {'unit' : 'kg/m^3'}
#         )
#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
                
#         return self.discreteProperties['mass']       

#     def _calculateRefractiveIndex(self):
        
#         if 'sizeCenter' not in self.discreteProperties.variables:
#             raise RuntimeError('Estimate maximum dimension first!')
        
#         self.discreteProperties['refractiveIndex'] = xr.DataArray(
#             np.zeros(self.discreteProperties['sizeCenter'].shape + (self.nFrequencies,),dtype=np.complex),
#             coords = self.discretePropertiesCoordsExt.values(),
#             dims = self.discretePropertiesCoordsExt.keys(),
#             attrs = {'unit' : '?'}
#         )
        
#         #some random number for now
#         self.discreteProperties['refractiveIndex'][:] = 1.77701000 +2.18657774e-05j

#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
                
#         return self.discreteProperties['refractiveIndex']        

#     def _calculateSingleScattering(self):
        
#         if 'sizeCenter' not in self.discreteProperties.variables:
#             raise RuntimeError('Estimate maximum dimension first!')
        
#         self.discreteProperties['backScatteringCrossSection'] = xr.DataArray(
#             np.zeros(self.discreteProperties['sizeCenter'].shape + (self.nFrequencies,)),
#             coords = self.discretePropertiesCoordsExt.values(),
#             dims = self.discretePropertiesCoordsExt.keys(),
#             attrs = {'unit' : 'm^2'}
#         )
        
#         #some random number for now
#         self.discreteProperties['backScatteringCrossSection'][:] = 1.77701000 +2.18657774e-05j

#         if np.any(np.isnan(self.discreteProperties['sizeCenter'].values)):
#             raise ValueError('found NAN in discreteProperties.sizeCenter')
                
#         return self.discreteProperties['backScatteringCrossSection']        


#     def addScatteringProperties(
#         refractiveIndex = None,
#         singleScattering = None,
#         ):

#         self.refractiveIndex = refractiveIndex
#         self.singleScattering = singleScattering

#         return

#     def calculateScatteringProperties(self):

#         if self.refractiveIndex is not None:
#             self._calculateRefractiveIndex()
#         if self.singleScattering is not None:
#             self._calculateSingleScattering()

#         return self.discreteProperties
    

#     def calculateProperties(self):

#         self._calculateHydroMaxDim()
#         self._calculateHydroPsd()
#         if self.aspectRatio is not None:
#             self._calculateHydroAR()
#         if self.crossSectionArea is not None:
#             self._calculateHydroArea()
#         if (self.density is not None) and (self.mass is not None):
#             self._calculateHydroMassDensity()

#         return self.discreteProperties
    
