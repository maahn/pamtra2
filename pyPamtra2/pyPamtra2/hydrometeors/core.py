# -*- coding: utf-8 -*-
from collections import OrderedDict
from copy import deepcopy
import numpy as np
import xarray as xr
import toolz

from .density import softEllipsoid

class properties(object):
    """
    class to store hydrometeor properties. 
    """
    def __init__(
        self,
        parent,
        name = None, #or None, then str(index)
        kind = None, #liquid, ice
        nBins = None,
        maximumDimension = None,
        psd = None,
        aspectRatio = None,
        mass = None,
        density = None,
        crossSectionArea = None,
        discreteProperties = None
        ):

        self.name =name
        self.kind = kind
        self.nBins = nBins
        self.maximumDimension =maximumDimension
        self.psd =psd
        self.aspectRatio = aspectRatio
        self.mass = mass
        self.density = density
        self.crossSectionArea =crossSectionArea
        
        self.index = np.where(parent.profile.hydrometeor.values==name)[0][0]
        self._parent = parent

        if discreteProperties is None:
            self.discretePropertiesCoords = OrderedDict(self._parent.profile.coords)
            self.discretePropertiesCoords.pop('hydrometeor')
            self.discretePropertiesCoords.update(
                OrderedDict(
                    sizeBin =range(self.nBins)
                )
            )
            self.discreteProperties = xr.Dataset(coords=self.discretePropertiesCoords)
        else:
            self.discreteProperties = discreteProperties
            self.discretePropertiesCoords = OrderedDict(discreteProperties.coords)

        return
    
    def sel(self,**indexers):
        """
        Select subset by name
        see xarray.Dataset.sel() for documentation
        """
        
        for k in list(indexers.keys()):
            if not hasattr(indexers[k], '__iter__'):
                indexers[k] = [indexers[k]]
        
        subSet = hydrometeor(
        self._parent,
        name = self.name, #or None, then str(index)
        kind = self.kind, #liquid, ice
        nBins = self.nBins,
        maximumDimension = self.maximumDimension,
        psd = self.psd,
        aspectRatio = self.aspectRatio,
        mass = self.mass,
        density = self.density,
        crossSectionArea = self.crossSectionArea,
        discreteProperties=self.discreteProperties.sel(**indexers),
        )  
        
        return subSet

    def isel(self,**indexers):
        """
        Select subset by index
        see xarray.Dataset.sel() for documentation
        """  
        
        for k in list(indexers.keys()):
            if not hasattr(indexers[k], '__iter__'):
                indexers[k] = [indexers[k]]
        
        subSet = hydrometeor(
        self._parent,
        name = self.name, #or None, then str(index)
        kind = self.kind, #liquid, ice
        nBins = self.nBins,
        maximumDimension = self.maximumDimension,
        psd = self.psd,
        aspectRatio = self.aspectRatio,
        mass = self.mass,
        density = self.density,
        crossSectionArea = self.crossSectionArea,
        discreteProperties=self.discreteProperties.isel(**indexers),
        )  
        
        return subSet    
    
    def _selectHydro(self,arr):
        '''
        helper function to get the wanted hydrometeo out of a xarray.
        '''
        try:
            name = arr.name
        except AttributeError:
            # print('is not DataArray')
            pass
        else:
            try: 
                arr = deepcopy(self._parent.profile[name])
            except KeyError:
                # print('key %s not found in DataArray'%name)
                pass
            else:
                try:
                    arr = arr.sel(hydrometeor=self.name,drop=True)
                    # print('managed to extract %s'%self.name)
                except ValueError:
                    # print('cannot extract hydrometeor from DataArray %s'%name)
                    pass
        return arr
    
    def _arrayOrFunc(self,thisDesription,*args):
        
        thisShape = tuple()
        for val in self.discretePropertiesCoords.values():
            try:
                thisShape = thisShape + (len(val),)
            except TypeError:
                pass
        
        thisProperty = np.zeros(thisShape) * np.nan
        if hasattr(thisDesription, '__iter__') and callable(thisDesription[0]):
            print('callable')
            try:
                func, kwargs = thisDesription
            except ValueError:
                func = thisDesription[0]
                kwargs = {}
            # print(func, kwargs)
            #by transposing we comply with the broadcasting rules and save a loop
            # use .T instead of np.transpose, becasue of xr
            args = list(map(lambda arr: np.asarray(self._selectHydro(arr)).T, args))
            # print(args)
            kwargs = toolz.valmap(lambda arr: np.asarray(self._selectHydro(arr)).T, kwargs)
            # print(kwargs)
            thisProperty[:] = func(*args, **kwargs).T
        else:
            print('not callable',thisDesription)
            thisProperty[:] = thisDesription

        return thisProperty
    
    def _calculateHydroMaxDim(self):
        maximumDimension = self._arrayOrFunc(self.maximumDimension,self.nBins)
        self.discreteProperties['maximumDimension'] = xr.DataArray(
            maximumDimension,
            coords = self.discretePropertiesCoords.values(),
            dims = self.discretePropertiesCoords.keys(),
            attrs = {'unit' : 'm'}
        )
        
        if np.any(np.isnan(self.discreteProperties['maximumDimension'].values)):
            raise ValueError('found NAN in discreteProperties.maximumDimension')
        
        return self.discreteProperties['maximumDimension']
    
    def _calculateHydroPsd(self):
        
        if 'maximumDimension' not in self.discreteProperties.variables:
            raise RuntimeError('Estimate maximum dimension first!')
        
        particleSizeDistribution = self._arrayOrFunc(self.psd,self.discreteProperties['maximumDimension'].values)
        self.discreteProperties['particleSizeDistribution'] = xr.DataArray(
            particleSizeDistribution,
            coords = self.discretePropertiesCoords.values(),
            dims = self.discretePropertiesCoords.keys(),
            attrs = {'unit' : 'm^(-4)'}
        )
        
        if np.any(np.isnan(self.discreteProperties['maximumDimension'].values)):
            raise ValueError('found NAN in discreteProperties.maximumDimension')
                
        return self.discreteProperties['particleSizeDistribution']    
    
    def _calculateHydroAR(self):
        
        if 'maximumDimension' not in self.discreteProperties.variables:
            raise RuntimeError('Estimate maximum dimension first!')
        
        aspectRatio = self._arrayOrFunc(self.aspectRatio,self.discreteProperties['maximumDimension'].values)
        self.discreteProperties['aspectRatio'] = xr.DataArray(
            aspectRatio,
            coords = self.discretePropertiesCoords.values(),
            dims = self.discretePropertiesCoords.keys(),
            attrs = {'unit' : 'm^(-4)'}
        )
        
        if np.any(np.isnan(self.discreteProperties['maximumDimension'].values)):
            raise ValueError('found NAN in discreteProperties.maximumDimension')
                
        return self.discreteProperties['aspectRatio']       
    

    def _calculateHydroArea(self):
        
        if 'maximumDimension' not in self.discreteProperties.variables:
            raise RuntimeError('Estimate maximum dimension first!')
        
        crossSectionArea = self._arrayOrFunc(self.crossSectionArea,self.discreteProperties['maximumDimension'].values)
        self.discreteProperties['crossSectionArea'] = xr.DataArray(
            crossSectionArea,
            coords = self.discretePropertiesCoords.values(),
            dims = self.discretePropertiesCoords.keys(),
            attrs = {'unit' : 'm^2'}
        )
        
        if np.any(np.isnan(self.discreteProperties['maximumDimension'].values)):
            raise ValueError('found NAN in discreteProperties.maximumDimension')
                
        return self.discreteProperties['crossSectionArea']        

    def _calculateHydroMassDensity(self):
        if 'maximumDimension' not in self.discreteProperties.variables:
            raise RuntimeError('Estimate maximum dimension first!')
        
        if self.density[0] is softEllipsoid:
            print ('softsphere!')
            mass = self._arrayOrFunc(
                self.mass,
                self.discreteProperties['maximumDimension'].values,
                )
            density = self._arrayOrFunc(
                self.density,
                self.discreteProperties['maximumDimension'].values,
                self.discreteProperties['aspectRatio'].values,
                mass,
                )
        else:
            print ('NOT softsphere!')
            density = self._arrayOrFunc(
                self.density,
                self.discreteProperties['maximumDimension'].values,
                )
            mass = self._arrayOrFunc(
                self.mass,
                self.discreteProperties['maximumDimension'].values,
                self.discreteProperties['aspectRatio'].values,
                density,
                )

        self.discreteProperties['mass'] = xr.DataArray(
            mass,
            coords = self.discretePropertiesCoords.values(),
            dims = self.discretePropertiesCoords.keys(),
            attrs = {'unit' : 'kg'}
        )
        self.discreteProperties['density'] = xr.DataArray(
            density,
            coords = self.discretePropertiesCoords.values(),
            dims = self.discretePropertiesCoords.keys(),
            attrs = {'unit' : 'kg/m^3'}
        )
        if np.any(np.isnan(self.discreteProperties['maximumDimension'].values)):
            raise ValueError('found NAN in discreteProperties.maximumDimension')
                
        return self.discreteProperties['mass']       

    
    def calculateProperties(self):

        self._calculateHydroMaxDim()
        self._calculateHydroPsd()
        if self.aspectRatio is not None:
            self._calculateHydroAR()
        if self.crossSectionArea is not None:
            self._calculateHydroArea()
        if (self.density is not None) and (self.mass is not None):
            self._calculateHydroMassDensity()

        return self.discreteProperties
    
