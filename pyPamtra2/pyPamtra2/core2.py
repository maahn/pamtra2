# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function
from builtins import super 

from collections import OrderedDict
from copy import deepcopy
import numpy as np
import xarray as xr
import json
from functools import lru_cache

import pyPamtraRadarSimulator
import pyPamtraRadarMoments
import refractive

from . import configuration
from . import helpers
from . import decorators
from . import constants
from . import sizeDistributions

__version__ = '0.1'


class profile(object):

  def  __init__(
    self,
    nLayer = 1,
    settings = None,
    additionalDims = {},
    dataset = xr.Dataset()
    ):


    self.additionalDims =OrderedDict(additionalDims)
    self.nLayer = nLayer

    self.hydrometeors = []
    self.nHydroBins = 0
    self.radarFrequencies=[]
    self.radarPolarisations = []
    self.nFFT = 0
    self.nPeaks = 0


    if settings is None:
      self.settings = configuration.Settings(self.radarFrequencies,self.hydrometeors)
    else:
      self.settings = configuration.Settings(settings)


    self._createXrObjects()
    self._fillDataObject()


    if not isinstance(self._dataset,xr.Dataset):
      self._dataset = xr.Dataset(self._dataset)
    self._dataset.update(dataset)


  def __repr__(self):
    repr = {
    'nLayer': self.nLayer,
    'settings': self.settings,
    'additionalDims': self.additionalDims,
    'dataset': self._dataset,
    }
    return str(repr)

  @property
  def height(self):
    """ height vector """
    return self._dataset['height']
  @height.setter
  def height(self,data):
    self._dataset['height'].values[:] = data
  @height.deleter
  def height(self):
    self._dataset['height'].values[:] = np.nan

  @property
  def temperature(self):
    """ temperature vector """
    return self._dataset['temperature']
  @temperature.setter
  def temperature(self,data):
    self._dataset['temperature'].values[:] = data
  @temperature.deleter
  def temperature(self):
    self._dataset['temperature'].values[:] = np.nan

  @property
  def pressure(self):
    """ pressure vector """
    return self._dataset['pressure']
  @pressure.setter
  def pressure(self,data):
    self._dataset['pressure'].values[:] = data
  @pressure.deleter
  def pressure(self):
    self._dataset['pressure'].values[:] = np.nan

  @property
  def relativeHumidity(self):
    """ relativeHumidity vector """
    return self._dataset['relativeHumidity']
  @relativeHumidity.setter
  def relativeHumidity(self,data):
    self._dataset['relativeHumidity'].values[:] = data
  @relativeHumidity.deleter
  def relativeHumidity(self):
    self._dataset['relativeHumidity'].values[:] = np.nan

  @property
  def horizontalWind(self):
    """ horizontalWind vector """
    return self._dataset['horizontalWind']
  @horizontalWind.setter
  def horizontalWind(self,data):
    self._dataset['horizontalWind'].values[:] = data
  @horizontalWind.deleter
  def horizontalWind(self):
    self._dataset['horizontalWind'].values[:] = np.nan

  @property
  def verticalWind(self):
    """ verticalWind vector """
    return self._dataset['verticalWind']
  @verticalWind.setter
  def verticalWind(self,data):
    self._dataset['verticalWind'].values[:] = data
  @verticalWind.deleter
  def verticalWind(self):
    self._dataset['verticalWind'].values[:] = np.nan

  @property
  def eddyDissipationRate(self):
    """ eddyDissipationRate vector """
    return self._dataset['eddyDissipationRate']
  @eddyDissipationRate.setter
  def eddyDissipationRate(self,data):
    self._dataset['eddyDissipationRate'].values[:] = data
  @eddyDissipationRate.deleter
  def eddyDissipationRate(self):
    self._dataset['eddyDissipationRate'].values[:] = np.nan

  @property
  def specBroadening(self):
    assert not np.all(np.isnan(self._dataset.eddyDissipationRate.values)), 'found only nan in eddyDissipationRate'
    assert not np.all(np.isnan(self._dataset.horizontalWind.values)), 'found only nan in horizontalWind'
    assert not np.any(np.isnan(self._dataset.height.values)), 'found nan in height'
    assert not len(self.radarFrequencies) == 0, 'define radarFrequencies first!'

    coordsND = self.additionalDimList + [self._dataset.layer, self._dataset.frequency]
    shapeND = tuple(map(len,coordsND))
    if 'specBroadening' not in self._dataset.variables:
      self._dataset['specBroadening'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'m/s'})

    for ff, freq in enumerate(self.radarFrequencies):

      specBroadening = pyPamtraRadarSimulator.calcSpectralBraodening(
        self._dataset.eddyDissipationRate.values[:],
        self._dataset.horizontalWind.values[:],
        self._dataset.height.values[:],
        self.settings['radarProperties'][repr(freq)]['fwhrBeamwidthDeg'],
        self.settings['radarProperties'][repr(freq)]['integrationTime'],
        freq,
      )
      self._dataset.specBroadening.values[...,ff] = specBroadening
    return self._dataset.specBroadening
  @specBroadening.setter
  def specBroadening(self,data):
    coordsND = self.additionalDimList + [self._dataset.layer, self._dataset.frequency]
    shapeND = tuple(map(len,coordsND))   
    if 'specBroadening' not in self._dataset.variables:
      self._dataset['specBroadening'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'m/s'})
    self._dataset['specBroadening'].values[:] = data
  @specBroadening.deleter
  def specBroadening(self):
    if 'specBroadening' not in self._dataset.variables:
      raise AttributeError('specBroadening')
    self._dataset['specBroadening'].values[:] = np.nan


  def getSpectralBroadening(self):

    assert not np.all(np.isnan(self._dataset.eddyDissipationRate.values)), 'found only nan in eddyDissipationRate'
    assert not np.all(np.isnan(self._dataset.horizontalWind.values)), 'found only nan in horizontalWind'
    assert not np.any(np.isnan(self._dataset.height.values)), 'found nan in height'
    assert not len(self.radarFrequencies) == 0, 'define radarFrequencies first!'

    coordsND = self.additionalDimList + [self._dataset.layer, self._dataset.frequency]
    shapeND = tuple(map(len,coordsND))
    if 'specBroadening' not in self._dataset.keys():
      self._dataset['specBroadening'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'m/s'})

    for ff, freq in enumerate(self.radarFrequencies):

      specBroadening = pyPamtraRadarSimulator.calcSpectralBraodening(
        self._dataset.eddyDissipationRate.values[:],
        self._dataset.horizontalWind.values[:],
        self._dataset.height.values[:],
        self.settings['radarProperties'][repr(freq)]['fwhrBeamwidthDeg'],
        self.settings['radarProperties'][repr(freq)]['integrationTime'],
        freq,
      )
      self._dataset.specBroadening.values[...,ff] = specBroadening

    @property
    def specBroadening(self):
      """ specBroadening vector """
      return self._dataset['specBroadening']

    return


  def _createXrObjects(self):

    coordsData = deepcopy(self.additionalDims)
    coordsData.update(OrderedDict(
      layer = range(self.nLayer),
      hydrometeor = self.hydrometeors,
      hydroBin = range(self.nHydroBins),
      frequency = self.radarFrequencies,
      activePolarisation = self.radarPolarisations,
      dopplerVelBin = range(self.nFFT),
      dopplerPeak = range(self.nPeaks),
      radarPolarisations = self.radarPolarisations,
      ))
    self._dataset = xr.Dataset(coords=coordsData)

    return

  def _fillDataObject(self):

    #for convinience, store additonal variables
    self.additionalDimList = []
    for cc in self.additionalDims.keys():
      self.additionalDimList.append((self._dataset[cc]))


    self.coordsGeometry = self.additionalDimList + [self._dataset.layer]
    self.coordsGeometryHydro = self.coordsGeometry + [self._dataset.hydrometeor]
    self.coordsGeometryHydroSpectrum = self.coordsGeometryHydro + [self._dataset.hydroBin]


    for var,unit,coords in [
      ('height','m', self.coordsGeometry),
      ('temperature','K', self.coordsGeometry),
      ('pressure','Pa', self.coordsGeometry),
      ('relativeHumidity','%', self.coordsGeometry),
      ('horizontalWind','m/s', self.coordsGeometry),
      ('verticalWind','m/s', self.coordsGeometry),
      ('eddyDissipationRate','m2/s3', self.coordsGeometry),
      # ('hydroSize','m',self.coordsGeometryHydroSpectrum),
      # # ('hydroSizeBinWidth','m',self.coordsGeometryHydroSpectrum),
      # ('hydroRho','kg/m3',self.coordsGeometryHydroSpectrum),
      # ('hydroMass','kg',self.coordsGeometryHydroSpectrum),
      # ('hydroCrossSectionArea','m2',self.coordsGeometryHydroSpectrum),
      # ('hydroPSD','1/m3',self.coordsGeometryHydroSpectrum),
      ]:
        thisShape = tuple(map(len,coords))
        self._dataset[var] = xr.DataArray(
          np.zeros(thisShape)*np.nan,
          coords=coords,
          attrs={'unit':unit},
          )

    self._dataset['hydro_name'] = xr.DataArray(np.arange(len(self._dataset.hydrometeor)).astype('S10'),coords=[self._dataset.hydrometeor])

    return

'''
add function to add radar properties
'''



class foo(object):
  def  __init__(
    self,
    nLayer = 1,
    hydrometeors = ['hydrometeor_1'],
    nHydroBins = 0,
    frequencies=[35.5],
    radarPolarisations = ['NN'],
    nFFT = 0,
    nPeaks = 1,
    settings = None,
    additionalDims = {}
    ):
    self.nLayer =nLayer 
    self.hydrometeors =hydrometeors
    self.nHydro =len(hydrometeors) 
    self.nHydroBins =nHydroBins 
    self.radarFrequencies=np.array(frequencies)
    self.radarPolarisations =radarPolarisations 
    self.nFFT =nFFT 
    self.nPeaks =nPeaks 
    self.additionalDims =OrderedDict(additionalDims)
    self.dopplerVelBins = range(nFFT)

    if settings is None:
      self.settings = configuration.Settings(self.radarFrequencies,self.hydrometeors)
    else:
      self.settings = configuration.Settings(settings)

    self._createXrObjects()
    self._fillDataObject()


    return





  def _fillDataObject(self):

    #for convinience, store additonal variables
    self.additionalDimList = []
    for cc in self.additionalDims.keys():
      self.additionalDimList.append((self._dataset[cc]))


    self.coordsGeometry = self.additionalDimList + [self._dataset.layer]
    self.coordsGeometryHydro = self.coordsGeometry + [self._dataset.hydrometeor]
    self.coordsGeometryHydroSpectrum = self.coordsGeometryHydro + [self._dataset.hydroBin]


    for var,unit,coords in [
      ('height','m', self.coordsGeometry),
      ('temperature','K', self.coordsGeometry),
      ('pressure','Pa', self.coordsGeometry),
      ('relativeHumidity','%', self.coordsGeometry),
      ('horizontalWind','m/s', self.coordsGeometry),
      ('verticalWind','m/s', self.coordsGeometry),
      ('eddyDissipationRate','m2/s3', self.coordsGeometry),
      ('hydroSize','m',self.coordsGeometryHydroSpectrum),
      # ('hydroSizeBinWidth','m',self.coordsGeometryHydroSpectrum),
      ('hydroRho','kg/m3',self.coordsGeometryHydroSpectrum),
      ('hydroMass','kg',self.coordsGeometryHydroSpectrum),
      ('hydroCrossSectionArea','m2',self.coordsGeometryHydroSpectrum),
      ('hydroPSD','1/m3',self.coordsGeometryHydroSpectrum),
      ]:
        thisShape = tuple(map(len,coords))
        self._dataset[var] = xr.DataArray(
          np.zeros(thisShape)*np.nan,
          coords=coords,
          attrs={'unit':unit},
          )

    self._dataset['hydro_name'] = xr.DataArray(np.arange(len(self._dataset.hydrometeor)).astype('S10'),coords=[self._dataset.hydrometeor])

    return


  def addHydrometeorBinDimension(self,nHydroBins):

    """
    In case you want to add it later due to memory concerns
    """

    self._dataset = self._dataset.reindex(hydroBin = range(nHydroBins))
    self.nHydroBins = nHydroBins
    return 

  # def addLiquidHydrometeor(self,hydroIndex):
  #   """
  #   we need something smarter here!
  #   """

  #   self._dataset['hydroRho'][...,hydroIndex,:] = constants.rhoWater
  #   self._dataset['hydroMass'][...,hydroIndex,:] = (self._dataset['hydroSize'][...,hydroIndex,:]/2.)**3 * np.pi * 4./3. * constants.rhoWater
  #   self._dataset['hydroCrossSectionArea'][...,hydroIndex,:] = (self._dataset['hydroSize'][...,hydroIndex,:]/2.)**2 * np.pi

  #   return

  def addHydrometeorProperties(
    self,
    hydroIndexName,
    diameters,
    psd,
    hydroType = 'liquid',
    useLogSizeSpacing = False,
    psdNormalized = False,
    psdFuncArgs = [],
    psdFuncKwArgs = {},
    hydrometeorProperties = {},
    ):
    """
    diameters : array_like
      if len == 2 from Dmin to Dmax
    mass and area default to properties of liquid water
    """

    assert self.nHydroBins >0, 'First use addHydrometeorBinDimension to create hydrometeor bin dimension.'
    assert hydroType in ['liquid','ice','snow']

    if isinstance(hydroIndexName,str):
      hydroName = hydroIndexName
      hydroIndex = self.hydrometeors.index(hydroIndexName)
    else:
      hydroIndex = hydroIndexName
      hydroName = self.hydrometeors[hydroIndexName]

    self.settings['hydrometeorProperties'][hydroName]['type'] = hydroType
    self.settings['hydrometeorProperties'][hydroName].update(
      configuration.DEFAULT_HYDROMETEOR_PROPERTIES_BY_TYPE[hydroType]
      )
    for key,value in hydrometeorProperties.items():
      self.settings['hydrometeorProperties'][hydroName][key] = value



    if len(diameters) ==2 and not useLogSizeSpacing:
      self._dataset.hydroSize.values[...,hydroIndex,:] = np.linspace(diameters[0],diameters[1],self.nHydroBins)
    elif len(diameters) ==2 and useLogSizeSpacing:
      self._dataset.hydroSize.values[...,hydroIndex,:] = np.logspace(np.log10(diameters[0]),np.log10(diameters[1]),self.nHydroBins)
    else:
      self._dataset.hydroSize.values[...,hydroIndex,:] = diameters
    # self._dataset.hydroSizeBinWidth.values[...,hydroIndex,:] = np.gradient(self._dataset.hydroSize.values[...,hydroIndex,:])


    #either explicit distribution or 
    if callable(psd):
      self._dataset.hydroPSD.values[...,hydroIndex,:] = psd(self,hydroIndex,*psdFuncArgs,**psdFuncKwArgs)
    else:
      self._dataset.hydroPSD.values[...,hydroIndex,:] = psd

    #roll back normalization
    if psdNormalized:
      self._dataset.hydroPSD.values[...,hydroIndex,:] = self._dataset.hydroPSD.values[...,hydroIndex,:]/np.gradient(self._dataset.hydroSize.values[...,hydroIndex,:],axis=-1)


    self._addHydrometeorDensityMassArea(hydroIndex)

    return

  def _addHydrometeorDensityMassArea(self,hydroIndex):

    hydroName = self.hydrometeors[hydroIndex]
    hydroType = self.settings['hydrometeorProperties'][hydroName]['type']

    if hydroType == 'liquid':
      self._dataset.hydroRho.values[...,hydroIndex,:] = constants.rhoWater
      self._dataset.hydroMass.values[...,hydroIndex,:] = (self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**3 * np.pi * 4./3. * constants.rhoWater
      self._dataset.hydroCrossSectionArea.values[...,hydroIndex,:] = (self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**2 * np.pi
    elif hydroType == 'ice':
      self._dataset.hydroRho.values[...,hydroIndex,:] = constants.rhoIce
      self._dataset.hydroMass.values[...,hydroIndex,:] = (self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**3 * np.pi * 4./3. * constants.rhoIce
      self._dataset.hydroCrossSectionArea.values[...,hydroIndex,:] = (self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**2 * np.pi
    elif hydroType == 'snow':
      self._dataset.hydroMass.values[...,hydroIndex,:] = (self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**self.settings['hydrometeorProperties'][hydroName]['massSizeB_snow'] * self.settings['hydrometeorProperties'][hydroName]['massSizeA_snow']
      self._dataset.hydroRho.values[...,hydroIndex,:] = self._dataset.hydroMass.values[...,hydroIndex,:] / ((self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**3 * np.pi * 4./3.)
      self._dataset.hydroCrossSectionArea.values[...,hydroIndex,:] = (self._dataset.hydroSize.values[...,hydroIndex,:]/2.)**self.settings['hydrometeorProperties'][hydroName]['areaSizeB_snow'] * self.settings['hydrometeorProperties'][hydroName]['areaSizeA_snow']
      self._dataset.hydroRho[...,hydroIndex,:].values[
        self._dataset.hydroRho.values[...,hydroIndex,:]>self.settings['hydrometeorProperties'][hydroName]['maxDensity_snow']
        ] = self.settings['hydrometeorProperties'][hydroName]['maxDensity_snow']
    else:
      ValueError("hydroType %s not in ['liquid','ice','snow']"%(hydroType))

    return

  def getSpectralBroadening(self,method='pyPamtraRadarSimulator'):

    assert not np.all(np.isnan(self._dataset.eddyDissipationRate.values)), 'found nan in eddyDissipationRate'
    assert not np.all(np.isnan(self._dataset.horizontalWind.values)), 'found nan in horizontalWind'
    assert not np.any(np.isnan(self._dataset.height.values)), 'found nan in height'
    assert method == 'pyPamtraRadarSimulator'

    coordsND = self.additionalDimList + [self._tmp.layer, self._tmp.frequency]
    shapeND = tuple(map(len,coordsND))
    if 'specBroadening' not in self._tmp.keys():
      self._tmp['specBroadening'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'m/s'})

    for ff, freq in enumerate(self.radarFrequencies):

      specBroadening = pyPamtraRadarSimulator.calcSpectralBraodening(
        self._dataset.eddyDissipationRate.values[:],
        self._dataset.horizontalWind.values[:],
        self._dataset.height.values[:],
        self.settings['radarProperties'][repr(freq)]['fwhrBeamwidthDeg'],
        self.settings['radarProperties'][repr(freq)]['integrationTime'],
        freq,
      )
      self._tmp.specBroadening.values[...,ff] = specBroadening
    return


  def getRefractiveIndex(self,):

    """
    Relies on the settings refractiveIndexModel and refractiveIndexMix_snow in hydrometeorProperties
    """

    coordsND = self.additionalDimList + [self._tmp.layer, self._tmp.frequency,self._tmp.hydrometeor, self._tmp.hydroBin]
    shapeND = tuple(map(len,coordsND))
    if 'K2' not in self._tmp.keys():
      self._tmp['K2'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'-'})
    if 'eps' not in self._tmp.keys():
      self._tmp['eps'] = xr.DataArray(np.zeros(shapeND).astype(np.complex)*np.nan,coords=coordsND,attrs={'unit':'-'})

    Temperatures  = self._dataset['temperature'].values
    Frequencies = self.radarFrequencies * 1e9

    for ff, freq in enumerate(Frequencies):
      for hh, hydroName in enumerate(self.hydrometeors):

        if self.settings['hydrometeorProperties'][hydroName]['type'] == 'liquid':
          eps = refractive.water.eps(
            Temperatures,
            freq,
            model=self.settings['hydrometeorProperties'][hydroName]['refractiveIndexModel'],
            )[...,np.newaxis]
        elif self.settings['hydrometeorProperties'][hydroName]['type'] == 'ice':
          eps = refractive.ice.eps(
            Temperatures,
            freq,
            model=self.settings['hydrometeorProperties'][hydroName]['refractiveIndexModel'],
            )[...,np.newaxis]
        elif self.settings['hydrometeorProperties'][hydroName]['type'] == 'snow':
          Densities = self._dataset['hydroRho'].sel(hydrometeor=hydroName).values
          eps = decorators.NDto2DtoND(
            referenceIn=2,
            noOfInDimsToKeep=1, 
            convertInputs=[0,2],
            convertOutputs=[0],
            verbosity=self.settings['general']['verbosity'],
            )(
            refractive.snow.eps
            )(
            Temperatures[...,np.newaxis], # new dimension becasue densities has hydro bin dims 
            freq,
            Densities,
            model_ice=self.settings['hydrometeorProperties'][hydroName]['refractiveIndexModel'],
            model_mix=self.settings['hydrometeorProperties'][hydroName]['refractiveIndexMix_snow'],
           )
        else:
          raise KeyError("self.settings['hydrometeorProperties'][%s]['type'] must be in ['liquid','ice','snow'] but is %s"%(hydroName,self.settings['hydrometeorProperties'][hydroName]['type']))

        self._tmp['eps'].isel(hydrometeor=hh,frequency=ff).values[:] = eps
    self._tmp['K2'].values[:] = refractive.utilities.K2(self._tmp['eps'].values)

    return

  def getScattering(self,method='Rayleigh'):
    """
    Add more sophistacted methods here
    """

    assert 'K2' in self._tmp.keys(), 'K2 not found, run getRefractiveIndex first'
    assert method == 'Rayleigh'


    coordsND = self.additionalDimList + [self._tmp.layer, self._tmp.frequency,self._tmp.activePolarisation,self._tmp.hydrometeor, self._tmp.hydroBin]
    shapeND = tuple(map(len,coordsND))
    if 'hydro_backSpec' not in self._tmp.keys():
      self._tmp['hydro_backSpec'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'?'})


    for ff, freq in enumerate(self.radarFrequencies):
      for pp, pol in enumerate(self.radarPolarisations):
        for hh in range(self.nHydro):
          if method == 'Rayleigh':
            assert pol == 'NN'

            back_spec = helpers.rayleigh(
              self._dataset['hydroSize'].values[...,hh,:], 
              self._tmp['K2'].values[...,ff,pp,hh], 
              freq
              )
            self._tmp['hydro_backSpec'].values[...,ff,pp,hh,:] = back_spec * self._dataset['hydroPSD'][...,hh,:].values

    return


  def getPIA(self):
    """
    Add more sophistacted methods here
    """

    coordsND = self.additionalDimList + [self._tmp.layer,self._tmp.frequency]
    shapeND = tuple(map(len,coordsND))
    if 'pathIntegratedAttenuation' not in self._tmp.keys():
      self.results['pathIntegratedAttenuation'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'?'})


      self.results['pathIntegratedAttenuation'].values[:] = 0.0


    return 

  
  def getRadarSpectrum(self):

    """
    Convert a spectrum of hydrometeor backscattering (per hydrometeor)
    as a function of size into a merged spectrum as a function of velocity. 
 
    Adds radar_spectrum to result object
    """
    assert 'pathIntegratedAttenuation' in self.results.keys()
    assert 'hydro_backSpec' in self._tmp.keys()
    assert len(self.dopplerVelBins) > 0

    if 'radar_spectrum' not in self.results.keys():

      coordsND = self.additionalDimList + [self._dataset.layer,  self.results.frequency,self.results.activePolarisation,self.results.dopplerVelBin]
      shapeND = tuple(map(len,coordsND))
      self.results['radar_spectrum'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'mm6/m3'})


      for ff, freq in enumerate(self.radarFrequencies):
        for pp, pol in enumerate(self.radarPolarisations):

            spec = pyPamtraRadarSimulator.simulateRadarSpectrum(
              self._dataset['hydroSize'].values[:],
              self._tmp['hydro_backSpec'].values[...,ff,pp,:,:],
              self._dataset.hydroMass.values[:],
              self._dataset.hydroRho.values[:],
              self._dataset.hydroCrossSectionArea.values[:],
              self.results['pathIntegratedAttenuation'].values[...,ff],
              self._dataset.height.values[:],
              self._dataset.temperature.values[:],
              self._dataset.pressure.values[:],
              self._dataset.verticalWind.values[:],
              self._tmp.specBroadening.values[...,ff],
              freq,
              fallVelocityRelation =self.settings['fallVelocityRelation'],
              radarMaxV =self.settings['radarProperties'][repr(freq)]['maxV'],
              radarMinV =self.settings['radarProperties'][repr(freq)]['minV'],
              radarAliasingNyquistInterv = self.settings['radarSimulator']['aliasingNyquistInterv'],
              radarNFFT = len(self.dopplerVelBins),
# to do: add non uniform beamfilling to pamtra again. 
              # radarAirmotion = False,
              # radarAirmotionModel = "step", #"constant","linear","step"
              # radarAirmotionVmin = -4.0,
              # radarAirmotionVmax = +4.0,
              # radarAirmotionLinearSteps = 30,
              # radarAirmotionStepVmin = 0.5,
              radarPNoise1000 = 10**(0.1*self.settings['radarProperties'][repr(freq)]['pNoise1000']),
              radarK2 =self.settings['radarProperties'][repr(freq)]['k2'], # dielectric constant |K|2 (always for liquid water by convention) for the radar equation
              radarNAve = self.settings['radarProperties'][repr(freq)]['nAve'],
              seed  = self.settings['radarSimulator']['randomSeed'],
              verbosity = self.settings['general']['verbosity']
                  )

            self.results.radar_spectrum[...,ff,pp,:] = spec
    return

  def getRadarMoments(self):

    assert 'radar_spectrum' in self.results.keys()


    if 'reflectivity' not in self.results.keys():

      coordsND = self.additionalDimList + [
        self._dataset.layer,  
        self.results.frequency,
        self.results.activePolarisation,
        self.results.dopplerPeak,
        ]
      shapeND = tuple(map(len,coordsND))


      for var, unit in [
        ('reflectivity','dBz'),
        ('meanDopplerVelocity','m/s'),
        ('dopplerSpectrumWidth','m/s'),
        ('skewness','-'),
        ('kurtosis','-'),
        ('leftSlope','dB s/m'),
        ('rightSlope','dB s/m'),
        ('leftEdge','m/s'),
        ('rightEdge','m/s'),
        ('signalToNoiseRatio','dB'),
        ]:
        self.results[var] = xr.DataArray(
        np.nan*np.ones(shapeND,dtype=np.float),
        coords=coordsND,
        attrs={'units':unit},
        )

      self.results['quality'] = xr.DataArray(
        -9999*np.ones(shapeND[:-1],dtype=int),
        coords=coordsND[:-1],
        attrs={
          'units':'bytes',
          'desciption':"1st byte: aliasing; 2nd byte: more peaks present; 7th: no peak found; 8th: principal peak isolated"
          },
      )


    for ff, freq in enumerate(self.radarFrequencies):

      output =  pyPamtraRadarMoments.calc_radarMoments(
        self.results.radar_spectrum.values[...,ff,:] ,
        maxV = self.settings['radarProperties'][repr(freq)]['maxV'],
        minV = self.settings['radarProperties'][repr(freq)]['minV'],
        noAve = self.settings['radarProperties'][repr(freq)]['nAve'],
        # noiseDistanceFactor = 0,
        specNoiseMax = None,
        specNoiseMean = None,
        nPeaks = len(self.results.dopplerPeak),
        # peakMinBins = 2,
        peakMinSnr = -10,
        smoothSpectrum = self.settings['radarProperties'][repr(freq)]['smoothSpectrum'],
        # useWiderPeak = False,
        verbose = self.settings['general']['verbosity'],
        )
      spectrum_out, moments, slope, edge, quality, noise = output

      # self.results['reflectivityFromSpectra'].values[...,ff,:,:] = 10*np.log10(np.sum(spectraCalib))
      self.results['reflectivity'].values[...,ff,:,:] = 10*np.log10(moments[...,0,:])
      self.results['meanDopplerVelocity'].values[...,ff,:,:] = moments[...,1,:]
      self.results['dopplerSpectrumWidth'].values[...,ff,:,:] = moments[...,2,:]
      self.results['skewness'].values[...,ff,:,:] = moments[...,3,:]
      self.results['kurtosis'].values[...,ff,:,:] = moments[...,4,:]
      self.results['leftSlope'].values[...,ff,:,:] = slope[...,0,:]
      self.results['rightSlope'].values[...,ff,:,:] = slope[...,1,:]
      self.results['leftEdge'].values[...,ff,:,:] = edge[...,0,:]
      self.results['rightEdge'].values[...,ff,:,:] = edge[...,1,:]
      self.results['quality'].values[...,ff,:] = quality

      #Store used settings in result object
      self.results.attrs['pyPamtra2 Settings'] = json.dumps(self.settings)

    return


# to do enable copyimg of settings
  # def sel(self,*args,**kwargs):
  #   """
  #   Select a subprofile based on labels. See xr.Dataset.sel for documentation
  #   """
  #   newSelf = deepcopy(self)
  #   newSelf._dataset = self._dataset.sel(*args,**kwargs)
  #   try:
  #     newSelf._tmp = self._tmp.sel(*args,**kwargs)
  #   except AttributeError:
  #     pass
  #   try:
  #     newSelf.results = self.results.sel(*args,**kwargs)
  #   except AttributeError:
  #     pass
  #   return newSelf

  # def isel(self,*args,**kwargs):
  #   """
  #   Select a subprofile based on index. See xr.Dataset.isel for documentation
  #   """
  #   newSelf = deepcopy(self)
  #   newSelf._dataset = self._dataset.isel(*args,**kwargs)
  #   try:
  #     newSelf._tmp = self._tmp.isel(*args,**kwargs)
  #   except AttributeError:
  #     pass
  #   try:
  #     newSelf.results = self.results.isel(*args,**kwargs)
  #   except AttributeError:
  #     pass
  #   return newSelf


