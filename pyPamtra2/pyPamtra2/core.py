# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function
from builtins import super 

from collections import OrderedDict
from copy import deepcopy
import numpy as np
import xarray as xr
import json

import pyPamtraRadarSimulator
import pyPamtraRadarMoments
import refractive

from . import configuration
from . import helpers
from . import decorators
from . import constants
from . import sizeDistributions

__version__ = '0.1'


class pyPamtra2(object):

  def  __init__(
    self,
    nLayer = 1,
    hydrometeors = ['hydrometeor_1'],
    nHydroBins = 0,
    frequencies=[35.5],
    activePolarisations = ['NN'],
    nFFT = 0,
    nPeaks = 1,
    settings = None,
    additionalDims = {}
    ):


    self.nLayer =nLayer 
    self.hydrometeors =hydrometeors
    self.nHydro =len(hydrometeors) 
    self.nHydroBins =nHydroBins 
    self.frequencies=np.array(frequencies)
    self.activePolarisations =activePolarisations 
    self.nFFT =nFFT 
    self.nPeaks =nPeaks 
    self.additionalDims =OrderedDict(additionalDims)
    self.dopplerVelBins = range(nFFT)

    if settings is None:
      self.settings = configuration.Settings(self.frequencies,self.hydrometeors)
    else:
      self.settings = configuration.Settings(settings)

    self._createXrObjects()
    self._fillDataObject()


    return


  def _createXrObjects(self):

    coordsData = deepcopy(self.additionalDims)
    coordsData.update(OrderedDict(
      layer = range(self.nLayer),
      hydrometeor = self.hydrometeors,
      hydroBin = range(self.nHydroBins),
      ))
    self.data = xr.Dataset(coords=coordsData)

    cordsTmp = coordsData.copy()
    cordsTmp.update(OrderedDict(
      frequency = self.frequencies,
      activePolarisation = self.activePolarisations
      ))
    self._tmp = xr.Dataset(coords=cordsTmp)

    cordsRes = cordsTmp.copy()
    cordsRes.update(OrderedDict(
      dopplerVelBin = range(self.nFFT),
      dopplerPeak = range(self.nPeaks),
      activePolarisations = self.activePolarisations,
      ))
    self.results = xr.Dataset(coords=cordsRes)

    return


  def _fillDataObject(self):

    #for convinience, store additonal variables
    self.additionalDimList = []
    for cc in self.additionalDims.keys():
      self.additionalDimList.append((self.data[cc]))


    self.coordsGeometry = self.additionalDimList + [self.data.layer]
    self.coordsGeometryHydro = self.coordsGeometry + [self.data.hydrometeor]
    self.coordsGeometryHydroSpectrum = self.coordsGeometryHydro + [self.data.hydroBin]


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
        self.data[var] = xr.DataArray(
          np.zeros(thisShape)*np.nan,
          coords=coords,
          attrs={'unit':unit},
          )

    self.data['hydro_name'] = xr.DataArray(np.arange(len(self.data.hydrometeor)).astype('S10'),coords=[self.data.hydrometeor])

    return


  def addHydrometeorBinDimension(self,nHydroBins):

    """
    In case you want to add it later due to memory concerns
    """

    self.data = self.data.reindex(hydroBin = range(nHydroBins))
    self.nHydroBins = nHydroBins
    return 

  # def addLiquidHydrometeor(self,hydroIndex):
  #   """
  #   we need something smarter here!
  #   """

  #   self.data['hydroRho'][...,hydroIndex,:] = constants.rhoWater
  #   self.data['hydroMass'][...,hydroIndex,:] = (self.data['hydroSize'][...,hydroIndex,:]/2.)**3 * np.pi * 4./3. * constants.rhoWater
  #   self.data['hydroCrossSectionArea'][...,hydroIndex,:] = (self.data['hydroSize'][...,hydroIndex,:]/2.)**2 * np.pi

  #   return

  def addHydrometeorProperties(
    self,
    hydroIndex,
    diameters,
    psd,
    hydroType = 'liquid',
    useLogSizeSpacing = False,
    psdNormalized = False,
    psdFuncArgs = [],
    psdFuncKwArgs = {},
    mass_size_a=np.pi * 4./3. * constants.rhoWater,
    mass_size_b=3.,
    area_size_a=np.pi,
    area_size_b=2.,
    hydrometeorProperties = {},
    ):
    """
    diameters : array_like
      if len == 2 from Dmin to Dmax
    mass and area default to properties of liquid water
    """

    assert self.nHydroBins >0, 'First use addHydrometeorBinDimension to create hydrometeor bin dimension.'
    assert hydroType in ['liquid','ice','snow']

    if isinstance(hydroIndex,str):
      hydroName = deepcopy(hydroIndex)
      hydroIndex = self.hydrometeors.index(hydroIndex)
    else:
      hydroIndex = deepcopy(hydroIndex)
      hydroName = self.hydrometeors[hydroIndex]


    if len(diameters) ==2 and not useLogSizeSpacing:
      self.data.hydroSize.values[...,hydroIndex,:] = np.linspace(diameters[0],diameters[1],self.nHydroBins)
    elif len(diameters) ==2 and useLogSizeSpacing:
      self.data.hydroSize.values[...,hydroIndex,:] = np.logspace(np.log10(diameters[0]),np.log10(diameters[1]),self.nHydroBins)
    else:
      self.data.hydroSize.values[...,hydroIndex,:] = diameters
    # self.data.hydroSizeBinWidth.values[...,hydroIndex,:] = np.gradient(self.data.hydroSize.values[...,hydroIndex,:])


    #either explicit distribution or 
    if callable(psd):
      self.data.hydroPSD.values[...,hydroIndex,:] = psd(self,hydroIndex,*psdFuncArgs,**psdFuncKwArgs)
    else:
      self.data.hydroPSD.values[...,hydroIndex,:] = psd

    #roll back normalization
    if psdNormalized:
      self.data.hydroPSD.values[...,hydroIndex,:] = self.data.hydroPSD.values[...,hydroIndex,:]/np.gradient(self.data.hydroSize.values[...,hydroIndex,:],axis=-1)



    self.data.hydroRho.values[...,hydroIndex,:] = constants.rhoWater
    self.data.hydroMass.values[...,hydroIndex,:] = (self.data.hydroSize.values[...,hydroIndex,:]/2.)**3 * np.pi * 4./3. * constants.rhoWater
    self.data.hydroCrossSectionArea.values[...,hydroIndex,:] = (self.data.hydroSize.values[...,hydroIndex,:]/2.)**2 * np.pi

    self.settings['hydrometeorProperties'][hydroName]['type'] = hydroType
    self.settings['hydrometeorProperties'][hydroName].update(
      configuration.DEFAULT_HYDROMETEOR_PROPERTIES_BY_TYPE[hydroType]
      )
    for key,value in hydrometeorProperties.items():
      self.settings['hydrometeorProperties'][hydroName][key] = value

    return


  def getSpectralBroadening(self,method='pyPamtraRadarSimulator'):

    assert not np.all(np.isnan(self.data.eddyDissipationRate.values)), 'found nan in eddyDissipationRate'
    assert not np.all(np.isnan(self.data.horizontalWind.values)), 'found nan in horizontalWind'
    assert not np.any(np.isnan(self.data.height.values)), 'found nan in height'
    assert method == 'pyPamtraRadarSimulator'

    coordsND = self.additionalDimList + [self._tmp.layer, self._tmp.frequency]
    shapeND = tuple(map(len,coordsND))
    if 'specBroadening' not in self._tmp.keys():
      self._tmp['specBroadening'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'m/s'})

    for ff, freq in enumerate(self.frequencies):

      specBroadening = pyPamtraRadarSimulator.calcSpectralBraodening(
        self.data.eddyDissipationRate.values[:],
        self.data.horizontalWind.values[:],
        self.data.height.values[:],
        self.settings['radarProperties'][freq]['fwhrBeamwidthDeg'],
        self.settings['radarProperties'][freq]['integrationTime'],
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

    Temperatures  = self.data['temperature'].values
    Frequencies = self.frequencies * 1e9

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
          Densities = self.data['hydroRho'].sel(hydrometeor=hydroName).values
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


    for ff, freq in enumerate(self.frequencies):
      for pp, pol in enumerate(self.activePolarisations):
        for hh in range(self.nHydro):
          if method == 'Rayleigh':
            assert pol == 'NN'

            back_spec = helpers.rayleigh(
              self.data['hydroSize'].values[...,hh,:], 
              self._tmp['K2'].values[...,ff,pp,hh], 
              freq
              )
            self._tmp['hydro_backSpec'].values[...,ff,pp,hh,:] = back_spec * self.data['hydroPSD'][...,hh,:].values

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

      coordsND = self.additionalDimList + [self.data.layer,  self.results.frequency,self.results.activePolarisation,self.results.dopplerVelBin]
      shapeND = tuple(map(len,coordsND))
      self.results['radar_spectrum'] = xr.DataArray(np.zeros(shapeND)*np.nan,coords=coordsND,attrs={'unit':'mm6/m3'})


      for ff, freq in enumerate(self.frequencies):
        for pp, pol in enumerate(self.activePolarisations):

            spec = pyPamtraRadarSimulator.simulateRadarSpectrum(
              self.data['hydroSize'].values[:],
              self._tmp['hydro_backSpec'].values[...,ff,pp,:,:],
              self.data.hydroMass.values[:],
              self.data.hydroRho.values[:],
              self.data.hydroCrossSectionArea.values[:],
              self.results['pathIntegratedAttenuation'].values[...,ff],
              self.data.height.values[:],
              self.data.temperature.values[:],
              self.data.pressure.values[:],
              self.data.verticalWind.values[:],
              self._tmp.specBroadening.values[...,ff],
              freq,
              fallVelocityRelation =self.settings['fallVelocityRelation'],
              radarMaxV =self.settings['radarProperties'][freq]['maxV'],
              radarMinV =self.settings['radarProperties'][freq]['minV'],
              radarAliasingNyquistInterv = self.settings['radarSimulator']['aliasingNyquistInterv'],
              radarNFFT = len(self.dopplerVelBins),
# to do: add non uniform beamfilling to pamtra again. 
              # radarAirmotion = False,
              # radarAirmotionModel = "step", #"constant","linear","step"
              # radarAirmotionVmin = -4.0,
              # radarAirmotionVmax = +4.0,
              # radarAirmotionLinearSteps = 30,
              # radarAirmotionStepVmin = 0.5,
              radarPNoise1000 = 10**(0.1*self.settings['radarProperties'][freq]['pNoise1000']),
              radarK2 =self.settings['radarProperties'][freq]['k2'], # dielectric constant |K|2 (always for liquid water by convention) for the radar equation
              radarNAve = self.settings['radarProperties'][freq]['nAve'],
              seed  = self.settings['radarSimulator']['randomSeed'],
              verbosity = self.settings['general']['verbosity']
                  )

            self.results.radar_spectrum[...,ff,pp,:] = spec
    return

  def getRadarMoments(self):

    assert 'radar_spectrum' in self.results.keys()


    if 'reflectivity' not in self.results.keys():

      coordsND = self.additionalDimList + [
        self.data.layer,  
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


    for ff, freq in enumerate(self.frequencies):

      output =  pyPamtraRadarMoments.calc_radarMoments(
        self.results.radar_spectrum.values[...,ff,:] ,
        maxV = self.settings['radarProperties'][freq]['maxV'],
        minV = self.settings['radarProperties'][freq]['minV'],
        noAve = self.settings['radarProperties'][freq]['nAve'],
        # noiseDistanceFactor = 0,
        specNoiseMax = None,
        specNoiseMean = None,
        nPeaks = len(self.results.dopplerPeak),
        # peakMinBins = 2,
        peakMinSnr = -10,
        smoothSpectrum = self.settings['radarProperties'][freq]['smoothSpectrum'],
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
  #   newSelf.data = self.data.sel(*args,**kwargs)
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
  #   newSelf.data = self.data.isel(*args,**kwargs)
  #   try:
  #     newSelf._tmp = self._tmp.isel(*args,**kwargs)
  #   except AttributeError:
  #     pass
  #   try:
  #     newSelf.results = self.results.isel(*args,**kwargs)
  #   except AttributeError:
  #     pass
  #   return newSelf


