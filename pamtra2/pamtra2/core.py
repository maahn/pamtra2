# -*- coding: utf-8 -*-

from collections import OrderedDict
from copy import deepcopy
import numpy as np
import xarray as xr

from . import helpers
from . import hydrometeors
# import json
# from functools import lru_cache

__version__ = 0.2

class profile (xr.Dataset):

  def __init__(self, 
    nLayer,
    hydrometeors,
    hydrometeorBulkProperties = ['waterContent'],
    additionalDims={}
    ):
    

    coordsGeom = {**additionalDims, 'layer': range(nLayer)}

    coordsGeom = helpers.concatDicts(additionalDims,{'layer' : range(nLayer)})
    cordsHydro = OrderedDict(hydrometeor = hydrometeors)
    coordsGeomHydro = helpers.concatDicts(coordsGeom,cordsHydro)
    coordsAll = helpers.concatDicts(coordsGeomHydro,OrderedDict())

    super().__init__(coords=coordsAll)
    
    arrayVars = [
      ('height','m',coordsGeom,np.float64),
      ('temperature','K',coordsGeom,np.float64),
      ('pressure','Pa',coordsGeom,np.float64),
      ('relativeHumidity','%',coordsGeom,np.float64),
      ('horizontalWind','m/s',coordsGeom,np.float64),
      ('verticalWind','m/s',coordsGeom,np.float64),
      # ('waterContent','kg/m^3',coordsGeomHydro,np.float64),
      # ('effectiveRadius','m',coordsGeomHydro,np.float64),
      # ('hydrometeorNtot','1/m^3',coordsGeomHydro,np.float64),
      # ('hydrometeor','-',cordsHydro,'S128'),
    ]
    for hydrometeorBulkPropertiy in hydrometeorBulkProperties:
      arrayVars.append(
        (hydrometeorBulkPropertiy,'SI',coordsGeomHydro,np.float64)
        )

    for var,unit,coords,dtype in arrayVars:
        thisShape = tuple(map(len,coords.values()))
        self[var] = xr.DataArray(
              (np.zeros(thisShape)*np.nan).astype(dtype),
              coords=coords.values(),
              dims = coords.keys(),
              attrs={'unit':unit},
              )

    return 


class pamtra2(object):

  def __init__(
    self,
    nLayer,
    hydrometeors,
    hydrometeorBulkProperties = ['waterContent'],
    additionalDims={}
    ):

    self.profile = profile(
      nLayer,
      hydrometeors,
      hydrometeorBulkProperties,
      additionalDims
      )
    self.additionalDims = additionalDims
    self.hydrometeorBulkProperties = hydrometeorBulkProperties
    self.hydrometeors = OrderedDict()
    for hh in hydrometeors:
      self.hydrometeors[hh] = None

    self.instruments = OrderedDict()

    return

  @property
  def nHydrometeors(self):
    return len(self.hydrometeors)

  @property
  def nInstruments(self):
    return len(self.instruments)

  @property
  def nLayer(self):
    return len(self.profile.layer)


  def describeHydrometeor(
    self,
    hydrometeorConstructor,
    **kwargs,
  ):

    name = kwargs['name']
    self.hydrometeors[name] = hydrometeorConstructor(
      self,
      **kwargs
      )

    self.hydrometeors[name].calculateProperties()

    return self.hydrometeors[name]

  def addInstrument(
    name,
    frequencies = [],
    ):

    self.frequencies = frequencies
    self.instrumentName = name

    for hh in self.hydrometeors.keys():
      self.hydrometeors[hh].frequencies = frequencies

