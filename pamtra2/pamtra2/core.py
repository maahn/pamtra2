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

  def __init__(self, nLayer,hydrometeors,additionalDims={}):
    

    coordsGeom = {**additionalDims, 'layer': range(nLayer)}

    coordsGeom = helpers.concatDicts(additionalDims,{'layer' : range(nLayer)})
    cordsHydro = OrderedDict(hydrometeor = hydrometeors)
    coordsGeomHydro = helpers.concatDicts(coordsGeom,cordsHydro)
    coordsAll = helpers.concatDicts(coordsGeomHydro,OrderedDict())

    super().__init__(coords=coordsAll)
    
    for var,unit,coords,dtype in [
      ('height','m',coordsGeom,np.float64),
      ('temperature','K',coordsGeom,np.float64),
      ('pressure','Pa',coordsGeom,np.float64),
      ('relativeHumidity','%',coordsGeom,np.float64),
      ('horizontalWind','m/s',coordsGeom,np.float64),
      ('verticalWind','m/s',coordsGeom,np.float64),
      ('hydrometeorWaterContent','kg/m^3',coordsGeomHydro,np.float64),
      ('hydrometeorEffectiveRadius','m',coordsGeomHydro,np.float64),
      ('hydrometeorNtot','1/m^3',coordsGeomHydro,np.float64),
      # ('hydrometeor','-',cordsHydro,'S128'),
    ]:
        thisShape = tuple(map(len,coords.values()))
        self[var] = xr.DataArray(
              (np.zeros(thisShape)*np.nan).astype(dtype),
              coords=coords.values(),
              dims = coords.keys(),
              attrs={'unit':unit},
              )

    return 


class pamtra2(object):

  def __init__(self,nLayer,hydrometeors,additionalDims={}):
    self.profile = profile(nLayer,hydrometeors,additionalDims)
    self.additionalDims = additionalDims

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
    name,
    kind = None, #liquid, ice
    nBins = None,
    sizeCenter = None,
    sizeDistribution = None,
    aspectRatio = None,
    mass = None,
    density = None,
    crossSectionArea = None,
  ):

    self.hydrometeors[name] = hydrometeors.properties(
        self,
        name = name, #or None, then str(index)
        kind = kind, #liquid, ice
        nBins = nBins,
        sizeCenter = sizeCenter,
        sizeDistribution = sizeDistribution,
        aspectRatio = aspectRatio,
        mass = mass,
        density = density,
        crossSectionArea = crossSectionArea,
        )

    return self.hydrometeors[name]

  def addInstrument(
    name,
    frequencies = [],
    ):

    self.frequencies = frequencies
    self.instrumentName = name

    for hh in self.hydrometeors.keys():
      self.hydrometeors[hh].frequencies = frequencies

