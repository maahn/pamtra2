# -*- coding: utf-8 -*-
# from __future__ import division, absolute_import, print_function
# from builtins import super 

from copy import deepcopy
import json



class SettingsSection(dict):

    """
    Like a dict, but with .freez() method and .__frozen attribute. If frozen, no more
    keys can be added.
    """

    def __init__(self,*args,**kwargs):
      super().__init__(*args,**kwargs)
      self.__frozen = False

    def freeze (self):
        self.__frozen = True

    def unfreeze (self):
        self.__frozen = False

    def __setitem__ (self, key, value):
      if self.__frozen and key not in self:
        raise ValueError('Frozen dictionary does not contain key %s'%key)
      super().__setitem__(key, value)


DEFAULT_RADAR_SIMULATOR = SettingsSection({
  'randomSeed' : 0,
  "aliasingNyquistInterv" : 1,
  })

DEFAULT_RADAR_PROPERTIES =  SettingsSection({
  'integrationTime'   :  1.4,  # MMCR Barrow during ISDAC
  'fwhrBeamwidthDeg' : 0.31,  # full width half radiation beam width MMCR Barrow 
  'k2' : 0.92,  # dielectric constant |K|² (always for liquid water by convention) for the radar equation
  "nAve" : 150,#: number of average spectra for noise variance reduction, typical range [1 150]
  "maxV" : 7.885,#: MinimumNyquistVelocity in m/sec
  "minV" : -7.885,#: MaximumNyquistVelocity in m/sec
  "pNoise1000": -32.23, # dBz in 1 km height, mean value for BArrow MMCR during ISDAC
  "smoothSpectrum": True, # Smooth spectrum before estimating moments
  })

DEFAULT_GENERAL = SettingsSection({
  'verbosity' : 0,
  })


DEFAULT_SETTINGS = SettingsSection({
  'general' : DEFAULT_GENERAL,
  'fallVelocityRelation' : 'khvorostyanov01_spheres',#'heymsfield10_particles', #move
  'radarProperties':SettingsSection({}),
  'radarSimulator':DEFAULT_RADAR_SIMULATOR,
})



class Settings(SettingsSection):
    def __init__(self,arg):
      """
      Create PyPamtra2 Settings

      if type(arg) is str: load settings from JSON file
      if type(arg) is dict: get settings from arguments
      if type(arg) is lift of frequencies: populate settings with default values including radar defualt for each frequency.
      """
      if isinstance(arg, dict):
        super().__init__(arg)
      elif isinstance(arg, str):
        with open(arg, 'r') as jsonfile:
          settings = json.load(jsonfile,object_hook=SettingsSection)
        super().__init__(settings)
      else:
        super().__init__(deepcopy(DEFAULT_SETTINGS))
        for freq in arg:
          self['radarProperties'][freq] = deepcopy(DEFAULT_RADAR_PROPERTIES)
          self['radarProperties'][freq].freeze()

      for k in self.keys():
        try:
          self[k].freeze()
        except AttributeError:
          print (k)
          pass

    def __setitem__ (self, key, value):
      raise ValueError('Only children of this dictionary can be modified')


    def to_json(self,fname):
      with open(fname, 'w') as jsonfile:
        json.dump(self,jsonfile)



