# -*- coding: utf-8 -*-

from collections import OrderedDict
import numpy as np
import xarray as xr

from . import helpers
from . import units
from . import dimensions

__version__ = 0.2


class customProfile (xr.Dataset):
    """ """
    def __init__(
      self,
      profileVars=[],
      nLayer=None,
      hydrometeors=None,
      frequencies=None,
      additionalDims={},
    ):

        if isinstance(profileVars, xr.Dataset):
            self = profileVars
            return
        else:

            coords = {}
            coords[dimensions.ADDITIONAL] = OrderedDict(additionalDims)
            coords[dimensions.LAYER] = OrderedDict(layer=range(nLayer))
            coords[dimensions.HYDROMETEOR] = OrderedDict(
              hydrometeor=hydrometeors)
            coords[dimensions.FREQUENCY] = OrderedDict(frequency=frequencies)

            coordsAll = helpers.concatDicts(
              coords[dimensions.ADDITIONAL],
              coords[dimensions.LAYER],
              coords[dimensions.HYDROMETEOR],
              coords[dimensions.FREQUENCY],
              )

            super().__init__(coords=coordsAll)

            for var, coord, dtype in profileVars:
                coord = helpers.concatDicts(*map(lambda x: coords[x], coord))
                thisShape = tuple(map(len, coord.values()))
                self[var] = xr.DataArray(
                    (np.zeros(thisShape)*np.nan).astype(dtype),
                    coords=coord.values(),
                    dims=coord.keys(),
                    attrs={'unit': units.units[var]},
                )

            return


class pamtra2(object):
    """ """

    def __init__(
      self,
      profileVars=[
        (
          'height',
          [dimensions.ADDITIONAL, dimensions.LAYER],
          np.float64
          ),
        (
          'temperature',
          [dimensions.ADDITIONAL, dimensions.LAYER],
          np.float64
          ),
        (
          'pressure',
          [dimensions.ADDITIONAL, dimensions.LAYER],
          np.float64
          ),
        (
          'relativeHumidity',
          [dimensions.ADDITIONAL, dimensions.LAYER],
          np.float64
          ),
        (
          'horizontalWind',
          [dimensions.ADDITIONAL, dimensions.LAYER],
          np.float64
          ),
        (
          'verticalWind',
          [dimensions.ADDITIONAL, dimensions.LAYER],
          np.float64
          ),
        (
          'waterContent',
          [dimensions.ADDITIONAL, dimensions.LAYER, dimensions.HYDROMETEOR],
          np.float64
          ),
      ],
      nLayer=None,
      hydrometeors=None,
      frequencies=None,
      additionalDims={},
    ):

        self.profile = customProfile(
          profileVars=profileVars,
          nLayer=nLayer,
          hydrometeors=hydrometeors,
          frequencies=frequencies,
          additionalDims=additionalDims
        )
        self.additionalDims = additionalDims
        self.hydrometeors = OrderedDict()
        for hh in hydrometeors:
            self.hydrometeors[hh] = None

        self.instruments = OrderedDict()

        return

    @property
    def nHydrometeors(self):
        """ """
        return len(self.hydrometeors)

    @property
    def nInstruments(self):
        """ """
        return len(self.instruments)

    @property
    def nLayer(self):
        """ """
        return len(self.profile.layer)

    def describeHydrometeor(
        self,
        hydrometeorClass,
        **kwargs,
    ):
        """
        Add hydrometeor properties for one hydrometeor. Hydrometeor is added
        to self.hydrometeors[name]

        Parameters
        ----------
        hydrometeorClass :
            hydrometeor class
        **kwargs : dict
            Arguments handed over to hydrometeorClass to initialize class

        Returns
        -------
        hydrometeorClass :
            Evaluated hydrometeor class
        """

        name = kwargs['name']
        self.hydrometeors[name] = hydrometeorClass(
            self,
            **kwargs
        )

        self.hydrometeors[name].calculateProperties()

        return self.hydrometeors[name]

    def addInstrument(
        self,
        name,
        frequencies=[],
    ):
        """

        Parameters
        ----------
        name :

        frequencies :
             (Default value = [])

        Returns
        -------

        """

        self.frequencies = frequencies
        self.instrumentName = name

        for hh in self.hydrometeors.keys():
            self.hydrometeors[hh].frequencies = frequencies
