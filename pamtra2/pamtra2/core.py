# -*- coding: utf-8 -*-

from collections import OrderedDict
import numpy as np
import xarray as xr

from . import helpers
from . import units
from . import dimensions
from . import constants

__version__ = 0.2


class customProfile (xr.Dataset):
    """ """
    def __init__(
      self,
      parent,
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
            super().__init__(coords=parent.coords['all'])

            for var, coord, dtype in profileVars:
                coord = helpers.concatDicts(
                  *map(lambda x: parent.coords[x], coord)
                  )
                thisShape = tuple(map(len, coord.values()))
                self[var] = xr.DataArray(
                    (np.zeros(thisShape)*np.nan).astype(dtype),
                    coords=coord.values(),
                    dims=coord.keys(),
                    attrs={'unit': units.units[var]},
                )

            self['wavelength'] = constants.speedOfLight/self.frequency

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
          'eddyDissipationRate',
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

        self.coords = {}
        self.coords[dimensions.ADDITIONAL] = OrderedDict(additionalDims)
        self.coords[dimensions.LAYER] = OrderedDict(layer=range(nLayer))
        self.coords[dimensions.HYDROMETEOR] = OrderedDict(
          hydrometeor=hydrometeors)
        self.coords[dimensions.FREQUENCY] = OrderedDict(
          frequency=frequencies)

        self.coords['all'] = helpers.concatDicts(
          self.coords[dimensions.ADDITIONAL],
          self.coords[dimensions.LAYER],
          self.coords[dimensions.HYDROMETEOR],
          self.coords[dimensions.FREQUENCY],
          )


        self.profile = customProfile(
          self,
          profileVars=profileVars,
          nLayer=nLayer,
          hydrometeors=hydrometeors,
          frequencies=frequencies,
          additionalDims=additionalDims
        )
        self.additionalDims = additionalDims
        self.hydrometeors = helpers.AttrDict()
        for hh in hydrometeors:
            self.hydrometeors[hh] = None

        self.instruments = OrderedDict()

        return

    def getProfileAllBroadcasted(self, variables=None, sel={}):
        if variables is None:
            return xr.broadcast(self.profile.sel(**sel))[0]
        else:
            return xr.broadcast(self.profile.sel(**sel)[variables])[0]

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
