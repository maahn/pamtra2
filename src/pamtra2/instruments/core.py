# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr

from .. import helpers, units
from ..libs import pamgasabs


class instrument(object):
    '''Base instrument simulator class

    This is the base class for all instrument simulators providing
    the most important data structures.

    Parameters
    ----------
    parent : {pamtra2 class}
        Calling parent object
    frequencies : {list} or {'all'}, optional
        Use this instrument for the frequencies indicated in the list or
        'all' frequencies (the default is 'all', which means all frequencies)
    **settings : {dict}
        Additional settings for the instrument

    Attributes
    ----------
    frequencies : {list}
        List of used frequencies
    settings : {dict}
        Instrument specific settings
    parent : {pamtra2 class}
        Calling parent object
    profile : {xr.Dataset}
        Calling parent's profile
    hydrometeorProfiles : {pamtra2.helpers.AttrDict}
         Calling parent's hydrometeor profiles
    results : {xr.Dataset}
        Instrument simulator results
    '''

    def __init__(
        self,
        parent=None,
        frequencies='all',
        name='instrument',
        **settings
    ):
        if not hasattr(frequencies, '__iter__'):
            frequencies = [frequencies]
        self.frequencies = frequencies
        self.settings = settings
        self.parent = parent
        self.name = name

        self.results = xr.Dataset()

    def _link_parent(self):

        if self.parent is None:
            raise AttributeError('set .parent attribute to pamtra2 object')

        if (type(self.frequencies) is str) and (self.frequencies == 'all'):
            self.frequencies = self.parent.profile.frequency

        self.profile = self.parent.profile.sel(frequency=self.frequencies)
        self.hydrometeorProfiles = helpers.AttrDict()
        for hydro in self.parent.hydrometeors.keys():
            self.hydrometeorProfiles[hydro] = self.parent.hydrometeors[
                hydro].profile.sel(frequency=self.frequencies)

    def to_netcdf(
        self,
        fname,
        resultVariables='all',
        profileVariables=['height'],
        **kwargs
    ):
        '''Store result as netcdf file

        Parameters
        ----------

        fname : str
            filename
        resultVariables : 'all' or list of str, optional
            variales to store in result (default 'all')
        profileVariables : 'all' or list of str, optional
            profile variables to include in file (default ['height'])
        **kwargs : kwargs to pass to xarray's to_netcdf function.

        '''
        if resultVariables == 'all':
            results = self.results
        else:
            results = self.results[resultVariables]

        if profileVariables == 'all':
            profile = self.parent.profile
        else:
            profile = self.parent.profile[profileVariables]

        results.merge(profile).to_netcdf(fname, **kwargs)


class microwaveInstrument(instrument):
    def __init__(
        self,
        parent=None,
        frequencies='all',
        gaseousAttenuationModel='Rosenkranz98',
        **settings
    ):
        super().__init__(
            parent=parent,
            frequencies=frequencies,
            gaseousAttenuationModel=gaseousAttenuationModel,
            **settings,
        )

    def _calcHydrometeorAbsorption(self):
        '''Calculate hydrometeor absorption

        Returns
        -------
        xr.DataArray
            absorption coefficient
        '''
        if self.parent.nHydrometeors>0:

            hydroAbs = self.parent.getIntegratedScatteringCrossSections(
                crossSections=['extinctionCrossSection'],
                frequencies=self.frequencies,
            )['extinctionCrossSection']

        else:
            hydroAbs = 0

        return hydroAbs

    def _calcGaseousAbsorption(self):
        '''Calculate gaseous absorption

        Returns
        -------
        xr.DataArray
            absorption coefficient
        '''
        coords = helpers.concatDicts(
            self.parent.coords['additional'],
            self.parent.coords['layer'],
            self.parent.coords['frequency']
        )
        thisProf = self.parent.profile.sel(
            frequency=self.frequencies
            ).stack(merged=coords)

        kwargs = {}

        args = [thisProf.frequency,
                thisProf.temperature,
                thisProf.waterVaporPressure,
                thisProf.pressure
                ]
        args = xr.broadcast(*args)

        if self.settings['gaseousAttenuationModel'] == 'Rosenkranz98':
            kwargs['sumResults'] = True
            func = pamgasabs.calculate_gas_absorption_rosenkranz98
        elif self.settings['gaseousAttenuationModel'] == 'Liebe93':
            func = pamgasabs.calculate_gas_absorption_liebe93
        else:
            raise ValueError('Do not recognize gaseousAttenuationModel: %s' %
                             self.settings['gaseousAttenuationModel'])

        gasAbs = xr.apply_ufunc(
            func,
            *args,
            kwargs=kwargs,
            dask='parallelized',
            output_dtypes=[args[1].dtype],
        )

        gasAbs = xr.Dataset({'gasAbs': gasAbs})
        gasAbs = helpers.xrFastUnstack(gasAbs, 'merged').gasAbs

        return gasAbs
