# -*- coding: utf-8 -*-
import numpy as np
import collections
from copy import deepcopy
import xarray as xr
import inspect
from functools import wraps
from collections import OrderedDict

from . import decorators


@decorators.NDto2DtoND(
    referenceIn=0, noOfInDimsToKeep=1, convertInputs=[0, 1],
    convertOutputs=[0], verbosity=10
)
def rayleigh(diameter, K2, frequency):
    """
    To do: move this routine into separte module together with t-matrix,
    mie etc code.
    """
    C = 299792458.

    K2 = np.asarray(K2)
    diameter = np.asarray(diameter)

    wavelength = C / (frequency*1e9)
    prefactor = np.pi**5 * K2 / wavelength**4
    back_spec = prefactor[:, np.newaxis] * diameter**6

    return back_spec


def provideArgKwargNames(func):
    '''Provide all arguments and keyword arguments of a function. Does not
    work with *arg and **kwarg
    '''
    argNames = []
    kwargNames = []
    insp = inspect.signature(func)
    for p in insp.parameters:
        default = insp.parameters[p].default
        if default is insp.empty:
            argNames.append(p)
        else:
            kwargNames.append(p)
    return argNames, kwargNames


def concatDicts(*dicts):
    """
    concatenate (ordered) dicts
    """
    dMerged = deepcopy(dicts[0])
    for dd in range(1, len(dicts)):
        dMerged.update(dicts[dd])
    return dMerged


def swapListItems(li, i1, i2):
    '''
    Swap i1 and i2 of li.
    '''
    li = deepcopy(li)
    a = li.index(i1)
    b = li.index(i2)
    li[b], li[a] = li[a], li[b]

    return li


class AttrDict(collections.OrderedDict):
    """Dictionary accesible through attributes"""

    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def dimensionToVariables(darray, dimension, variables):
    profile = {}
    for ii, var in enumerate(variables):
        profile[var] = darray.isel(**{dimension: ii})
    profile = xr.Dataset(profile)
    return profile


def getInputCoreDims(args, coreDims):
    input_core_dims = []
    for arg in args:
        this_core_dims = []
        for coreDim in coreDims:
            if coreDim in arg.dims:
                this_core_dims.append(coreDim)
        input_core_dims.append(this_core_dims)
    return input_core_dims


def mergeFlattResults(
    core_shapes=None
):
    """Decorator to merge multiple outputs of a function's result
    into one joined array. Dimensions NOT in core_shapes are preserved.
    """
    def mergeFlattResults_decorator(func):
        @wraps(func)
        def flatt_func(*args, **kwargs):
            result = func(*args, **kwargs)
            if type(result) is not tuple:
                return result
            else:
                result = list(result)
                for ii in range(len(result)):
                    res_shape = list(result[ii].shape)
                    for sh in core_shapes[ii]:
                        res_shape.remove(sh)
                    new_shape = tuple(
                        res_shape) + tuple([np.product(
                            core_shapes[ii]).astype(int)])
                    result[ii] = result[ii].reshape(new_shape)
                result = np.concatenate(result, axis=-1)
            return result
        return flatt_func
    return mergeFlattResults_decorator


def apply_ufunc_extended(
    func,
    *args,
    **kwargs
):
    """Extended version of xarray's `xr.apply_ufunc` which can handle with
    multiple output of a functions and dask. New keyword output_names
    required with a list of the names of the returned variables. Note that
    different types in the results are not preserved.
    """

    input_core_dims = kwargs.get('input_core_dims', [])
    output_sizes = kwargs.pop('output_sizes', {})
    output_names = kwargs.pop('output_names', [])
    output_core_dims = kwargs.pop('output_core_dims', [])

    OrderedDict(zip(output_names, output_core_dims))

    all_input_core_dims = []
    for input_core_dim in input_core_dims:
        all_input_core_dims.append(*input_core_dim)
    all_input_core_dims = list(set(all_input_core_dims))

    non_core_dims = []
    for arg in args:
        non_core_dims += [d for d in arg.dims if d not in all_input_core_dims]
    non_core_dims = list(set(non_core_dims))
    arg0 = args[0]._to_temp_dataset()

    non_core_shape = tuple([arg0.dims[key] for key in non_core_dims])

    for ii in range(len(output_core_dims)):
        if type(output_core_dims[ii]) is not tuple:
            output_core_dims[ii] = tuple([output_core_dims[ii]])

    output_shapes = OrderedDict()
    output_core_shapes = []
    for key, val in zip(output_names, output_core_dims):
        output_shapes[key] = non_core_shape + \
            tuple([output_sizes[ii] for ii in val])
        output_core_shapes.append(tuple([output_sizes[ii] for ii in val])
                                  )

    len_tmp_dim = np.sum([np.prod(vv) for vv in output_shapes.values()])

    flat_func = mergeFlattResults(output_core_shapes)(func)

    results = xr.apply_ufunc(
        flat_func,
        *args,
        output_core_dims=[('<TEMP_DIM>',)],
        output_sizes={'<TEMP_DIM>': len_tmp_dim},
        **kwargs,
    )

    ii = 0
    this_result = OrderedDict()
    for kk, (key, this_shape) in enumerate(
            zip(output_names, output_core_shapes)):
        new_shape = non_core_shape + this_shape
        this_len = np.prod(this_shape).astype(int)
        this_result[key] = results.data[..., ii:ii+this_len].reshape(new_shape)
        this_ccords = [arg0.coords[k] for k in non_core_dims] + [
            xr.DataArray(
                range(output_sizes[cc]),
                name=cc,
                dims=[cc]
                ) for cc in output_core_dims[kk]]

        this_result[key] = xr.DataArray(this_result[key], coords=this_ccords)
        ii += this_len
    return xr.Dataset(this_result)


def xrGradient(data, dimension=None):
    '''
    Wrapper for np.gradient which is not available in xarray
    '''
    if dimension is None:
        axis = 0
    else:
        axis = data.get_axis_num(dimension)
    return xr.DataArray(np.gradient(data, axis=axis), coords=data.coords)
