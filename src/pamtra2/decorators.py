# -*- coding: utf-8 -*-
import pickle
from copy import deepcopy
from functools import wraps

import numpy as np


class MemoizeMutable:
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}

    def __call__(self, *args, **kwds):
        str = pickle.dumps(args, 1)+pickle.dumps(kwds, 1)
        if not self.memo.has_key(str):
            print("miss")  # DEBUG INFO
            self.memo[str] = self.fn(*args, **kwds)
        else:
            print("hit")  # DEBUG INFO
        return self.memo[str]


def NDto2DtoND(
  referenceIn=0, noOfInDimsToKeep=1, convertInputs=[0], convertOutputs=[0],
  verbosity=0
  ):
    """
    Decorator to turn the Pamtra functions which expect only one dimension
    (usually height) into function accepting arbritrary shapes.

    Parameters
    ----------

    referenceIn : int, optional
        Input argument used as reference for reshaping (default 0)
    noOfInDimsToKeep : int, optional
        number of dimensions to preserve of referenceIn (counting from the
        back) (default 1)
    convertInputs : list of int, optional
        list of indices of input arguments which should be treated. (default
        [0])
    convertOutputs : list of int, optional
        list of indices of output arguments which should be treated. (default
        [0])
    verbosity : int, optional
        verbosity level (default 0)
    Returns
    -------

    function : function
      decorated function

    """
    def NDto2DtoND_decorator(func):
        @wraps(func)
        def inner(*args, **kwargs):
            if verbosity > 0:
                print('decorating %s' % func.__name__, referenceIn,
                      noOfInDimsToKeep, convertInputs, convertOutputs)
            referenceShape = None
            addDims = False
            args = list(args)
            # make sure we start with referenceIn
            convertInputsHere = deepcopy(convertInputs)
            convertInputsHere.remove(referenceIn)
            referenceIteration = True
            for ii in [referenceIn] + convertInputsHere:
                args[ii] = np.asarray(args[ii])
                inShape = args[ii].shape
                # special case we have to add dimensions
                if (
                  addDims or
                  (referenceIteration and (len(inShape) == noOfInDimsToKeep))
                  ):
                    inShapeNew = (1,) + inShape
                    addDims = True
                # we remove dimensions
                else:
                    # first loop
                    if referenceIteration:
                        noOfInDimsToKeepHere = noOfInDimsToKeep
                    # be more careful to remove not too many dimensions
                    else:
                        noOfInDimsToKeepHere = len(
                            inShape) - len(referenceShape)

                    # nothing to keep, just add
                    if noOfInDimsToKeepHere == 0:
                        inShapeKeep = ()
                    # find out which needs to be untouched
                    else:
                        inShapeKeep = inShape[-noOfInDimsToKeepHere:]
                    inShapeFlatten = inShape[:len(
                        inShape)-noOfInDimsToKeepHere]
                    # the reference to reshpe the output
                    if (ii == referenceIn):
                        referenceShape = deepcopy(inShapeFlatten)
                    inShapeNew = (
                        np.prod(inShapeFlatten).astype(int),) + inShapeKeep
                if verbosity > 0:
                    print('in', ii, inShape, inShapeNew,
                          referenceShape, addDims)
                args[ii] = args[ii].reshape(inShapeNew)
                referenceIteration = False
            # finally, run the function
            result = func(*args, **kwargs)
            # make sure output containes more than variable, otherwise make
            # it iterable
            if not isinstance(result, tuple):
                result = (result,)
            result = list(result)
            for oo in convertOutputs:
                result[oo] = np.asarray(result[oo])
                outShape = result[oo].shape
                # special case we had to add dimensions, now remove it
                if addDims:
                    outShapeNew = outShape[1:]
                else:
                    outShapeNew = referenceShape+outShape[1:]
                if verbosity > 0:
                    print('out', oo, outShape, outShapeNew, referenceShape)
                result[oo] = result[oo].reshape(outShapeNew)
            result = tuple(result)
            if len(result) == 1:
                result = result[0]
            return result
        return inner
    return NDto2DtoND_decorator
