# -*- coding: utf-8 -*-
import numpy as np


def linspace(nBins, Dmin=None, Dmax=None):

    return np.linspace(Dmin, Dmax, nBins)


def logspace(nBins, Dmin=None, Dmax=None):

    return np.logspace(np.log10(Dmin), np.log10(Dmax), nBins)


# from .core2 import *
