from collections import OrderedDict

import numpy as np
import pamtra2
import xarray as xr


class TestConcatDicts(object):
    def setUp(self):
        self.d1 = OrderedDict({'a': 1, 'b': 2})
        self.d2 = OrderedDict({'c': 3})

    def test_order(self):
        self.setUp()
        dm = pamtra2.helpers.concatDicts(self.d2, self.d1)
        assert list(dm.keys())[0] == 'c'
        assert list(dm.values())[0] == 3

    def test_len(self):
        self.setUp()
        dm = pamtra2.helpers.concatDicts(self.d2, self.d1)
        assert len(dm) == 3
        assert list(dm.values())[0] == 3


class TestSwapListItems(object):
    def test_len(self):
        l1 = [1, 2, 3, 4]
        l2 = pamtra2.helpers.swapListItems(l1, 2, 3)
        assert len(l2) == len(l1)

    def test_sum(self):
        l1 = [1, 2, 3, 4]
        l2 = pamtra2.helpers.swapListItems(l1, 2, 3)
        assert np.sum(l2) == np.sum(l1)

    def test_swap(self):
        l1 = [1, 2, 3, 4]
        l2 = pamtra2.helpers.swapListItems(l1, 2, 3)
        assert l2[1] == l1[2]
        assert l2[2] == l1[1]


class TestAttrDict(object):
    def test_equal(self):
        d1 = {'a': 1, 'b': 2}
        d2 = pamtra2.helpers.AttrDict(d1)
        assert np.all(np.array(list(d1.values())) ==
                      np.array(list(d2.values())))
        assert np.all(np.array(list(d1.keys())) == np.array(list(d2.keys())))

    def test_attr(self):
        d1 = {'a': 1, 'b': 2}
        d2 = pamtra2.helpers.AttrDict(d1)
        assert hasattr(d2, 'a')
        assert hasattr(d2, 'b')


class TestDimensionToVariables(object):
    def test_sum(self):
        self.arr = xr.DataArray(np.arange(6).reshape((2, 3)), dims=['x', 'y'])
        keys = ['a', 'b', 'c']
        arr2 = pamtra2.helpers.dimensionToVariables(self.arr, 'y', keys)
        assert (arr2.a.sum() + arr2.b.sum() + arr2.c.sum()) == self.arr.sum()


class TestGetInputCoreDims(object):
    def test_coreDims(self):
        arr1 = xr.DataArray(np.arange(6).reshape((2, 3)), dims=['x', 'y'])
        arr2 = xr.DataArray(np.arange(3), dims=['y'])
        self.ds = xr.Dataset({'a': arr1, 'b': arr2})
        args = [self.ds.a, self.ds.b]
        coreDims = ['x']
        res = pamtra2.helpers.getInputCoreDims(args, coreDims)
        assert res == [['x'], []]
