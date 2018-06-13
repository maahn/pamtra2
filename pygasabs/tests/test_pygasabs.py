import pygasabs
import numpy as np

'''
Values are NOT compared with any reference...
'''


def test_pygasabs():

    absair, abswv = pygasabs.calculate_gas_absorption(
        100e9, 300, 0.001, 1000*100)

    assert np.allclose(absair, 6.0573921e-06)
    assert np.allclose(abswv, 9.7943904e-06)


def test_pygasabs_vec():

    absair, abswv = pygasabs.calculate_gas_absorption(
        np.array([100e9, 200e9]),
        np.array([300, 200]),
        np.array([0.001, 0]),
        np.array([1000*100, 900*100]))

    assert np.allclose(absair[0], 6.0573921e-06)
    assert np.allclose(abswv[0], 9.7943904e-06)
    assert np.allclose(absair[1], 1.08597465e-05)
    assert np.allclose(abswv[1], 0)
