import pygasabs
import numpy as np

'''
Values are NOT compared with any reference...
'''


def test_pygasabs_rosenkranz98():

    absair, abswv = pygasabs.calculate_gas_absorption_rosenkranz98(
        100e9, 300, 138.45, 1000*100)

    assert np.allclose(absair, 6.0573921e-06)
    assert np.allclose(abswv, 9.7943904e-06)


def test_pygasabs_vec_rosenkranz98():

    absair, abswv = pygasabs.calculate_gas_absorption_rosenkranz98(
        np.array([100e9, 200e9]),
        np.array([300, 200]),
        np.array([138.45, 0]),
        np.array([1000*100, 900*100]))

    assert np.allclose(absair[0], 6.0573921e-06)
    assert np.allclose(abswv[0], 9.7943904e-06)
    assert np.allclose(absair[1], 1.08597465e-05)
    assert np.allclose(abswv[1], 0)


def test_pygasabs_liebe93():

    absair = pygasabs.calculate_gas_absorption_liebe93(
        100e9, 300, 138.45, 1000*100)

    assert np.allclose(absair, 1.52462334e-05)


def test_pygasabs_vec_liebe93():

    absair = pygasabs.calculate_gas_absorption_liebe93(
        np.array([100e9, 200e9]),
        np.array([300, 200]),
        np.array([138.45, 0]),
        np.array([1000*100, 900*100]))

    assert np.allclose(absair, np.array([1.52462334e-05, 1.10830357e-05]))
