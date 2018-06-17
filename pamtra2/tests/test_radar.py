import collections
import numpy as np
import pytest
import pamtra2
import pamtra2.libs.refractiveIndex as refractiveIndex
import xarray as xr


@pytest.fixture()
def create_simple_cloud_creator():
    def simple_cloud_creator(
        nHeights=1,
        edr=1e-3,
        Ntot=[1],
        scattering=pamtra2.hydrometeors.scattering.Rayleigh,
        radarPNoise1000=-30,
        instrument='simple',
        dask=False,
        additionalDims=collections.OrderedDict(),
        size=0.001
    ):

        additionalDims = additionalDims

        pam2 = pamtra2.pamtra2(
            nLayer=nHeights,
            hydrometeors=['cloud'],
            additionalDims=additionalDims,
            frequencies=[3e9],
        )

        pam2.profile.height[:] = 1000
        pam2.profile.temperature[:] = 250
        pam2.profile.relativeHumidity[:] = 90
        pam2.profile.pressure[:] = 100000
        pam2.profile.eddyDissipationRate[:] = edr
        pam2.profile.horizontalWind[:] = 10
        pam2.profile['heightBinDepth'] = 10
        pam2.addMissingVariables()

        pam2.describeHydrometeor(
            pamtra2.hydrometeors.softEllipsoidFixedDensity,
            name='cloud',  # or None, then str(index)
            nBins=2,
            sizeBounds=pamtra2.hydrometeors.size.linspaceBounds,
            sizeCenter=pamtra2.hydrometeors.size.boundsToMid,
            sizeBoundsWidth=pamtra2.hydrometeors.size.boundsWidth,
            sizeDistribution=pamtra2.hydrometeors.sizeDistribution.\
            monoDisperse,
            aspectRatio=1.0,
            mass=pamtra2.hydrometeors.mass.ellipsoid,
            density=pamtra2.hydrometeors.density.water,
            crossSectionArea=pamtra2.hydrometeors.crossSectionArea.sphere,
            relativePermittivity=refractiveIndex.water.turner_kneifel_cadeddu,
            scattering=scattering,
            fallVelocity=pamtra2.hydrometeors.fallVelocity.\
            khvorostyanov01_drops,
            Dmin=size - .5e-10,
            Dmax=size + .5e-10,
            Ntot=xr.DataArray(Ntot, coords=[pam2.profile.layer]),
            useFuncArgDefaults=False,
        )

        pam2.profile['pathIntegratedAtenuattion'] = xr.zeros_like(
            pam2.hydrometeors.cloud.profile.backscatterCrossSection.isel(
                sizeBin=0))

        if dask:
            pam2.profile = pam2.profile.chunk({'layer': 1, 'frequency': 1})

        if instrument == 'simple':
            results = pam2.addInstrument(
                pamtra2.instruments.radar.simpleRadar,
                name='simple',
                frequencies=3e9,
            )
        elif instrument == 'spectral':
            results = pam2.addInstrument(
                pamtra2.instruments.radar.dopplerRadarPamtra,
                name='spectral',
                frequencies=3e9,
                momentsNPeaks=1,
                seed=11,
                radarAliasingNyquistInterv=0,
                radarPNoise1000=radarPNoise1000
            )

        if dask:
            results.results.load()

        return results
    return simple_cloud_creator


def test_rayleigh_mie(create_simple_cloud_creator):
    ray = create_simple_cloud_creator(
        scattering=pamtra2.hydrometeors.scattering.Rayleigh
    ).results.radarReflectivity.values
    mie = create_simple_cloud_creator(
        scattering=pamtra2.hydrometeors.scattering.Mie
    ).results.radarReflectivity.values

    assert np.allclose(ray, mie, rtol=1e-01, atol=1e-01)


def test_rayleigh_scale_N(create_simple_cloud_creator):
    ray = create_simple_cloud_creator(
        Ntot=[0.1, 1, 10],
        nHeights=3,
    ).results.radarReflectivity.values.flatten()

    assert np.allclose(ray[0], -10, rtol=1e-01, atol=1e-01)
    assert np.allclose(ray[1], 0, rtol=1e-01, atol=1e-01)
    assert np.allclose(ray[2], 10, rtol=1e-01, atol=1e-01)


def test_rayleigh_scale_N_dask(create_simple_cloud_creator):
    ray = create_simple_cloud_creator(
        Ntot=[0.1, 1, 10],
        nHeights=3,
        instrument='spectral',
        dask=True,
    ).results.radarReflectivity.values.flatten()

    assert np.allclose(ray[0], -10, rtol=1e-01, atol=2e-01)
    assert np.allclose(ray[1], 0, rtol=1e-01, atol=2e-01)
    assert np.allclose(ray[2], 10, rtol=1e-01, atol=2e-01)


def test_rayleigh_scale_size(create_simple_cloud_creator):
    ray_1mm = create_simple_cloud_creator(
        size=0.001,
    ).results.radarReflectivity.values.flatten()
    ray_1cm = create_simple_cloud_creator(
        size=0.01,
    ).results.radarReflectivity.values.flatten()
    ray_1um = create_simple_cloud_creator(
        size=0.0001,
    ).results.radarReflectivity.values.flatten()

    assert np.allclose(ray_1um, -60, rtol=1e-01, atol=1e-01)
    assert np.allclose(ray_1mm, 0, rtol=1e-01, atol=1e-01)
    assert np.allclose(ray_1cm, 60, rtol=1e-01, atol=1e-01)


def test_simple_spectra(create_simple_cloud_creator):
    simple = create_simple_cloud_creator(
    ).results.radarReflectivity.values.flatten()
    spectral = create_simple_cloud_creator(
        instrument='spectral',
    ).results.radarReflectivity.values.flatten()

    assert np.allclose(simple, spectral, rtol=1e-01, atol=1e-01)


def test_dask(create_simple_cloud_creator):
    dask = create_simple_cloud_creator(
        instrument='spectral',
        # additionalDims={'lat': np.arange(10)},
        dask=True,
        nHeights=100,
        Ntot=[100]*100,
    )
    nodask = create_simple_cloud_creator(
        instrument='spectral',
        # additionalDims={'lat': np.arange(10)},
        nHeights=100,
        Ntot=[100]*100,
    )
    # assert 0
    assert np.allclose(dask.results.radarReflectivity.values.flatten(),
                       nodask.results.radarReflectivity.values.flatten(),
                       rtol=1e-01, atol=2e-01)


@pytest.mark.parametrize("noise", [
    -30, -20, -10,
])
def test_scale_noise(create_simple_cloud_creator, noise):
    m10 = create_simple_cloud_creator(
        instrument='spectral',
        radarPNoise1000=noise,
    ).results.noiseMean.values.flatten()
    assert np.allclose(m10, noise, rtol=1e-01, atol=1e-01)


@pytest.mark.parametrize(("edr", 'result'), [
    (1e-3, 0.4), (1e-4, 0.4/2.1), (1e-5, 0.4/2.1**2),
])
def test_scale_edr(create_simple_cloud_creator, edr, result):
    m10 = create_simple_cloud_creator(
        instrument='spectral',
        edr=edr,
    ).results.spectrumWidth.values.flatten()
    assert np.allclose(m10, result, rtol=2e-01, atol=2e-01)


# def test_attenuation2pia():
#     arr = xr.DataArray(np.ones(4), coords={'layer': range(4)}, dims=['layer'])
#     PIA_bottomup, PIA_topdown = pamtra2.instruments.radar._attenuation2pia(arr)
#     assert np.all(PIA_bottomup.values == PIA_topdown.values[::-1])
#     assert np.all(PIA_bottomup.values == np.array([1., 3., 5., 7.]))
