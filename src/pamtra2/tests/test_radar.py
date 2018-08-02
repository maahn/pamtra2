import collections

import numpy as np
import pamtra2
import pamtra2.libs.refractiveIndex as refractiveIndex
import pytest
import xarray as xr


@pytest.fixture()
def create_simple_cloud_creator():
    def simple_cloud_creator(
        nHeights=1,
        edr=1e-3,
        Ntot=[1],
        scattering=pamtra2.hydrometeors.scattering.Rayleigh,
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        water_turner_kneifel_cadeddu,
        radarPNoise1000=-30,
        temperature=250,
        instrument='simple',
        dask=False,
        additionalDims=collections.OrderedDict(),
        size=0.001,
        hydrometeor='cloud',
        hydrometeorContent=0.0001,
        **kwargs
    ):

        additionalDims = additionalDims

        pam2 = pamtra2.pamtra2(
            nLayer=nHeights,
            hydrometeors=['hydrometeor'],
            additionalDims=additionalDims,
            frequencies=[3e9],
        )

        pam2.profile.height[:] = 1000
        pam2.profile.temperature[:] = temperature
        pam2.profile.relativeHumidity[:] = 90
        pam2.profile.pressure[:] = 100000
        pam2.profile.eddyDissipationRate[:] = edr
        pam2.profile.horizontalWind[:] = 10
        pam2.profile.hydrometeorContent[:] = hydrometeorContent

        pam2.profile['heightBinDepth'] = 10
        pam2.addMissingVariables()

        if hydrometeor == 'cloud':
            pam2.describeHydrometeor(
                pamtra2.hydrometeors.softEllipsoidFixedDensity,
                name='hydrometeor',  # or None, then str(index)
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
                relativePermittivity=relativePermittivity,
                scattering=scattering,
                fallVelocity=pamtra2.hydrometeors.fallVelocity.\
                khvorostyanov01_drops,
                Dmin=size - .5e-10,
                Dmax=size + .5e-10,
                Ntot=xr.DataArray(Ntot, coords=[pam2.profile.layer]),
                checkTemperatureForRelativePermittivity=False,
                useFuncArgDefaults=False,
                **kwargs,
            )
        elif hydrometeor == 'snow':
            pam2.describeHydrometeor(
                pamtra2.hydrometeors.softEllipsoidMassSize,
                name='hydrometeor',  # or None, then str(index)
                nBins=10,
                sizeBounds=pamtra2.hydrometeors.size.linspaceBounds,
                sizeCenter=pamtra2.hydrometeors.size.boundsToMid,
                sizeBoundsWidth=pamtra2.hydrometeors.size.boundsWidth,
                sizeDistribution=pamtra2.hydrometeors.sizeDistribution.
                exponentialFieldWC,
                aspectRatio=1.0,
                mass=pamtra2.hydrometeors.mass.powerLaw,
                density=pamtra2.hydrometeors.density.softEllipsoid,
                crossSectionArea=pamtra2.hydrometeors.crossSectionArea.sphere,
                relativePermittivity=relativePermittivity,
                scattering=scattering,
                fallVelocity=pamtra2.hydrometeors.fallVelocity.
                heymsfield10_particles,
                Dmin=size,
                Dmax=size + 0.01,
                Ntot=xr.DataArray(Ntot, coords=[pam2.profile.layer]),
                checkTemperatureForRelativePermittivity=False,
                useFuncArgDefaults=False,
                massSizeA=0.0121,
                massSizeB=1.9,
                minDensity=100,
                maxDensity=pamtra2.constants.rhoIce,
                **kwargs,
            )

        pam2.profile['pathIntegratedAtenuattion'] = xr.zeros_like(
            pam2.hydrometeors.hydrometeor.profile.backscatterCrossSection.isel(
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


def test_refractiveIndex_liquid(create_simple_cloud_creator):
    turner_kneifel_cadeddu = create_simple_cloud_creator(
        nHeights=2,
        Ntot=[1, 2],
        temperature=273.15
    ).results.radarReflectivity.values
    ellison = create_simple_cloud_creator(
        nHeights=2,
        Ntot=[1, 2],
        temperature=273.15,
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        water_ellison
    ).results.radarReflectivity.values
    assert np.allclose(turner_kneifel_cadeddu, ellison, rtol=1e-01, atol=1e-01)


def test_refractiveIndex_ice(create_simple_cloud_creator):
    ice_warren_brandt_2008 = create_simple_cloud_creator(
        nHeights=1,
        Ntot=[1],
        temperature=260.15,
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        ice_warren_brandt_2008
    ).results.radarReflectivity.values
    ice_matzler_2006 = create_simple_cloud_creator(
        nHeights=1,
        Ntot=[1],
        temperature=260.15,
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        ice_matzler_2006
    ).results.radarReflectivity.values
    ice_iwabuchi_yang_2011 = create_simple_cloud_creator(
        nHeights=1,
        Ntot=[1],
        temperature=[260.15],
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        ice_iwabuchi_yang_2011
    ).results.radarReflectivity.values
    assert np.allclose(ice_warren_brandt_2008,
                       ice_matzler_2006, rtol=1e-01, atol=1e-01)
    assert np.allclose(ice_iwabuchi_yang_2011,
                       ice_matzler_2006, rtol=1e-01, atol=1e-01)


def test_refractiveIndex_mixing(create_simple_cloud_creator):

    mixing_maxwell_garnett = create_simple_cloud_creator(
        nHeights=2,
        Ntot=[1, 2],
        temperature=273.15,
        hydrometeor='snow',
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        mixing_maxwell_garnett,
        relativePermittivityIce=pamtra2.hydrometeors.relativePermittivity.
        ice_iwabuchi_yang_2011,

    ).results.radarReflectivity.values

    mixing_bruggeman = create_simple_cloud_creator(
        nHeights=2,
        Ntot=[1, 2],
        temperature=273.15,
        hydrometeor='snow',
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        mixing_bruggeman,
        relativePermittivityIce=pamtra2.hydrometeors.relativePermittivity.
        ice_iwabuchi_yang_2011,

    ).results.radarReflectivity.values

    mixing_sihvola = create_simple_cloud_creator(
        nHeights=2,
        Ntot=[1, 2],
        temperature=273.15,
        hydrometeor='snow',
        relativePermittivity=pamtra2.hydrometeors.relativePermittivity.
        mixing_sihvola,
        relativePermittivityIce=pamtra2.hydrometeors.relativePermittivity.
        ice_iwabuchi_yang_2011,

    ).results.radarReflectivity.values
    assert np.allclose(mixing_maxwell_garnett,
                       mixing_bruggeman, rtol=1e-01, atol=1e-01)
    assert np.allclose(mixing_maxwell_garnett,
                       mixing_sihvola, rtol=1e-01, atol=1e-01)


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
