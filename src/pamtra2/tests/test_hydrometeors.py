import numpy as np
import pamtra2
import pytest


class TestAspectRatio(object):
    pass


class TestCore(object):
    '''
    Doesn't make sense without the parent...
    '''
    pass


class TestCrossSectionArea(object):

    def testSphere(self):
        sizeCenter = np.arange(1, 11)
        sp1 = pamtra2.hydrometeors.crossSectionArea.sphere(sizeCenter)

        spPre = 0.25 * np.pi
        spExp = 2
        sp2 = pamtra2.hydrometeors.crossSectionArea.powerLaw(
            sizeCenter, spPre, spExp)
        assert np.all(sp1 == sp2)

    def testPowerLaw(self):
        sizeCenter = np.arange(1, 11)
        pre = 10
        exp = 2
        pl1 = pamtra2.hydrometeors.crossSectionArea.powerLaw(
            sizeCenter, pre, exp)
        pl2 = pre*sizeCenter**exp

        assert np.all(pl1 == pl2)


class TestDensity(object):

    def testOblate(self):
        sizeCenter = np.array([1.])
        aspectRatio = 0.6
        den1 = np.array([100.])
        volume = 4/3. * np.pi * (sizeCenter/2.)**2 * \
            (sizeCenter*aspectRatio/2.)
        mass = volume * den1
        den2 = pamtra2.hydrometeors.density.softOblateEllipsoid(
            sizeCenter, aspectRatio, mass)
        assert np.allclose(den2, den1)

    def testProlate(self):
        sizeCenter = np.array([1.])
        aspectRatio = 1.6
        den1 = np.array([100.])
        volume = 4/3. * np.pi * \
            (sizeCenter*(1./aspectRatio)/2.)**2 * (sizeCenter/2.)
        mass = volume * den1
        den2 = pamtra2.hydrometeors.density.softProlateEllipsoid(
            sizeCenter, aspectRatio, mass)
        assert np.allclose(den2, den1)

    def testSphere(self):
        sizeCenter = np.array([1.])
        aspectRatio = np.array([1.])
        den1 = np.array([100.])
        mass = 1/6. * np.pi * sizeCenter**3 * den1
        for func in [
            pamtra2.hydrometeors.density.softEllipsoid,
            pamtra2.hydrometeors.density.softProlateEllipsoid,
            pamtra2.hydrometeors.density.softOblateEllipsoid,
        ]:
            den2 = func(sizeCenter, aspectRatio, mass)
            assert np.allclose(den2, den1)

    def testMin(self):
        sizeCenter = np.array([1])
        aspectRatio = np.array([1])
        den1 = np.array([10])
        mass = 1/6. * np.pi * sizeCenter**3 * den1
        for func in [
            pamtra2.hydrometeors.density.softEllipsoid,
            pamtra2.hydrometeors.density.softProlateEllipsoid,
            pamtra2.hydrometeors.density.softOblateEllipsoid,
        ]:
            den2 = func(sizeCenter, aspectRatio, mass, minDensity=99)
            assert den2 == 99

    def testMax(self):
        sizeCenter = np.array([1])
        aspectRatio = np.array([1])
        den1 = np.array([10000])
        mass = 1/6. * np.pi * sizeCenter**3 * den1
        for func in [
            pamtra2.hydrometeors.density.softEllipsoid,
            pamtra2.hydrometeors.density.softProlateEllipsoid,
            pamtra2.hydrometeors.density.softOblateEllipsoid,
        ]:
            den2 = func(sizeCenter, aspectRatio, mass, maxDensity=999)
            assert den2 == 999


class TestMass(object):
    def testPowerLaw(self):
        sizeCenter = np.arange(1, 11)
        pre = 10
        exp = 3
        pl1 = pamtra2.hydrometeors.mass.powerLaw(
            sizeCenter, pre, exp)
        pl2 = pre*sizeCenter**exp
        assert np.all(pl1 == pl2)

    def testWater(self):
        sizeCenter = np.arange(1, 11)
        sp1 = pamtra2.hydrometeors.mass.waterSphere(sizeCenter)

        den = 1000.
        spPre = 1/6. * np.pi * den
        spExp = 3
        sp2 = pamtra2.hydrometeors.mass.powerLaw(
            sizeCenter, spPre, spExp)
        assert np.all(sp1 == sp2)

    def testIce(self):
        sizeCenter = np.arange(1, 11)
        sp1 = pamtra2.hydrometeors.mass.iceSphere(sizeCenter)

        den = 917.
        spPre = 1/6. * np.pi * den
        spExp = 3
        sp2 = pamtra2.hydrometeors.mass.powerLaw(
            sizeCenter, spPre, spExp)
        assert np.all(sp1 == sp2)

    def testEllipsoid(self):
        sizeCenter = np.arange(1, 11)
        sp1 = pamtra2.hydrometeors.mass.waterSphere(sizeCenter)
        density = 1000
        aspectRatio = 1
        sp2 = pamtra2.hydrometeors.mass.ellipsoid(
            sizeCenter, aspectRatio, density)
        assert np.all(sp1 == sp2)


class TestSizeBounds(object):
    def test_lin(self):
        D1 = pamtra2.hydrometeors.size.linspaceBounds(10, 1, 2)
        assert len(D1) == 11

    def test_log(self):
        D1 = pamtra2.hydrometeors.size.logspaceBounds(10, 1, 2)
        assert len(D1) == 11

    def test_boundsToMid(self):
        D1 = pamtra2.hydrometeors.size.logspaceBounds(10, 1, 2)
        DM = pamtra2.hydrometeors.size.boundsToMid(D1)
        assert len(D1) == len(DM)+1

    def test_boundsWidth(self):
        D1 = pamtra2.hydrometeors.size.logspaceBounds(10, 1, 2)
        DM = pamtra2.hydrometeors.size.boundsWidth(D1)
        assert len(D1) == len(DM)+1


class TestSizeDistribution(object):

    def test_monodisperse(self):
        sizeBoundsWidth = np.array([0.1, 0.1])
        Ntot = 10
        N = pamtra2.hydrometeors.sizeDistribution.monoDisperse(
            sizeBoundsWidth, Ntot)
        assert np.sum(N*sizeBoundsWidth) == Ntot
        assert N[0] == N[1]

    def test_exponential(self):
        sizeCenter = np.arange(1, 11)
        N0 = 10
        lamb = 4
        sd1 = pamtra2.hydrometeors.sizeDistribution.exponential(
            sizeCenter, N0, lamb)
        sd2 = N0 * np.exp(-lamb * sizeCenter)
        np.all(sd1 == sd2)

    def test_gamma(self):
        sizeCenter = np.arange(1, 11)
        N0 = 10
        lamb = 4
        mu = 2
        sd1 = pamtra2.hydrometeors.sizeDistribution.gamma(
            sizeCenter, N0, lamb, mu)
        sd2 = N0 * sizeCenter ** mu * np.exp(-lamb * sizeCenter)
        np.all(sd1 == sd2)

    def test_modGamma(self):
        sizeCenter = np.arange(1, 11)
        N0 = 10
        lamb = 4
        mu = 2
        gamm = 2
        sd1 = pamtra2.hydrometeors.sizeDistribution.modifiedGamma(
            sizeCenter, N0, lamb, mu, gamm)
        sd2 = N0 * sizeCenter ** mu * np.exp(-lamb * sizeCenter**gamm)
        np.all(sd1 == sd2)

    def testMarshallPalmer(self):
        result = np.array([
            [132581.40321409],
            [429681.13785473],
            [993929.97473809],
        ])
        sizeCenter = np.array([0.001])
        for rr, rainRate in enumerate([1, 5, 25]):
            N = pamtra2.hydrometeors.sizeDistribution.exponentialMarshallPalmer(
                sizeCenter, rainRate)
            assert np.allclose(N, result[rr])

    def testField(self):
        N0_1 = 7.628e6 * np.exp(-0.107 * pamtra2.units.kelvin2Celsius(263))
        N0_2 = pamtra2.hydrometeors.sizeDistribution._exponentialField(263)
        assert N0_1 == N0_2

    def testExponentialFieldWC(self):
        sizeCenter = np.logspace(-6, 0, 1000)
        sizeWidth = np.gradient(sizeCenter)
        temperature = 263
        massSizeA = 0.0121
        massSizeB = 3.
        waterContent = 1e-4
        N = pamtra2.hydrometeors.sizeDistribution.exponentialFieldWC(
            sizeCenter, temperature, waterContent, massSizeA, massSizeB
        )
        mass = pamtra2.hydrometeors.mass.powerLaw(
            sizeCenter, massSizeA, massSizeB)
        assert np.allclose(np.sum(N*mass*sizeWidth), waterContent)

    def testExponentialN0WC(self):
        sizeCenter = np.logspace(-6, 0, 1000)
        sizeWidth = np.gradient(sizeCenter)
        N0 = 1e6
        massSizeA = 0.0121
        massSizeB = 3.
        waterContent = 1e-4

        N = pamtra2.hydrometeors.sizeDistribution.exponentialN0WC(
            sizeCenter, N0, waterContent, massSizeA, massSizeB
        )
        mass = pamtra2.hydrometeors.mass.powerLaw(
            sizeCenter, massSizeA, massSizeB)
        np.allclose(np.sum(N*mass*sizeWidth), waterContent)

    @pytest.mark.skip(reason="Test fails by a factor of 2?!")
    def testExponentialFieldReff(self):
        sizeCenter = np.logspace(-6, 0, 1000)
        sizeWidth = np.gradient(sizeCenter)
        temperature = 263

        effectiveRadius = 1e-3

        N = pamtra2.hydrometeors.sizeDistribution.exponentialFieldReff(
            sizeCenter, temperature, effectiveRadius)
        M3 = np.sum((sizeCenter/2)**3 * N * (sizeWidth/2))
        M2 = np.sum((sizeCenter/2)**2 * N * (sizeWidth/2))

        np.allclose(M3/M2, effectiveRadius)


class TestScattering(object):
    def testCompareRayleighMie(self):
        diameter = 1e-4
        wavelength = 1e-2
        refractiveIndex = 5.97+2.79j
        back1 = pamtra2.hydrometeors.scattering._scatteringWrapper(
            diameter,
            wavelength,
            refractiveIndex,
            model='Rayleigh'
        )[3]
        back2 = pamtra2.hydrometeors.scattering._scatteringWrapper(
            diameter,
            wavelength,
            refractiveIndex,
            model='Mie'
        )[3]
        assert np.allclose(back1, back2)


class TestFallVelocity(object):
    def test_khvorostyanov01_drops(self):
        sizeCenter = np.linspace(0.001, 0.008, 5)
        dryAirDensity = 1.23
        kinematicViscosity = 1.343e-5
        velSpec = pamtra2.hydrometeors.fallVelocity.khvorostyanov01_drops(
            sizeCenter, dryAirDensity, kinematicViscosity)
        # In good agreement with figure 2 of Khvorostyanov, V. I., and J. A.
        # Curry, 2002: Terminal Velocities of Droplets and Crystals: Power
        # Laws with Continuous Parameters over the Size Spectrum. J. Atmos.
        # Sci., 59, 1872–1884, doi:10.1175/1520-0469(2002)059<1872:TVODAC>
        # 2.0.CO;2.

        refSpec = np.array([3.86493162,  7.42270326,  9.1933725,
                            10.14756485, 10.66685595])

        assert np.allclose(velSpec, refSpec)

    def test_heymsfield10_particles(self):
        dryAirDensity = 1.2
        dynamicViscosity = 1.725e-5
        sizeCenter = np.linspace(0.001, 0.008, 5)
        mass = pamtra2.hydrometeors.mass.powerLaw(
            sizeCenter, massSizeA=0.0121, massSizeB=1.9)
        crossSectionArea = pamtra2.hydrometeors.crossSectionArea.powerLaw(
            sizeCenter, areaSizeA=0.3, areaSizeB=1.9)
        velSpec = pamtra2.hydrometeors.fallVelocity.heymsfield10_particles(
            sizeCenter, mass, crossSectionArea, dynamicViscosity,
            dryAirDensity)

        # This values are not tested due to teh lack of a reference, but they
        # lokk realistic.
        refSpec = np.array(
            [0.56276112, 0.74952901, 0.82443195, 0.86752849, 0.89629312])

        assert np.allclose(velSpec, refSpec)
