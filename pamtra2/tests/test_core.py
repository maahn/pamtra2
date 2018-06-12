import collections
import numpy as np
import pytest
import pandas as pn
import pamtra2
import refractiveIndex


@pytest.fixture
def initilize_pamtra():
    additionalDims = collections.OrderedDict()
    additionalDims['time'] = pn.date_range(
        '2016-01-01', '2016-01-05', freq='D')[:1]
    additionalDims['lat'] = np.arange(70, 80)
    nHeights = 100

    pam2 = pamtra2.pamtra2(
        nLayer=nHeights,
        hydrometeors=['rain', 'snow'],
        additionalDims=additionalDims,
        frequencies=[35e9, 94e9],
    )
    return pam2


@pytest.fixture
def add_data(initilize_pamtra):
    pam2 = initilize_pamtra
    pam2.profile.height[:] = np.linspace(0, 1000, pam2.nLayer)
    pam2.profile.temperature[:] = 250
    pam2.profile.relativeHumidity[:] = 90
    pam2.profile.pressure[:] = 10000
    pam2.profile.eddyDissipationRate[:] = 1e-4
    pam2.profile.horizontalWind[:] = 0

    pam2.profile.waterContent.values[:] = 0
    # rain
    pam2.profile.waterContent.values[..., 20:40, 0] = 1e-4
    # snow
    pam2.profile.waterContent.values[..., 20:40, 1] = 2e-4

    return pam2


@pytest.fixture
def add_rain(add_data):
    pam2 = add_data
    pam2.describeHydrometeor(
        pamtra2.hydrometeors.softEllipsoidFixedDensity,
        name='rain',  # or None, then str(index)
        kind='liquid',  # liquid, ice
        nBins=40,
        sizeCenter=pamtra2.hydrometeors.sizeCenter.linspace,
        sizeDistribution=pamtra2.hydrometeors.sizeDistribution.exponentialN0WC,
        aspectRatio=1.0,
        mass=pamtra2.hydrometeors.mass.ellipsoid,
        density=pamtra2.hydrometeors.density.water,
        crossSectionArea=pamtra2.hydrometeors.crossSectionArea.sphere,
        # replace with refractiveIndex.water.Turner.n
        relativePermittivity=refractiveIndex.water.turner_kneifel_cadeddu,
        scattering=pamtra2.hydrometeors.scattering.Rayleigh,
        Dmin=1e-6,
        Dmax=1e-2,
        N0=8e6,
        model='Turner',
        useFuncArgDefaults=True,
    )


pam2.profile.waterContent.values[:] = 0
# rain
pam2.profile.waterContent.values[..., 20:40, 0] = 1e-4
# snow
pam2.profile.waterContent.values[..., 20:40, 1] = 2e-4


class InsufficientAmount(Exception):
    pass


class Wallet(object):

    def __init__(self, initial_amount=0):
        self.balance = initial_amount

    def spend_cash(self, amount):
        if self.balance < amount:
            raise InsufficientAmount(
                'Not enough available to spend {}'.format(amount))
        self.balance -= amount

    def add_cash(self, amount):
        self.balance += amount


@pytest.fixture
def empty_wallet():
    '''Returns a Wallet instance with a zero balance'''
    return Wallet()


@pytest.fixture
def wallet(empty_wallet):
    '''Returns a Wallet instance with a balance of 20'''
    empty_wallet.add_cash(20)
    return empty_wallet


def test_default_initial_amount(empty_wallet):
    assert empty_wallet.balance == 0


def test_setting_initial_amount(wallet):
    assert wallet.balance == 20


def test_wallet_add_cash(wallet):
    wallet.add_cash(80)
    assert wallet.balance == 100


def test_wallet_spend_cash(wallet):
    wallet.spend_cash(10)
    assert wallet.balance == 10


def test_wallet_spend_cash_raises_exception_on_insufficient_amount(empty_wallet):
    with pytest.raises(InsufficientAmount):
        empty_wallet.spend_cash(100)
