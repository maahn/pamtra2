import pamtra2


def test_kelvin():
    assert 0 == pamtra2.units.kelvin2Celsius(273.15)


def test_kelvin2():
    assert -273.15 == pamtra2.units.kelvin2Celsius(0)
