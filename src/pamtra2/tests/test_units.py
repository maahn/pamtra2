import pamtra2


def test_kelvin():
    assert 0 == pamtra2.units.kelvin2Celsius(273.15)


def test_kelvin2():
    assert -273.15 == pamtra2.units.kelvin2Celsius(0)


def test_celsius():
    assert 0 == pamtra2.units.celsius2Kelvin(-273.15)


def test_celsius2():
    assert 273.15 == pamtra2.units.celsius2Kelvin(0)
