from astropy.time import Time


def test_julian_century():
    from approximate_coords.time import delta_julian_century

    assert delta_julian_century(Time('J2000.0', scale='ut1')) == 0
    assert delta_julian_century(Time('2100-01-01T12:00:00', scale='ut1')) == 1
