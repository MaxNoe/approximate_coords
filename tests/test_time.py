from astropy.time import Time
import numpy as np
import erfa


def test_julian_century():
    from approximate_coords.time import delta_julian_century

    assert delta_julian_century(Time('J2000.0', scale='tt')) == 0
    assert delta_julian_century(Time('2100-01-01T12:00:00', scale='tt')) == 1


def test_era():
    from approximate_coords.time import earth_rotation_angle
    from astropy.coordinates.builtin_frames.utils import get_jd12

    time = Time(['2010-01-01', '2020-01-01'])
    assert np.allclose(earth_rotation_angle(time), erfa.era00(*get_jd12(time, scale='ut1')))
