from astropy.coordinates import SkyCoord, AltAz, CIRS, ICRS, EarthLocation
from astropy.time import Time

crab = SkyCoord.from_name('Crab Nebula')
time = Time('2020-01-01T23:00')
location = EarthLocation.of_site('Roque de los Muchachos')


def test_trafo():
    from approximate_coords.change_trafo import use_transformation
    from approximate_coords.almanac_chp_6 import altaz_to_cirs_only_precession

    # this uses normal astropy conversion
    crab_altaz = crab.transform_to(AltAz(obstime=time, location=location))

    with use_transformation(AltAz, CIRS, altaz_to_cirs_only_precession):
        crab_approx = crab_altaz.transform_to(ICRS)

    assert crab_approx.separation(crab).arcsecond < 30
