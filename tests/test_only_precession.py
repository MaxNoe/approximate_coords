from astropy.coordinates import SkyCoord, AltAz, ICRS, EarthLocation
from astropy.time import Time

crab = SkyCoord.from_name('Crab Nebula')
time = Time('2020-01-01T23:00')
location = EarthLocation.of_site('Roque de los Muchachos')


def test_trafo():
    from approximate_coords.almanac_chp_6 import approximate_precesssion

    # this uses normal astropy conversion
    crab_altaz = crab.transform_to(AltAz(obstime=time, location=location))

    with approximate_precesssion():
        crab_approx = crab_altaz.transform_to(ICRS)

    assert crab_approx.separation(crab).arcsecond < 30
