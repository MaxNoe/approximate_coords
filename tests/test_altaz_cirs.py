from astropy.coordinates.angle_utilities import angular_separation
from astropy.coordinates import SkyCoord, AltAz, CIRS, EarthLocation
import numpy as np


def test_altaz_to_cirs():
    from approximate_coords.transform import altaz_to_cirs

    alts = [-45, 0, 25]
    azs = np.arange(0, 360, 22.5)

    alt = np.deg2rad(np.repeat(alts, len(azs)))
    az = np.deg2rad(np.tile(azs, len(alts)))

    frame = AltAz(
        obstime='2010-01-01T20:00',
        location=EarthLocation.of_site('Roque de los Muchachos'),
    )
    coord_altaz = SkyCoord(alt=alt, az=az, unit='rad', frame=frame)
    coord_cirs = coord_altaz.transform_to(CIRS)


    ra_i, dec_i = altaz_to_cirs(alt, az, frame.obstime, frame.location)

    sep = 3600 * np.rad2deg(angular_separation(ra_i, dec_i, coord_cirs.ra.rad, coord_cirs.dec.rad))
    assert np.all(sep < 40)
