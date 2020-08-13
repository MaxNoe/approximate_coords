import astropy.units as u
from astropy.time import Time
import numpy as np


J2000_TT = Time('J2000.0', scale='tt')
J2000_UT1 = Time('J2000.0', scale='ut1')

DAYS_IN_CENTURY = 36525
TWO_PI = 2 * np.pi


def delta_julian_century(time, epoch=J2000_TT):
    '''
    Convert an astropy time object to julian centuries from the J2000 epoch

    This is denoted with `T` in the expl. suppl.
    '''
    return modified_julian_date(time, epoch) / DAYS_IN_CENTURY


def modified_julian_date(time, epoch=J2000_TT):
    '''
    Get the modified julian date (number of days since the specified epoch)

    This is denoted with `D_U` in the expl. suppl.
    '''
    return (time - epoch).to(u.day).value


def earth_rotation_angle(time):
    '''
    Precise formula for the ERA from (7.157), p. 299

    This is denoted with `θ` in the expl. suppl.
    '''
    # we need the difference in UT1 days, which are not
    # 86400 SI Seconds but depend on the orientation of Earth as measured
    # by radio telescopes.
    # UTC would also be ok as it is < 1 second from UT1 at all times.
    # But astropy does not support UTC timedeltas.
    DU = modified_julian_date(time.ut1, epoch=J2000_UT1)

    era = TWO_PI * (0.7790_5727_32640 + 1.0027_3781_1911_35448 * DU)
    return era % TWO_PI


def local_hour_angle(ra, time, location):
    '''
    Calculate local hour angle. right column of (7.6), p. 255

    This is denoted with `h` in the expl. suppl.
    '''
    λ = location.lon.rad
    return earth_rotation_angle(time) + λ - ra
