import astropy.units as u
from astropy.time import Time
import numpy as np


J2000 = Time('J2000.0', scale='tt')

DAYS_IN_CENTURY = 36525
TWO_PI = 2 * np.pi


def delta_julian_century(time, epoch=J2000):
    '''
    Convert an astropy time object to julian centuries from the J2000 epoch

    This is denoted with `T` in the expl. suppl.
    '''
    return modified_julian_date(time, epoch) / DAYS_IN_CENTURY


def modified_julian_date(time, epoch=J2000):
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

    DU = modified_julian_date(time, epoch=Time('J2000.0', scale='ut1'))
    era = TWO_PI * (0.7790_5727_32640 + 1.0027_3781_1911_35448 * DU)
    return era % TWO_PI


def local_hour_angle(ra, time, location):
    '''
    Calculate local hour angle. right column of (7.6), p. 255

    This is denoted with `h` in the expl. suppl.
    '''
    λ = location.lon.rad
    return earth_rotation_angle(time) + λ - ra
