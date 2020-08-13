'''
This is an implementation of the Formulas in
The Explanatory Supplement to the Astronomical Almanac,
Sean E. Urban and P. Kenneth Seidelmann,
3rd Edition,
Section 6.6.2.4, pp 220ff

Notation here follows the formulas in the book exactly.

It is said to be accurate to 30" between 1900 and 2100.
'''

import numpy as np
from astropy.coordinates import (
    FunctionTransform,
    SphericalRepresentation,
    AltAz,
    ICRS,
    CIRS,
    UnitSphericalRepresentation,
)
import astropy.units as u

from .time import delta_julian_century
from .transform import altaz_to_radec


def calc_M(T):
    '''
    Calculate precession angle M to a precision of 1", line 1 of (6.34), p. 221
    '''

    return np.deg2rad(
        1.281_16 * T
        + 0.000_39 * T**2
        + 0.000_01 * T**3
    )


def calc_N(T):
    '''
    Calculate precession angle N to a precision of 1", line 2 of (6.34), p. 221
    '''
    return np.deg2rad(
        0.556_72 * T
        - 0.000_12 * T**2
        - 0.000_01 * T**3
    )


def mean_to_j2000(alpha, delta, M, N):
    '''
    Calculate mean right ascension and declination for conversion from
    coordinates using the mean equinox and equator of date to J2000,
    line 1 and 2 of (6.33)
    '''

    alpha_m = alpha - 0.5 * (M + N * np.sin(alpha) * np.tan(delta))
    delta_m = delta - 0.5 * N * np.cos(alpha_m)

    return alpha_m, delta_m


def to_j2000(alpha, delta, T):
    '''
    Apply precession correction for conversion from
    coordinates using the mean equinox and equator of date to j2000
    according to section 6.6.2.4, formulas
    at the top of page 221.
    '''
    M = calc_M(T)
    N = calc_N(T)

    alpha_m, delta_m = mean_to_j2000(alpha, delta, M, N)

    alpha0 = alpha - M - N * np.sin(alpha_m) * np.tan(delta_m)
    delta0 = delta - N * np.cos(alpha_m)

    return alpha0, delta0


def mean_from_j2000(alpha0, delta0, M, N):
    '''
    Calculate mean right ascension and declination for conversion from
    J2000 to coordinates using the mean equinox and equator of date to J2000,
    line 3 and 4 of (6.33)
    '''
    alpha_m = alpha0 + 0.5 * (M + N * np.sin(alpha0) * np.tan(delta0))
    delta_m = delta0 + 0.5 * N * np.cos(alpha_m)

    return alpha_m, delta_m


def from_j2000(alpha0, delta0, T):
    '''
    Apply precession correction for conversion from J2000
    to coordinates using the mean equinox and equator of date
    according to section 6.6.2.4,
    formulas at the top of page 221.
    '''
    M = calc_M(T)
    N = calc_N(T)

    alpha_m, delta_m = mean_from_j2000(alpha0, delta0, M, N)
    alpha = alpha0 + M + N * np.sin(alpha_m) * np.tan(delta_m)
    delta = delta0 + N * np.cos(alpha_m)

    return alpha, delta


def altaz_to_radec_only_precession(alt, az, time, location):
    T = delta_julian_century(time)

    ra, dec = altaz_to_radec(
        alt=alt,
        az=az,
        time=time,
        location=location,
    )
    ra, dec = to_j2000(ra, dec, T)

    ra = u.Quantity(ra, u.rad, copy=False)
    dec = u.Quantity(dec, u.rad, copy=False)

    return ra, dec


def altaz_to_cirs(alt, az, time, location):
    T = delta_julian_century(time)

    ra, dec = altaz_to_radec(
        alt=alt,
        az=az,
        time=time,
        location=location,
    )

    ra = u.Quantity(ra, u.rad, copy=False)
    dec = u.Quantity(dec, u.rad, copy=False)

    return ra, dec


def _altaz_to_icrs_only_precession(fromcoord, toframe):
    T = delta_julian_century(fromcoord.obstime)

    ra, dec = altaz_to_radec(
        alt=fromcoord.alt.rad,
        az=fromcoord.az.rad,
        time=fromcoord.obstime,
        location=fromcoord.location,
    )
    ra, dec = to_j2000(ra, dec, T)

    ra = u.Quantity(ra, u.rad, copy=False)
    dec = u.Quantity(dec, u.rad, copy=False)

    if fromcoord.distance is None:
        rep = UnitSphericalRepresentation(lon=ra, lat=dec)
    else:
        rep = SphericalRepresentation(
            lon=ra, lat=dec, distance=fromcoord.distance
        )

    return toframe.realize_frame(rep)


def _cirs_to_icrs_only_precession(fromcoord, toframe):
    T = delta_julian_century(fromcoord.obstime)

    ra, dec = to_j2000(fromcoord.ra.rad, fromcoord.dec.rad, T)

    ra = u.Quantity(ra, u.rad, copy=False)
    dec = u.Quantity(np.clip(dec, -np.pi/2, np.pi/2), u.rad, copy=False)

    if isinstance(fromcoord.data, UnitSphericalRepresentation):
        rep = UnitSphericalRepresentation(lon=ra, lat=dec)
    else:
        rep = SphericalRepresentation(
            lon=ra, lat=dec, distance=fromcoord.distance
        )

    return toframe.realize_frame(rep)


altaz_to_icrs_only_precession = FunctionTransform(
    func=_altaz_to_icrs_only_precession,
    fromsys=AltAz,
    tosys=ICRS,
)

cirs_to_icrs_only_precession = FunctionTransform(
    func=_cirs_to_icrs_only_precession,
    fromsys=CIRS,
    tosys=ICRS,
)
