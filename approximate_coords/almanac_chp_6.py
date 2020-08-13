'''
This is an implementation of the Formulas in
The Explanatory Supplement to the Astronomical Almanac,
Sean E. Urban and P. Kenneth Seidelmann,
3rd Edition,
Section 6.6.2.4, pp 220ff

Notation here follows the formulas in the book exactly.

It is claimed to be accurate to 30" between 1900 and 2100.
Which is supported by the benchmarks that compare to astropy.

Not that the all formulas are missing the M quantity, the precession
offset in right ascension, as this is already included in the definition
of the Earth Rotation Angle, see p 78 for the explanation.
'''
from contextlib import contextmanager
import numpy as np
from astropy.coordinates import (
    FunctionTransform,
    SphericalRepresentation,
    ICRS,
    CIRS,
    UnitSphericalRepresentation,
)
import astropy.units as u

from .time import delta_julian_century
from .change_trafo import use_transformation


def calc_N(T):
    '''
    Calculate precession angle N to a precision of 1", line 2 of (6.34), p. 221
    '''
    return np.deg2rad(
        0.556_72 * T
        - 0.000_12 * T**2
        - 0.000_01 * T**3
    )


def mean_to_j2000(alpha, delta, N):
    '''
    Calculate mean right ascension and declination for conversion from
    coordinates using the mean equinox and equator of date to J2000,
    line 1 and 2 of (6.33)
    '''

    alpha_m = alpha - 0.5 * (N * np.sin(alpha) * np.tan(delta))
    delta_m = delta - 0.5 * N * np.cos(alpha_m)

    return alpha_m, delta_m


def to_j2000(alpha, delta, T):
    '''
    Apply precession correction for conversion from
    coordinates using the mean equinox and equator of date to j2000
    according to section 6.6.2.4, formulas
    at the top of page 221.

    The M part is not necessary, as it is included in the earth rotation angle,
    so this correction is already applied to a coordinate in CIRS.
    '''
    N = calc_N(T)

    alpha_m, delta_m = mean_to_j2000(alpha, delta, N)

    alpha0 = alpha - N * np.sin(alpha_m) * np.tan(delta_m)
    delta0 = delta - N * np.cos(alpha_m)

    return alpha0, delta0


def mean_from_j2000(alpha0, delta0, N):
    '''
    Calculate mean right ascension and declination for conversion from
    J2000 to coordinates using the mean equinox and equator of date to J2000,
    line 3 and 4 of (6.33)
    '''
    alpha_m = alpha0 + 0.5 * N * np.sin(alpha0) * np.tan(delta0)
    delta_m = delta0 + 0.5 * N * np.cos(alpha_m)

    return alpha_m, delta_m


def from_j2000(alpha0, delta0, T):
    '''
    Apply precession correction for conversion from J2000
    to coordinates using the mean equinox and equator of date
    according to section 6.6.2.4,
    formulas at the top of page 221.
    '''
    N = calc_N(T)

    alpha_m, delta_m = mean_from_j2000(alpha0, delta0, N)
    alpha = alpha0 + N * np.sin(alpha_m) * np.tan(delta_m)
    delta = delta0 + N * np.cos(alpha_m)

    return alpha, delta


def _cirs_to_icrs_only_precession(fromcoord, toframe):
    T = delta_julian_century(fromcoord.obstime)

    ra, dec = fromcoord.ra.rad, fromcoord.dec.rad
    ra, dec = to_j2000(ra, dec, T)

    ra = u.Quantity(ra, u.rad, copy=False)
    dec = u.Quantity(np.clip(dec, -np.pi/2, np.pi/2), u.rad, copy=False)

    if isinstance(fromcoord.data, UnitSphericalRepresentation):
        rep = UnitSphericalRepresentation(lon=ra, lat=dec)
    else:
        rep = SphericalRepresentation(
            lon=ra, lat=dec, distance=fromcoord.distance
        )

    return toframe.realize_frame(rep)


def _icrs_to_cirs_only_precession(fromcoord, toframe):
    T = delta_julian_century(fromcoord.obstime)

    ra, dec = fromcoord.ra.rad, fromcoord.dec.rad
    ra, dec = from_j2000(ra, dec, T)

    ra = u.Quantity(ra, u.rad, copy=False)
    dec = u.Quantity(np.clip(dec, -np.pi/2, np.pi/2), u.rad, copy=False)

    if isinstance(fromcoord.data, UnitSphericalRepresentation):
        rep = UnitSphericalRepresentation(lon=ra, lat=dec)
    else:
        rep = SphericalRepresentation(
            lon=ra, lat=dec, distance=fromcoord.distance
        )

    return toframe.realize_frame(rep)


cirs_to_icrs_approximate_precession = FunctionTransform(
    func=_cirs_to_icrs_only_precession,
    fromsys=CIRS,
    tosys=ICRS,
)

icrs_to_cirs_approximate_precession = FunctionTransform(
    func=_cirs_to_icrs_only_precession,
    fromsys=CIRS,
    tosys=ICRS,
)


@contextmanager
def approximate_precesssion():
    '''
    Context manager to use the transformations between CIRS and ICRS
    that only use an approximated precession and no nutation information.

    Accurate to about 30 arcseconds.
    '''

    try:
        with use_transformation(ICRS, CIRS, icrs_to_cirs_approximate_precession):
            with use_transformation(CIRS, ICRS, cirs_to_icrs_approximate_precession):
                yield
    finally:
        pass
