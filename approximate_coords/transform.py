from numpy import sin, cos, arctan2, arcsin, pi
from astropy.coordinates import (
    FunctionTransform, UnitSphericalRepresentation, SphericalRepresentation,
    AltAz, CIRS
)
import astropy.units as u

from .time import local_hour_angle, earth_rotation_angle


def is_unit_repr(coord):
    return (
        isinstance(coord.data, UnitSphericalRepresentation)
        or coord.cartesian.x.unit == u.one
    )


def cirs_to_altaz(ra, dec, time, location):
    '''
    Convert intermediate right ascension / declination to altitude / azimuth.
    Azimuth is measured from North (0°) through East (90°)
    '''

    # declare variables to match notation
    δ = dec
    ϕ = location.lat.rad
    h = local_hour_angle(ra, time, location)

    # precalculate terms we need more than once
    sinδ = sin(δ)
    cosδ = cos(δ)
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)
    sin_h = sin(h)
    cos_h = cos(h)

    # divide line1 by line 2 of (7.16), p. 256
    az = arctan2(
        - sinδ * sin_h,
        sinδ * cosϕ - cosδ * cos_h * sinϕ
    )

    alt = arcsin(sinδ * sinϕ + cosδ * cos_h * cosϕ)

    return alt, az


def altaz_to_cirs(alt, az, time, location):
    '''
    Convert altitude / azimuth to intermediate right ascension / declination.
    Azimuth is measured from North (0°) through East (90°)
    '''

    # declare variables to match notation
    λ, ϕ, _ = location.to_geodetic('WGS84')
    λ = λ.to_value(u.rad)
    ϕ = ϕ.to_value(u.rad)

    # precalculate terms we need more than once
    sin_a = sin(alt)
    cos_a = cos(alt)
    sin_Az = sin(az)
    cos_Az = cos(az)
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)

    # divide line1 by line 2 of (7.19), p. 257
    h = arctan2(
        - cos_a * sin_Az,
        sin_a * cosϕ - cos_a * cos_Az * sinϕ
    )

    dec_i = arcsin(sin_a * sinϕ + cos_a * cos_Az * cosϕ)
    ra_i = (earth_rotation_angle(time) + λ - h) % (2 * pi)

    return ra_i, dec_i


def _altaz_to_cirs(fromcoord, toframe):
    ra_i, dec_i = altaz_to_cirs(
        alt=fromcoord.alt.rad,
        az=fromcoord.az.rad,
        time=fromcoord.obstime,
        location=fromcoord.location,
    )

    ra_i = u.Quantity(ra_i, u.rad, copy=False)
    dec_i = u.Quantity(dec_i, u.rad, copy=False)

    if is_unit_repr(fromcoord):
        rep = UnitSphericalRepresentation(lon=ra_i, lat=dec_i)
    else:
        rep = SphericalRepresentation(
            lon=ra_i, lat=dec_i, distance=fromcoord.distance
        )

    return toframe.realize_frame(rep)


def _cirs_to_altaz(fromcoord, toframe):
    alt, az = cirs_to_altaz(
        ra=fromcoord.ra.rad,
        dec=fromcoord.dec.rad,
        time=fromcoord.obstime,
        location=fromcoord.location,
    )

    alt = u.Quantity(alt, u.rad, copy=False)
    az = u.Quantity(az, u.rad, copy=False)

    if is_unit_repr(fromcoord):
        rep = UnitSphericalRepresentation(lon=az, lat=alt)
    else:
        rep = SphericalRepresentation(
            lon=az, lat=alt, distance=fromcoord.distance
        )

    return toframe.realize_frame(rep)


transform_altaz_to_cirs = FunctionTransform(_altaz_to_cirs, AltAz, CIRS)
transform_cirs_to_altaz = FunctionTransform(_cirs_to_altaz, CIRS, AltAz)
