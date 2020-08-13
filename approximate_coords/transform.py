from numpy import sin, cos, arctan2, arcsin, pi

from .time import local_hour_angle, earth_rotation_angle


def radec_to_altaz(ra, dec, time, location):
    '''
    Convert right ascension / declination to altitude / azimuth.
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


def altaz_to_radec(alt, az, time, location):
    '''
    Convert altitude / azimuth to right ascension / declination.
    Azimuth is measured from North (0°) through East (90°)
    '''

    # declare variables to match notation
    ϕ = location.lat.rad
    λ = location.lon.rad

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

    dec = arcsin(sin_a * sinϕ + cos_a * cos_Az * cosϕ)
    ra = (earth_rotation_angle(time) + λ - h) % (2 * pi)

    return ra, dec
