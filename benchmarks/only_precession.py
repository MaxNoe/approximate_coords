from astropy.coordinates import SkyCoord, AltAz, ICRS, CIRS, EarthLocation
from astropy.time import Time
import astropy.units as u
import numpy as np

from approximate_coords.change_trafo import use_transformation
from approximate_coords.almanac_chp_6 import (
    cirs_to_icrs_only_precession,
    altaz_to_icrs_only_precession,
    altaz_to_radec_only_precession,
    altaz_to_cirs,
)
import warnings
from astropy.utils.exceptions import ErfaWarning
import matplotlib.pyplot as plt
from time import sleep

warnings.filterwarnings('ignore', category=ErfaWarning)


coord = SkyCoord.from_name('Mrk 501')
time = Time('2000-01-01T00:01') + u.Quantity(np.arange(0, 50 * 365, 7), u.day)
location = EarthLocation.of_site('Roque de los Muchachos')


# this uses normal astropy conversion
coord_altaz = coord.transform_to(AltAz(obstime=time, location=location))

ra, dec = altaz_to_cirs(alt=coord_altaz.alt.rad, az=coord_altaz.az.rad, time=time, location=location)
coord_approx = SkyCoord(ra=ra, dec=dec, frame='cirs', obstime=time)

coord = coord_altaz.transform_to(CIRS(obstime=time))
# coord_approx = coord_altaz.transform_to(ICRS)

# with use_transformation(CIRS, ICRS, cirs_to_icrs_only_precession):
#     coord_approx = coord_altaz.transform_to(ICRS)

separation = coord_approx.separation(coord).arcsecond
dra = (coord_approx.ra - coord.ra).wrap_at('180d').arcsecond
ddec = (coord_approx.dec - coord.dec).wrap_at('180d').arcsecond

axs = [
    ['time', 'time', 'time'],
    ['sep', 'dra', 'ddec'],
]

fig = plt.figure(constrained_layout=True)
axd = fig.subplot_mosaic(axs)

axd['time'].plot(time.to_datetime(), separation, '-', label='sep')
axd['time'].plot(time.to_datetime(), dra, label=r'$\Delta$ra')
axd['time'].plot(time.to_datetime(), ddec, label=r'$\Delta$dec')
axd['time'].legend()


axd['sep'].hist(separation, bins=100)
axd['dra'].hist(dra, bins=100)
axd['ddec'].hist(ddec, bins=100)

fig.savefig('build/only_precession.pdf')
