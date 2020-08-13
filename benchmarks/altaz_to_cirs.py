from astropy.coordinates import SkyCoord, AltAz, CIRS, EarthLocation
from astropy.time import Time
import astropy.units as u
import numpy as np

from approximate_coords.altaz_cirs import no_polar_motion

import warnings
from astropy.utils.exceptions import ErfaWarning, AstropyWarning
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore', category=ErfaWarning)
warnings.filterwarnings('ignore', category=AstropyWarning)


coord = SkyCoord.from_name('Crab Nebula')
time = Time('2000-01-01T00:01') + u.Quantity(np.arange(0, 50 * 365, 7), u.day)
location = EarthLocation.of_site('Roque de los Muchachos')


coord_altaz = coord.transform_to(AltAz(obstime=time, location=location))

# use standard astropy coordinate transform as reference
coord_cirs = coord_altaz.transform_to(CIRS(obstime=time))

# swap in the approximate transform
with no_polar_motion():
    coord_approx = coord_altaz.transform_to(CIRS(obstime=time))


separation = coord_approx.separation(coord_cirs).arcsecond
dra = (coord_approx.ra - coord_cirs.ra).wrap_at('180d').arcsecond
ddec = (coord_approx.dec - coord_cirs.dec).wrap_at('180d').arcsecond

axs = [
    ['time', 'time', 'time'],
    ['sep', 'dra', 'ddec'],
]

fig = plt.figure(constrained_layout=True)
axd = fig.subplot_mosaic(axs)

axd['time'].plot(time.to_datetime(), dra, label=r'$\Delta$ra')
axd['time'].plot(time.to_datetime(), ddec, label=r'$\Delta$dec')
axd['time'].plot(time.to_datetime(), separation, '-', label='sep', color='k')
axd['time'].legend()
axd['time'].set_ylabel('difference / as')


axd['sep'].hist(separation, bins=100, color='k')
axd['sep'].set_xlabel('separation / as')
axd['dra'].hist(dra, bins=100, color='C0')
axd['dra'].set_xlabel('$\Delta$ra / as')
axd['ddec'].hist(ddec, bins=100, color='C1')
axd['ddec'].set_xlabel('$\Delta$dec / as')

fig.savefig('build/altaz_to_cirs.pdf')
