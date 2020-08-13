from astropy.coordinates import SkyCoord, AltAz, ICRS, EarthLocation
from astropy.time import Time
import astropy.units as u
import numpy as np

from approximate_coords.almanac_chp_6 import approximate_precesssion
import warnings
from astropy.utils.exceptions import ErfaWarning, AstropyWarning
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore', category=ErfaWarning)
warnings.filterwarnings('ignore', category=AstropyWarning)


coord = SkyCoord(ra='90d', dec='45d')
time = Time('2000-01-01T00:01') + u.Quantity(np.arange(0, 50 * 365, 7), u.day)
location = EarthLocation.of_site('Roque de los Muchachos')

coord_altaz = coord.transform_to(AltAz(obstime=time, location=location))


with approximate_precesssion():
    coord_approx = coord_altaz.transform_to(ICRS)

separation = coord_approx.separation(coord).arcsecond
dra = (coord_approx.ra - coord.ra).wrap_at('180d').arcsecond
ddec = (coord_approx.dec - coord.dec).wrap_at('180d').arcsecond


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
axd['dra'].set_xlabel(r'$\Delta$ra / as')
axd['ddec'].hist(ddec, bins=100, color='C1')
axd['ddec'].set_xlabel(r'$\Delta$dec / as')

fig.savefig('build/only_precession.pdf')
