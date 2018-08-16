from astropy.io import fits
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u

#CFHT data-- using the old i filter; this has already been dust corrected
CFHT_hdu = fits.open('/Users/amandaquirk/Documents/M33/data/M33.GI.matchcat.extcorr.fits')
CFHT_data = CFHT_hdu[1].data
CFHT_RA = CFHT_data['RA2'] #degrees
CFHT_Dec = CFHT_data['DEC2'] #degrees
CFHT_i = CFHT_data['MAG2_AUTO_0']

CFHT_coords = SkyCoord(CFHT_RA, CFHT_Dec, unit=(u.deg, u.deg))
star1 = CFHT_coords[0]
star2 = CFHT_coords[1]

'''
test #1: astropy
'''
astropy_sep = star1.separation(star2) 
astropy_sep_unit = astropy_sep.arcsecond
print(astropy_sep_unit)

'''
test #2: formula
'''

def wiki_formula(point1, point2):
	ra1 = star1.ra.radian
	dec1 = star1.dec.radian
	ra2 = star2.ra.radian
	dec2 = star2.dec.radian
	angle = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(abs(ra1 - ra2))
	return np.arccos(angle) * 3600 * 180 / np.pi 

print(wiki_formula(star1, star2))

np.savetxt('/Users/amandaquirk/Desktop/CFHT_data.txt', np.c_[CFHT_RA, CFHT_Dec, CFHT_i], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), i (mag)')