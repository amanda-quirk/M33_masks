from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u
from shapely.geometry import Polygon
from shapely.geometry import Point 

'''
purpose: sort CFHT into target lists
first part of scripts goes through 3 criteria
'''

#CFHT data-- using the old i filter; this has already been dust corrected
CFHT_hdu = fits.open('/Users/amandaquirk/Documents/M33/data/M33.GI.matchcat.extcorr.fits')
CFHT_data = CFHT_hdu[1].data
CFHT_RA = CFHT_data['RA2'] #degrees
CFHT_Dec = CFHT_data['DEC2'] #degrees
CFHT_i = CFHT_data['MAG2_AUTO_0']
CFHT_g = CFHT_data['MAG1_AUTO_0']

#use below just to look for bugs in code
# idx = [0, 1, 2, 3, 4, 5]
# CFHT_RA = CFHT_RA[idx]
# CFHT_Dec = CFHT_Dec[idx]
# CFHT_i = CFHT_i[idx]
# CFHT_g = CFHT_g[idx]

'''
step 1) see if stars are isolated or not =================================================================
-reject a star if it has at least one neighbor that satisfies: I_neighbor < I_target - (d/0.8")^1.5 +2
'''
CFHT_coords = SkyCoord(CFHT_RA, CFHT_Dec, unit=(u.deg, u.deg))
ras = CFHT_coords.ra.degree
decs = CFHT_coords.dec.degree

indices_dec = np.arange(0, len(CFHT_coords)).astype(np.int)
indices_ra = np.arange(0, len(CFHT_coords)).astype(np.int)

decs_sorted, indices_sorted_dec = np.array(sorted(zip(decs,indices_dec))).T
ras_sorted, indices_sorted_ra = np.array(sorted(zip(ras,indices_ra))).T

indices_sorted_dec = indices_sorted_dec.astype(np.int)
indices_sorted_ra = indices_sorted_ra.astype(np.int)

window = 0.00277777777

isolated = np.zeros(len(CFHT_RA)) #will contain a 0 if the star doesn't pass the isolation criteria and a 1 if it does
for i in range(int(len(CFHT_coords)/1000)):
	star1 = CFHT_coords[i] #go through star by star
	I_tgt = CFHT_i[i]
	ra = star1.ra.degree
	dec = star1.dec.degree
	left_ra = ra - window
	left_dec = dec - window 
	right_ra = ra + window
	right_dec = dec + window
	left_dec_index = np.searchsorted(decs_sorted, left_dec)
	right_dec_index = np.searchsorted(decs_sorted, right_dec) - 1
	left_ra_index = np.searchsorted(ras_sorted, left_ra) 
	right_ra_index = np.searchsorted(ras_sorted, right_ra) - 1
	in_window_decs = indices_sorted_dec[left_dec_index: right_dec_index]
	in_window_ras = set(indices_sorted_ra[left_ra_index: right_ra_index])
	in_window_index = list(in_window_ras.intersection(in_window_decs))
	CFHT_coords_in_window = CFHT_coords[in_window_index]
	CFHT_i_in_window = CFHT_i[in_window_index]
	sep = star1.separation(CFHT_coords_in_window) #get distances from that star to all other stars in region
	reject = I_tgt > CFHT_i_in_window + (sep.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = CFHT_coords_in_window[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		isolated[i] = 1

# plt.scatter(CFHT_RA, CFHT_Dec, c=isolated)
# plt.xlim(23, 24)
# plt.ylim(30, 31)
# plt.colorbar()
# plt.savefig('/Users/amandaquirk/Desktop/checking_isolation_original.png')
# plt.close()
# print('done with isolation!')


