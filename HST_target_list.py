from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u

'''
purpose: sort HST into target lists
first part of scripts goes through 3 criteria
'''

#HST data-- using the old i filter; this has already been dust corrected
HST_hdu = fits.open('/Users/amandaquirk/Documents/M33/data/m33-merged-F475W-F814W.fits', memmap=True)
HST_data = HST_hdu[1].data 
HST_RA = HST_data['ra'] #degrees
HST_Dec = HST_data['dec'] #degrees
F814W = HST_data['f814w_vega']
F475W = HST_data['f475w_vega']

'''
step 1) see if stars are isolated or not =================================================================
'''
HST_coords = SkyCoord(HST_RA, HST_Dec, unit=(u.deg, u.deg))

isolated = np.ones((len(HST_RA))) #will contain a 0 if the star doesn't pass the isolation criteria and a 1 if it does
for i in range(len(HST_coords)):
	star1 = HST_coords[i] #go through star by star
	I_tgt = F814W[i]
	sep = star1.separation(HST_coords) #get distances from that star to all other stars in region
	close = (sep.arcsecond < 10) & (sep.arcsecond > 0) #narrow our neighbor search for closer stars to save time; don't want to include the target star
	close_coords = HST_coords[close]
	close_mag = F814W[close] 
	close_distances = sep[close]
	reject = I_tgt > close_mag + (close_distances.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = close_coords[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		isolated[i] = 1
print('done with isolation!')

'''
step 2) see if stars are in the masks and if so, which mask =================================================================
'''
mask1_center = SkyCoord('1h34m02.7303s', '+30d44m11.000s')
mask2_center = SkyCoord('1h34m07.0646s', '+30d48m07.179s')
mask3_center = SkyCoord('1h33m52.6102s', '+30d32m15.893s')
mask4_center = SkyCoord('1h33m57.4856s', '+30d40m12.844s')
mask5_center = SkyCoord('1h33m55.2257s', '+30d36m14.442s')

half_length = 984 / 2 / 3600 #degrees
half_height = 240 / 2 / 3600 #degrees

RA_min1 = mask1_center.ra.degree - half_length
RA_max1 = mask1_center.ra.degree + half_length
Dec_min1 = mask1_center.dec.degree - half_height
Dec_max1 = mask1_center.dec.degree + half_height

RA_min2 = mask2_center.ra.degree - half_length
RA_max2 = mask2_center.ra.degree + half_length
Dec_min2 = mask2_center.dec.degree - half_height
Dec_max2 = mask2_center.dec.degree + half_height

RA_min3 = mask3_center.ra.degree - half_length
RA_max3 = mask3_center.ra.degree + half_length
Dec_min3 = mask3_center.dec.degree - half_height
Dec_max3 = mask3_center.dec.degree + half_height

RA_min4 = mask4_center.ra.degree - half_length
RA_max4 = mask4_center.ra.degree + half_length
Dec_min4 = mask4_center.dec.degree - half_height
Dec_max4 = mask4_center.dec.degree + half_height

RA_min5 = mask5_center.ra.degree - half_length
RA_max5 = mask5_center.ra.degree + half_length
Dec_min5 = mask5_center.dec.degree - half_height
Dec_max5 = mask5_center.dec.degree + half_height

mask = np.zeros(len(HST_RA)) #will contain the mask number (as noted above) that the star is in, 0 if outside of masks
for i in range(len(HST_RA)):
	if (HST_RA[i] > RA_min1) & (HST_RA[i] < RA_max1) & (HST_Dec[i] > Dec_min1) & (HST_Dec[i] < Dec_max1):
		mask[i] = 1
	elif (HST_RA[i] > RA_min2) & (HST_RA[i] < RA_max2) & (HST_Dec[i] > Dec_min2) & (HST_Dec[i] < Dec_max2):
		mask[i] = 2
	elif (HST_RA[i] > RA_min3) & (HST_RA[i] < RA_max3) & (HST_Dec[i] > Dec_min3) & (HST_Dec[i] < Dec_max3):
		mask[i] = 3
	elif (HST_RA[i] > RA_min4) & (HST_RA[i] < RA_max4) & (HST_Dec[i] > Dec_min4) & (HST_Dec[i] < Dec_max4):
		mask[i] = 4
	elif (HST_RA[i] > RA_min5) & (HST_RA[i] < RA_max5) & (HST_Dec[i] > Dec_min5) & (HST_Dec[i] < Dec_max5):
		mask[i] = 5
print('done with mask assignment!')

'''
step 3) see if stars are bright
'''
mag_cut = np.zeros(len(HST_RA)) #will be 0 if the star doesn't pass the magnitude cut and 1 if it does
for i in range(len(HST_RA)):
	if (F814W[i] < 22) | (F475W[i] < 24):
		mag_cut[i] = 1

print('done with brightness check!')

'''
save the data =================================================================
'''
np.savetxt('/Users/amandaquirk/Desktop/HST_target_criteria.txt', np.c_[HST_RA, HST_Dec, F814W, F475W, isolated, mask, mag_cut], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), F814W (mag), F475W (mag), isolation criteria, mask, brightness criteria')