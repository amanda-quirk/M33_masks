from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u

'''
purpose: examine the offset between the coordinate systems of HST and CFHT

'''

#CFHT data-- using the old i filter; this has already been dust corrected
CFHT_hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/M33.GI.matchcat.extcorr.fits')
CFHT_data = CFHT_hdu[1].data
CFHT_RA = CFHT_data['RA2'] #degrees
CFHT_Dec = CFHT_data['DEC2'] #degrees
CFHT_i = CFHT_data['MAG2_AUTO_0']
CFHT_g = CFHT_data['MAG1_AUTO_0']

#HST data
HST_hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/m33-merged-F475W-F814W.fits', memmap=True)
HST_data = HST_hdu[1].data 
HST_RA = HST_data['ra'] #degrees
HST_Dec = HST_data['dec'] #degrees
F814W = HST_data['f814w_vega']
F475W = HST_data['f475w_vega']

#limit data to not dense field
half_width = 150 #arcsec, half of the width of the box

RA_center1 = 23.5375 #deg, center of box
Dec_center1 = 30.5286 #deg, center of box
RA_min1 = RA_center1 - half_width/3600
RA_max1 = RA_center1 + half_width/3600
Dec_min1 = Dec_center1 - half_width/3600
Dec_max1 = Dec_center1 + half_width/3600

RA_center2 = 23.3792 #deg, center of box
Dec_center2 = 30.5536 #deg, center of box
RA_min2 = RA_center2 - half_width/3600
RA_max2 = RA_center2 + half_width/3600
Dec_min2 = Dec_center2 - half_width/3600
Dec_max2 = Dec_center2 + half_width/3600

RA_center3 = 23.4667 #deg, center of box
Dec_center3 = 30.8111 #deg, center of box
RA_min3 = RA_center3 - half_width/3600
RA_max3 = RA_center3 + half_width/3600
Dec_min3 = Dec_center3 - half_width/3600
Dec_max3 = Dec_center3 + half_width/3600

RA_center4 = 23.6167 #deg, center of box
Dec_center4 = 30.7892 #deg, center of box
RA_min4 = RA_center4 - half_width/3600
RA_max4 = RA_center4 + half_width/3600
Dec_min4 = Dec_center4 - half_width/3600
Dec_max4 = Dec_center4 + half_width/3600

#divide the stars into regions
HST_in_region1 = (HST_RA > RA_min1) & (HST_RA < RA_max1) & (HST_Dec > Dec_min1) & (HST_Dec < Dec_max1)      
CFHT_in_region1 = (CFHT_RA > RA_min1) & (CFHT_RA < RA_max1) & (CFHT_Dec > Dec_min1) & (CFHT_Dec < Dec_max1) 
HST_RA_region1 = HST_RA[HST_in_region1]
HST_Dec_region1 = HST_Dec[HST_in_region1]
F814W_region1 = F814W[HST_in_region1]
F475W_region1 = F475W[HST_in_region1]
CFHT_RA_region1 = CFHT_RA[CFHT_in_region1]
CFHT_Dec_region1 = CFHT_Dec[CFHT_in_region1]
CFHT_i_region1 = CFHT_i[CFHT_in_region1]
CFHT_g_region1 = CFHT_g[CFHT_in_region1]

HST_in_region2 = (HST_RA > RA_min2) & (HST_RA < RA_max2) & (HST_Dec > Dec_min2) & (HST_Dec < Dec_max2) 
CFHT_in_region2 = (CFHT_RA > RA_min2) & (CFHT_RA < RA_max2) & (CFHT_Dec > Dec_min2) & (CFHT_Dec < Dec_max2) 
HST_RA_region2 = HST_RA[HST_in_region2]
HST_Dec_region2 = HST_Dec[HST_in_region2]
F814W_region2 = F814W[HST_in_region2]
F475W_region2 = F475W[HST_in_region2]
CFHT_RA_region2 = CFHT_RA[CFHT_in_region2]
CFHT_Dec_region2 = CFHT_Dec[CFHT_in_region2]
CFHT_i_region2 = CFHT_i[CFHT_in_region2]
CFHT_g_region2 = CFHT_g[CFHT_in_region2]

HST_in_region3 = (HST_RA > RA_min3) & (HST_RA < RA_max3) & (HST_Dec > Dec_min3) & (HST_Dec < Dec_max3) 
CFHT_in_region3 = (CFHT_RA > RA_min3) & (CFHT_RA < RA_max3) & (CFHT_Dec > Dec_min3) & (CFHT_Dec < Dec_max3) 
HST_RA_region3 = HST_RA[HST_in_region3]
HST_Dec_region3 = HST_Dec[HST_in_region3]
F814W_region3 = F814W[HST_in_region3]
F475W_region3 = F475W[HST_in_region3]
CFHT_RA_region3 = CFHT_RA[CFHT_in_region3]
CFHT_Dec_region3 = CFHT_Dec[CFHT_in_region3]
CFHT_i_region3 = CFHT_i[CFHT_in_region3]
CFHT_g_region3 = CFHT_g[CFHT_in_region3]

HST_in_region4 = (HST_RA > RA_min4) & (HST_RA < RA_max4) & (HST_Dec > Dec_min4) & (HST_Dec < Dec_max4) 
CFHT_in_region4 = (CFHT_RA > RA_min4) & (CFHT_RA < RA_max4) & (CFHT_Dec > Dec_min4) & (CFHT_Dec < Dec_max4) 
HST_RA_region4 = HST_RA[HST_in_region4]
HST_Dec_region4 = HST_Dec[HST_in_region4]
F814W_region4 = F814W[HST_in_region4]
F475W_region4 = F475W[HST_in_region4]
CFHT_RA_region4 = CFHT_RA[CFHT_in_region4]
CFHT_Dec_region4 = CFHT_Dec[CFHT_in_region4]
CFHT_i_region4 = CFHT_i[CFHT_in_region4]
CFHT_g_region4 = CFHT_g[CFHT_in_region4]

'''
============================================================================================================================
-below preforms the isolation criteria for the HST data using F814W and the most relaxed critera
-reject a star if it has at least one neighbor that satisfies: I_neighbor < I_target - (d/0.8")^1.5 +2
'''

#covert the RA and Dec into astropy coordinates to calculate distance
HST_coords1 = SkyCoord(HST_RA_region1, HST_Dec_region1, unit=(u.deg, u.deg))
HST_coords2 = SkyCoord(HST_RA_region2, HST_Dec_region2, unit=(u.deg, u.deg))
HST_coords3 = SkyCoord(HST_RA_region3, HST_Dec_region3, unit=(u.deg, u.deg))
HST_coords4 = SkyCoord(HST_RA_region4, HST_Dec_region4, unit=(u.deg, u.deg))

#do isolation criteria in a way similar to the smoothing technique; maybe have it break once one neighbor is filled; should make this into a function later but i don't have time now
HST_RA_region1_keep = []
HST_Dec_region1_keep = []
F814W_region1_keep = []
F475W_region1_keep = []
for i in range(len(HST_coords1)):
	star1 = HST_coords1[i] #go through star by star
	I_tgt = F814W_region1[i]
	sep = star1.separation(HST_coords1) #get distances from that star to all other stars in region
	close = (sep.arcsecond < 10) & (sep.arcsecond > 0) #narrow our neighbor search for closer stars to save time; don't want to include the target star
	close_coords = HST_coords1[close]
	close_mag = F814W_region1[close] 
	close_distances = sep[close]
	reject = I_tgt > close_mag + (close_distances.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = close_coords[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		HST_RA_region1_keep.append(HST_RA_region1[i])
		HST_Dec_region1_keep.append(HST_Dec_region1[i])
		F475W_region1_keep.append(F475W_region1[i])
		F814W_region1_keep.append(F814W_region1[i])
print('Done with region 1 isolation')

HST_RA_region2_keep = []
HST_Dec_region2_keep = []
F814W_region2_keep = []
F475W_region2_keep = []
for i in range(len(HST_coords2)):
	star1 = HST_coords2[i] #go through star by star
	I_tgt = F814W_region2[i]
	sep = star1.separation(HST_coords2) #get distances from that star to all other stars in region
	close = (sep.arcsecond < 10) & (sep.arcsecond > 0) #narrow our neighbor search for closer stars to save time; don't want to include the target star
	close_coords = HST_coords2[close]
	close_mag = F814W_region2[close] 
	close_distances = sep[close]
	reject = I_tgt > close_mag + (close_distances.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = close_coords[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		HST_RA_region2_keep.append(HST_RA_region2[i])
		HST_Dec_region2_keep.append(HST_Dec_region2[i])
		F475W_region2_keep.append(F475W_region2[i])
		F814W_region2_keep.append(F814W_region2[i])
print('Done with region 2 isolation')

HST_RA_region3_keep = []
HST_Dec_region3_keep = []
F814W_region3_keep = []
F475W_region3_keep = []
for i in range(len(HST_coords3)):
	star1 = HST_coords3[i] #go through star by star
	I_tgt = F814W_region3[i]
	sep = star1.separation(HST_coords3) #get distances from that star to all other stars in region
	close = (sep.arcsecond < 10) & (sep.arcsecond > 0) #narrow our neighbor search for closer stars to save time; don't want to include the target star
	close_coords = HST_coords3[close]
	close_mag = F814W_region3[close] 
	close_distances = sep[close]
	reject = I_tgt > close_mag + (close_distances.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = close_coords[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		HST_RA_region3_keep.append(HST_RA_region3[i])
		HST_Dec_region3_keep.append(HST_Dec_region3[i])
		F475W_region3_keep.append(F475W_region3[i])
		F814W_region3_keep.append(F814W_region3[i])
print('Done with region 3 isolation')

HST_RA_region4_keep = []
HST_Dec_region4_keep = []
F814W_region4_keep = []
F475W_region4_keep = []
for i in range(len(HST_coords4)):
	star1 = HST_coords4[i] #go through star by star
	I_tgt = F814W_region4[i]
	sep = star1.separation(HST_coords4) #get distances from that star to all other stars in region
	close = (sep.arcsecond < 10) & (sep.arcsecond > 0) #narrow our neighbor search for closer stars to save time; don't want to include the target star
	close_coords = HST_coords4[close]
	close_mag = F814W_region4[close] 
	close_distances = sep[close]
	reject = I_tgt > close_mag + (close_distances.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = close_coords[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		HST_RA_region4_keep.append(HST_RA_region4[i])
		HST_Dec_region4_keep.append(HST_Dec_region4[i])
		F475W_region4_keep.append(F475W_region4[i])
		F814W_region4_keep.append(F814W_region4[i])
print('Done with region 4 isolation')

#make numpy arrays
HST_RA_region1_keep = np.array((HST_RA_region1_keep))
HST_Dec_region1_keep = np.array((HST_Dec_region1_keep))
F814W_region1_keep = np.array((F814W_region1_keep))
F475W_region1_keep = np.array((F475W_region1_keep))

HST_RA_region2_keep = np.array((HST_RA_region2_keep))
HST_Dec_region2_keep = np.array((HST_Dec_region2_keep))
F814W_region2_keep = np.array((F814W_region2_keep))
F475W_region2_keep = np.array((F475W_region2_keep))

HST_RA_region3_keep = np.array((HST_RA_region3_keep))
HST_Dec_region3_keep = np.array((HST_Dec_region3_keep))
F814W_region3_keep = np.array((F814W_region3_keep))
F475W_region3_keep = np.array((F475W_region3_keep))

HST_RA_region4_keep = np.array((HST_RA_region4_keep))
HST_Dec_region4_keep = np.array((HST_Dec_region4_keep))
F814W_region4_keep = np.array((F814W_region4_keep))
F475W_region4_keep = np.array((F475W_region4_keep))

'''
============================================================================================================================
only keep bright stars
'''

HST_bright1 = (F814W_region1_keep < 22) | (F475W_region1_keep < 24) 
CFHT_bright1 = (CFHT_i_region1 < 22) | (CFHT_g_region1 < 24)
HST_RA_region1 = HST_RA_region1_keep[HST_bright1]
HST_Dec_region1 = HST_Dec_region1_keep[HST_bright1]
F814W_region1 = F814W_region1_keep[HST_bright1]
CFHT_RA_region1 = CFHT_RA_region1[CFHT_bright1]
CFHT_Dec_region1 = CFHT_Dec_region1[CFHT_bright1]
CFHT_i_region1 = CFHT_i_region1[CFHT_bright1]

HST_bright2 = (F814W_region2_keep < 22) | (F475W_region2_keep < 24) 
CFHT_bright2 = (CFHT_i_region2 < 22) | (CFHT_g_region2 < 24)
HST_RA_region2 = HST_RA_region2_keep[HST_bright2]
HST_Dec_region2 = HST_Dec_region2_keep[HST_bright2]
F814W_region2 = F814W_region2_keep[HST_bright2]
CFHT_RA_region2 = CFHT_RA_region2[CFHT_bright2]
CFHT_Dec_region2 = CFHT_Dec_region2[CFHT_bright2]
CFHT_i_region2 = CFHT_i_region2[CFHT_bright2]

HST_bright3 = (F814W_region3_keep < 22) | (F475W_region3_keep < 24) 
CFHT_bright3 = (CFHT_i_region3 < 22) | (CFHT_g_region3 < 24)
HST_RA_region3 = HST_RA_region3_keep[HST_bright3]
HST_Dec_region3 = HST_Dec_region3_keep[HST_bright3]
F814W_region3 = F814W_region3_keep[HST_bright3]
CFHT_RA_region3 = CFHT_RA_region3[CFHT_bright3]
CFHT_Dec_region3 = CFHT_Dec_region3[CFHT_bright3]
CFHT_i_region3 = CFHT_i_region3[CFHT_bright3]

HST_bright4 = (F814W_region4_keep < 22) | (F475W_region4_keep < 24) 
CFHT_bright4 = (CFHT_i_region4 < 22) | (CFHT_g_region4 < 24)
HST_RA_region4 = HST_RA_region4_keep[HST_bright4]
HST_Dec_region4 = HST_Dec_region4_keep[HST_bright4]
F814W_region4 = F814W_region4_keep[HST_bright4]
CFHT_RA_region4 = CFHT_RA_region4[CFHT_bright4]
CFHT_Dec_region4 = CFHT_Dec_region4[CFHT_bright4]
CFHT_i_region4 = CFHT_i_region4[CFHT_bright4]

'''
============================================================================================================================
-do astropy coordinate matching
-only keep matches that are within +/-0.5 mag and are less than 1 arcsec in separation
'''

#convert to coordinates
HST_coords1 = SkyCoord(HST_RA_region1, HST_Dec_region1, unit=(u.deg, u.deg))
HST_coords2 = SkyCoord(HST_RA_region2, HST_Dec_region2, unit=(u.deg, u.deg))
HST_coords3 = SkyCoord(HST_RA_region3, HST_Dec_region3, unit=(u.deg, u.deg))
HST_coords4 = SkyCoord(HST_RA_region4, HST_Dec_region4, unit=(u.deg, u.deg))

CFHT_coords1 = SkyCoord(CFHT_RA_region1, CFHT_Dec_region1, unit=(u.deg, u.deg))
CFHT_coords2 = SkyCoord(CFHT_RA_region2, CFHT_Dec_region2, unit=(u.deg, u.deg))
CFHT_coords3 = SkyCoord(CFHT_RA_region3, CFHT_Dec_region3, unit=(u.deg, u.deg))
CFHT_coords4 = SkyCoord(CFHT_RA_region4, CFHT_Dec_region4, unit=(u.deg, u.deg))

#assign each CFHT data point the closest HST point
idx_1, d2d_1, d3d_1 = CFHT_coords1.match_to_catalog_sky(HST_coords1)
idx_2, d2d_2, d3d_2 = CFHT_coords2.match_to_catalog_sky(HST_coords2)
idx_3, d2d_3, d3d_3 = CFHT_coords3.match_to_catalog_sky(HST_coords3)
idx_4, d2d_4, d3d_4 = CFHT_coords4.match_to_catalog_sky(HST_coords4)

#calculate difference in the matched point's CFHT i' and HST F814W
mag_diff_region1 = abs(CFHT_i_region1 - F814W_region1[idx_1])
mag_diff_region2 = abs(CFHT_i_region2 - F814W_region2[idx_2])
mag_diff_region3 = abs(CFHT_i_region3 - F814W_region3[idx_3])
mag_diff_region4 = abs(CFHT_i_region4 - F814W_region4[idx_4])

#we only care about the HST data that has been matched
HST_RA_region1 = HST_RA_region1[idx_1]
HST_Dec_region1 = HST_Dec_region1[idx_1]
F814W_region1 = F814W_region1[idx_1]

HST_RA_region2 = HST_RA_region2[idx_2]
HST_Dec_region2 = HST_Dec_region2[idx_2]
F814W_region2 = F814W_region2[idx_2]

HST_RA_region3 = HST_RA_region3[idx_3]
HST_Dec_region3 = HST_Dec_region3[idx_3]
F814W_region3 = F814W_region3[idx_3]

HST_RA_region4 = HST_RA_region4[idx_4]
HST_Dec_region4 = HST_Dec_region4[idx_4]
F814W_region4 = F814W_region4[idx_4]

#determine which matches are good
good_match_1 = (d2d_1.arcsecond < 1) & (mag_diff_region1 < 0.5)
HST_RA_region1 = HST_RA_region1[good_match_1]
HST_Dec_region1 = HST_Dec_region1[good_match_1]
F814W_region1 = F814W_region1[good_match_1]
CFHT_RA_region1 = CFHT_RA_region1[good_match_1]
CFHT_Dec_region1 = CFHT_Dec_region1[good_match_1]
CFHT_i_region1 = CFHT_i_region1[good_match_1]

good_match_2 = (d2d_2.arcsecond < 1) & (mag_diff_region2 < 0.5)
HST_RA_region2 = HST_RA_region2[good_match_2]
HST_Dec_region2 = HST_Dec_region2[good_match_2]
F814W_region2 = F814W_region2[good_match_2]
CFHT_RA_region2 = CFHT_RA_region2[good_match_2]
CFHT_Dec_region2 = CFHT_Dec_region2[good_match_2]
CFHT_i_region2 = CFHT_i_region2[good_match_2]

good_match_3 = (d2d_3.arcsecond < 1) & (mag_diff_region3 < 0.5)
HST_RA_region3 = HST_RA_region3[good_match_3]
HST_Dec_region3 = HST_Dec_region3[good_match_3]
F814W_region3 = F814W_region3[good_match_3]
CFHT_RA_region3 = CFHT_RA_region3[good_match_3]
CFHT_Dec_region3 = CFHT_Dec_region3[good_match_3]
CFHT_i_region3 = CFHT_i_region3[good_match_3]

good_match_4 = (d2d_4.arcsecond < 1) & (mag_diff_region4 < 0.5)
HST_RA_region4 = HST_RA_region4[good_match_4]
HST_Dec_region4 = HST_Dec_region4[good_match_4]
F814W_region4 = F814W_region4[good_match_4]
CFHT_RA_region4 = CFHT_RA_region4[good_match_4]
CFHT_Dec_region4 = CFHT_Dec_region4[good_match_4]
CFHT_i_region4 = CFHT_i_region4[good_match_4]

# '''
# ============================================================================================================================
# save the data
# '''

np.savetxt('/Users/amandaquirk/Documents/M33/Data/CFHT_region1.txt', np.c_[CFHT_RA_region1, CFHT_Dec_region1, CFHT_i_region1], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), i (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/HST_region1.txt', np.c_[HST_RA_region1, HST_Dec_region1, F814W_region1], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), F814W (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/CFHT_region2.txt', np.c_[CFHT_RA_region2, CFHT_Dec_region2, CFHT_i_region2], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), i (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/HST_region2.txt', np.c_[HST_RA_region2, HST_Dec_region2, F814W_region2], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), F814W (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/CFHT_region3.txt', np.c_[CFHT_RA_region3, CFHT_Dec_region3, CFHT_i_region3], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), i (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/HST_region3.txt', np.c_[HST_RA_region3, HST_Dec_region3, F814W_region3], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), F814W (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/CFHT_region4.txt', np.c_[CFHT_RA_region4, CFHT_Dec_region4, CFHT_i_region4], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), i (mag)')
np.savetxt('/Users/amandaquirk/Documents/M33/Data/HST_region4.txt', np.c_[HST_RA_region4, HST_Dec_region4, F814W_region4], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), F814W (mag)')

'''
============================================================================================================================
plotting 
'''
def matching_plot(RA1, RA2, Dec1, Dec2, region_name):
	def deprojected_dRA(RA1, RA2): #deg
		M33_Dec = 30.6603 * np.pi / 180 #radians
		deg = (RA1 - RA2) * np.cos(M33_Dec)
		return deg * 3600 #arcsec

	def dDEC(Dec1, Dec2): #deg
		deg = Dec1 - Dec2
		return deg * 3600 #arcsec

	x = deprojected_dRA(RA1, RA2)
	y = dDEC(Dec1, Dec2)

	f, axes = plt.subplots(1, figsize=(4,4))
	circle= plt.Circle((0,0), 1, color='k', fill=False, linestyle='--')
	axes.add_artist(circle)
	circle_small = plt.Circle((0,0), .1, color='r', fill=False, linestyle='--')
	axes.add_artist(circle_small)
	circle2 = plt.Circle((0,0), .5, color='g', fill=False, linestyle='--')
	axes.add_artist(circle2)
	plt.scatter(x, y, c='b', alpha=0.4)
	plt.xlabel(r'$\rm\Delta RA\times\ DEC\ (arcsec)$')
	plt.ylabel(r'$\rm\Delta DEC\ (arcsec)$')
	plt.xlim(-1.1, 1.1)
	plt.ylim(-1.1, 1.1)
	f.subplots_adjust(left=0.18)
	plt.savefig('/Users/amandaquirk/Desktop/matching_zoom_{}.png'.format(region_name), bbox='tight')
	plt.close()

matching_plot(HST_RA_region1, CFHT_RA_region1, HST_Dec_region1, CFHT_Dec_region1, 'region1')
matching_plot(HST_RA_region2, CFHT_RA_region2, HST_Dec_region2, CFHT_Dec_region2, 'region2')
matching_plot(HST_RA_region3, CFHT_RA_region3, HST_Dec_region3, CFHT_Dec_region3, 'region3')
matching_plot(HST_RA_region4, CFHT_RA_region4, HST_Dec_region4, CFHT_Dec_region4, 'region4')


