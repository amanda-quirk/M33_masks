from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u
import h5py

'''
purpose: correct and examine the offset between the coordinate systems of HST and CFHT

'''

#CFHT data-- using the old i filter; this has already been dust corrected
CFHT_hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/M33.GI.matchcat.extcorr.fits')
CFHT_data = CFHT_hdu[1].data
CFHT_RA = CFHT_data['RA2'] #degrees
CFHT_Dec = CFHT_data['DEC2'] #degrees
CFHT_i = CFHT_data['MAG2_AUTO_0']
CFHT_g = CFHT_data['MAG1_AUTO_0']
CFHT_isolated_tag = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/CFHT_julia_isolated.txt') #1 is NOT isolated

#HST data
f = h5py.File('/Users/amandaquirk/Documents/M33/Data/m33_initsample.hdf5', 'r')
dset = f['df']
dset0 = dset['block0_values']
HST_RA = dset0[:, 90]
HST_Dec = dset0[:, 5]
f475w_vega = dset0[:, 72]
f475crowdmag = dset0[:, 59]
f814w_vega = dset0[:, 86]

input_file = h5py.File('/Users/amandaquirk/Documents/M33/Data/HST_anil_strict_julia_isolated.hdf5', 'r')
HST_isolated_tag = input_file["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated

dset_flag = dset['block1_values']
flag110 = dset_flag[:,0]
flag160 = dset_flag[:,1]
flag275 = dset_flag[:,2]
flag336 = dset_flag[:,3]
flag475 = dset_flag[:,4]
flag814 = dset_flag[:,5]

'''
===================================================================================
for a sanity check, let's see if Anil's and mine isolation is similar
Anil's keeps only ~8,000 stars while my isolation keeps ~20,000
'''

# me_isolated = HST_isolated_tag != 1
# anil_isolated_4 = (f475w_vega < f475crowdmag)

# print(sum(me_isolated))
# print(sum(anil_isolated_4))
# print(sum(anil_isolated_8))
# print(len(ra), len(HST_isolated_tag)) 

'''
===================================================================================
we only want to look at isolated and bright stars
-I will use Anil's criteria for now since it results in fewer stars
-HST has already been cut by magnitude
'''

CFHT_good = (CFHT_isolated_tag != 1) & ((CFHT_i < 22) | (CFHT_g < 24))

HST_good = (f475w_vega < f475crowdmag) & (flag475 == 0) & (flag814 == 0) #((flag110 == 0) | (flag160 == 0) | (flag275 == 0) | (flag336 == 0) | (flag475 == 0) | (flag814 == 0))

CFHT_RA = CFHT_RA[CFHT_good]
CFHT_Dec = CFHT_Dec[CFHT_good]
CFHT_i = CFHT_i[CFHT_good]

HST_RA = HST_RA[HST_good]
HST_Dec = HST_Dec[HST_good]
f814w_vega = f814w_vega[HST_good]

'''
====================================================================================
-do astropy coordinate matching
-only keep matches that are within +/-0.5 mag and are less than 1 arcsec in separation
'''

HST_coords = SkyCoord(HST_RA, HST_Dec, unit=(u.deg, u.deg))
CFHT_coords =SkyCoord(CFHT_RA, CFHT_Dec, unit=(u.deg, u.deg))

#assign each HST data point the closest CFHT point
idx, d2d, d3d = HST_coords.match_to_catalog_sky(CFHT_coords)

#we only care about CFHT data that has been matched
CFHT_RA = CFHT_RA[idx]
CFHT_Dec = CFHT_Dec[idx]
CFHT_i = CFHT_i[idx]

#want to keep matches that are close in mag and distance
mag_diff = abs((f814w_vega - CFHT_i))
good_match = (d2d.arcsecond < 1) & (mag_diff < 0.5)

HST_RA = HST_RA[good_match]
HST_Dec = HST_Dec[good_match]
CFHT_RA = CFHT_RA[good_match]
CFHT_Dec = CFHT_Dec[good_match]

print(len(HST_RA))
'''
====================================================================================
-plotting
'''

def matching_plot(RA1, RA2, Dec1, Dec2, region_name):
	def deprojected_dRA(RA1, RA2): #deg
		#M33_Dec = 30.6603 * np.pi / 180 #radians
		avg_dec = (Dec1 + Dec2) * np.pi / 360 #radians
		deg = (RA1 - RA2) * np.cos(avg_dec)#M33_Dec)
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
	plt.xlabel(r'$\rm\Delta RA\times\ cos(DEC)\ (arcsec)$')
	plt.ylabel(r'$\rm\Delta DEC\ (arcsec)$')
	plt.xlim(-1.1, 1.1)
	plt.ylim(-1.1, 1.1)
	f.subplots_adjust(left=0.18)
	plt.savefig('/Users/amandaquirk/Desktop/matching_coordinates_{}.png'.format(region_name), bbox='tight')
	plt.close()

#examine different regions
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
HST_in_region2 = (HST_RA > RA_min2) & (HST_RA < RA_max2) & (HST_Dec > Dec_min2) & (HST_Dec < Dec_max2) 
CFHT_in_region2 = (CFHT_RA > RA_min2) & (CFHT_RA < RA_max2) & (CFHT_Dec > Dec_min2) & (CFHT_Dec < Dec_max2) 
HST_in_region3 = (HST_RA > RA_min3) & (HST_RA < RA_max3) & (HST_Dec > Dec_min3) & (HST_Dec < Dec_max3) 
CFHT_in_region3 = (CFHT_RA > RA_min3) & (CFHT_RA < RA_max3) & (CFHT_Dec > Dec_min3) & (CFHT_Dec < Dec_max3) 
HST_in_region4 = (HST_RA > RA_min4) & (HST_RA < RA_max4) & (HST_Dec > Dec_min4) & (HST_Dec < Dec_max4) 
CFHT_in_region4 = (CFHT_RA > RA_min4) & (CFHT_RA < RA_max4) & (CFHT_Dec > Dec_min4) & (CFHT_Dec < Dec_max4) 

HST_RA_region1 = HST_RA[HST_in_region1]
HST_Dec_region1 = HST_Dec[HST_in_region1]
CFHT_RA_region1 = CFHT_RA[CFHT_in_region1]
CFHT_Dec_region1 = CFHT_Dec[CFHT_in_region1]

HST_RA_region2 = HST_RA[HST_in_region2]
HST_Dec_region2 = HST_Dec[HST_in_region2]
CFHT_RA_region2 = CFHT_RA[CFHT_in_region2]
CFHT_Dec_region2 = CFHT_Dec[CFHT_in_region2]

HST_RA_region3 = HST_RA[HST_in_region3]
HST_Dec_region3 = HST_Dec[HST_in_region3]
CFHT_RA_region3 = CFHT_RA[CFHT_in_region3]
CFHT_Dec_region3 = CFHT_Dec[CFHT_in_region3]

HST_RA_region4 = HST_RA[HST_in_region4]
HST_Dec_region4 = HST_Dec[HST_in_region4]
CFHT_RA_region4 = CFHT_RA[CFHT_in_region4]
CFHT_Dec_region4 = CFHT_Dec[CFHT_in_region4]

# matching_plot(HST_RA_region1, CFHT_RA_region1, HST_Dec_region1, CFHT_Dec_region1, 'region1')
# matching_plot(HST_RA_region2, CFHT_RA_region2, HST_Dec_region2, CFHT_Dec_region2, 'region2')
# matching_plot(HST_RA_region3, CFHT_RA_region3, HST_Dec_region3, CFHT_Dec_region3, 'region3')
# matching_plot(HST_RA_region4, CFHT_RA_region4, HST_Dec_region4, CFHT_Dec_region4, 'region4')

#===============
# for i in range(len(HST_RA)):
# 	plt.plot([HST_RA[i], CFHT_RA[i]], [HST_Dec[i], CFHT_Dec[i]], c='b')
# plt.xlabel(r'$\rm\ RA\ (deg)$')
# plt.ylabel(r'$\rm\ DEC\ (deg)$')
# plt.show()
# plt.close()

#=================
def deprojected_dRA(RA1, RA2, Dec1, Dec2): #deg
	#M33_Dec = 30.6603 * np.pi / 180 #radians
	avg_dec = (Dec1 + Dec2) * np.pi / 360 #radians
	deg = (RA1 - RA2) * np.cos(avg_dec)#M33_Dec)
	return deg * 3600 #arcsec

def dDEC(Dec1, Dec2): #deg
	deg = Dec1 - Dec2
	return deg * 3600 #arcsec

def arrow_plot(RA1, RA2, Dec1, Dec2, description):

	fig, ax=plt.subplots(1)

	dRA_array = deprojected_dRA(RA1, RA2, Dec1, Dec2)
	dDec_array = dDEC(Dec1, Dec2)

	for i in range(len(RA1)):
		plt.arrow(RA1[i], Dec1[i], dRA_array[i] * 1000 / 3600, dDec_array[i] * 1000 / 3600, color='b', alpha=0.4)
	plt.xlim(23.25, 23.75)
	plt.ylim(30.4, 30.9)
	plt.xlabel(r'$\rm\ RA\ (deg)$')
	plt.ylabel(r'$\rm\ DEC\ (deg)$')
	plt.title(r'$\rm {}\ Correction:\Delta RA\times\ cos(DEC)\ by\ \Delta DEC\ (deg\times\ 1000)$'.format(description))
	ax.annotate('med dRA= {}"'.format(round(np.median(dRA_array),4)), xy=(23.725,30.42), horizontalalignment='right', fontsize=11)
	ax.annotate('med dDec= {}"'.format(round(np.median(dDec_array), 4)), xy=(23.725,30.44), horizontalalignment='right', fontsize=11)

	if description == 'After':
		c = plt.Circle((23.31, 30.47), 200 / 3600, color='k', fill=False)
		ax.add_artist(c)

	plt.savefig('/Users/amandaquirk/Desktop/{}_correction.png'.format(description))
	plt.close()

#arrow_plot(HST_RA, CFHT_RA, HST_Dec, CFHT_Dec, 'Before')

#=================
delta_RA = deprojected_dRA(HST_RA, CFHT_RA, HST_Dec, CFHT_Dec)
delta_Dec = dDEC(HST_Dec, CFHT_Dec)

# plt.scatter(HST_RA, HST_Dec, s=7, c=delta_RA, cmap='Dark2')
# plt.xlabel(r'$\rm\ RA\ (deg)$')
# plt.ylabel(r'$\rm\ DEC\ (deg)$')
# cbr = plt.colorbar()
# cbr.set_label(r'$\rm\Delta RA\times\ cos(DEC)\ (arcsec)$')
# plt.savefig('/Users/amandaquirk/Desktop/delta_RA_map.png', bbox='tight')
# plt.close()

# plt.scatter(HST_RA, HST_Dec, s=7, c=delta_Dec, cmap='Dark2')
# plt.xlabel(r'$\rm\ RA\ (deg)$')
# plt.ylabel(r'$\rm\ DEC\ (deg)$')
# cbr = plt.colorbar()
# cbr.set_label(r'$\rm\Delta DEC\ (arcsec)$')
# plt.savefig('/Users/amandaquirk/Desktop/delta_Dec_map.png', bbox='tight')
# plt.close()

'''
====================================================================================
applying corrections to HST data based on neighbors
if going with this system, will need to save delta_RA and delta_Dec and HST_RA and HST_Dec to output file
'''

adjusted_HST_RA = np.zeros_like(HST_RA)
adjusted_HST_Dec = np.zeros_like(HST_Dec)
c = SkyCoord(ra=HST_RA, dec=HST_Dec, unit=(u.deg, u.deg))
for i in range(len(HST_RA)):
	#print(i)
	c1 = c[i]
	sep = c1.separation(c)
	radius = 200 #arcseconds; initial smoothing size
	close_distances = sep.arcsecond < radius #get everything in this circle
	if sum(close_distances) > 15: #want there to be at least 15 neighbors
		dRAs = delta_RA[close_distances] #arcsecond
		dDecs = delta_Dec[close_distances] #arcsecond
		adjusted_HST_RA[i] = -(np.median(dRAs) / 3600) + HST_RA[i] #degree
		adjusted_HST_Dec[i] = -(np.median(dDecs) / 3600) + HST_Dec[i] #degree
		#print(np.median(dRAs), np.median(dDecs))
	else:
		while sum(close_distances) < 15: #expand the circle until there are 15 neighbors
			radius = radius + 10
			close_distances = sep.arcsecond < radius
		dRAs = delta_RA[close_distances]
		dDecs = delta_Dec[close_distances]
		adjusted_HST_RA[i] = -(np.median(dRAs) / 3600) + HST_RA[i]
		adjusted_HST_Dec[i] = -(np.median(dDecs) / 3600) + HST_Dec[i]
		#print(np.median(dRAs), np.median(dDecs))

#arrow_plot(adjusted_HST_RA, CFHT_RA, adjusted_HST_Dec, CFHT_Dec, 'After')

np.savetxt('/Users/amandaquirk/Desktop/matched_HST_adjustments.txt', np.c_[HST_RA, HST_Dec, delta_RA, delta_Dec], delimiter='  ', header='HST RA (deg), HST Dec (deg), delta RA (arcsecond), delta Dec (arcsecond)') 







