from astropy.io import fits
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
import h5py

'''
purpose: create a list of science targets for the 2018B Keck season
'''

'''
========================================================================================================================
CFHT data
using the old i filter; this has already been dust corrected
'''

hdu = fits.open('/Users/amandaquirk/Documents/M33/data/M33.GI.matchcat.extcorr.fits')
data = hdu[1].data
CFHT_RA = data['RA2'] #degrees
CFHT_Dec = data['DEC2'] #degrees
CFHT_i_mag = data['MAG2_AUTO_0']
CFHT_g_mag = data['MAG1_AUTO_0']
CFHT_isolated_tag = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/CFHT_julia_isolated.txt') #1 is NOT isolated
CFHT_input_file = h5py.File('/Users/amandaquirk/Documents/M33/Data/CFHT_strict_isolated.hdf5', 'r') #use 1) saving_HST.py and 2) find_neigh.jl to create this file
CFHT_strict_isolated_tag = CFHT_input_file["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated 
CFHT_color = CFHT_g_mag - CFHT_i_mag
 
CFHT_ID = []
for i in range(len(CFHT_RA)):
	CFHT_ID.append('CFHT{}'.format(int(i)))

'''
========================================================================
ELIMINATE not crowded stars; only want to keep CFHT stars that could be alignment stars or that pass the color cut
'''

bright = ((((CFHT_i_mag > 15) & (CFHT_i_mag < 19.5)) | ((CFHT_g_mag > 15) & (CFHT_g_mag < 19.5))) | (CFHT_color > 0.5)) & (CFHT_isolated_tag == 0) & (CFHT_i_mag < 24)
CFHT_RA = CFHT_RA[bright]
CFHT_Dec = CFHT_Dec[bright]
CFHT_i_mag = CFHT_i_mag[bright]
CFHT_g_mag = CFHT_g_mag[bright]
CFHT_strict_isolated_tag = CFHT_strict_isolated_tag[bright]
CFHT_ID = np.array((CFHT_ID))[bright]
CFHT_color = CFHT_color[bright]

'''
========================================================================
ISOLATION
list 0 will be guide/alignment stars, 1 is for HST stars that based the isolation criteria, 2 is for CFHT stars that passed the strict isolation, 3 is HST data that did NOT pass isolation, and 4 is for CFHT stars that did not pass the strict isolation
'''

CFHT_list_assignment = np.zeros_like(CFHT_RA) 

for i in range(len(CFHT_list_assignment)):
	# if (CFHT_strict_isolated_tag[i] == 1) | (CFHT_color[i] < 0.5): #super deprioritize stars that are not isolated and do not pass the color cut
	# 	CFHT_list_assignment[i] = 4
	# else:
	if ((CFHT_i_mag[i] > 15) & (CFHT_i_mag[i] < 19.5)) | ((CFHT_g_mag[i] > 15) & (CFHT_g_mag[i] < 19.5)): #alignment/guide stars
		CFHT_list_assignment[i] = 0
	else:
		CFHT_list_assignment[i] = 2

'''
========================================================================
PRIORITY
will contain the priority of stars: for CFHT based on brightness, for HST is based on CMD space
HST: (see how this compares to Anil's priorities)
-200 for specifically targeted stars (like the HMXB)
-10 to 1 for rarest stars to RGB stars respectively 
-0 for MWFG stars
CFHT:
-6: 19 < i < 20.5 or 20 < g < 22.5 
-4: 20.5 < i < 21. or 22.5 < g < 23 
-2: 21 < i < 21.5 or 23 < g > 23.5 
-0: i > 22 or g > 24
(-1): guidestar
(-2): alignment star
'''
CFHT_priority = np.zeros_like(CFHT_RA) 

for i in range(len(CFHT_list_assignment)): 
	if CFHT_list_assignment[i] == 0: #alignment/guide stars 
		CFHT_priority[i] = -2
	else:#elif CFHT_color[i] > 0.5: #want stars that pass the color cut
		if ((CFHT_i_mag[i] > 21) & (CFHT_i_mag[i] < 22)):
			CFHT_priority[i] = 60
		elif ((CFHT_i_mag[i] > 20.5) & (CFHT_i_mag[i] < 21)):
			CFHT_priority[i] = 10
		elif ((CFHT_i_mag[i] > 22) & (CFHT_i_mag[i] < 22.5)):
			CFHT_priority[i] = 10
		else:
			CFHT_priority[i] = 2
	#else:
	#	CFHT_priority[i] = 1
'''
========================================================================
COORDINATES
need to have the coordinates in HA format for dsimulator
'''

CFHT_coords = SkyCoord(CFHT_RA, CFHT_Dec, unit=(u.deg, u.deg))
CFHT_formated_coords = CFHT_coords.to_string('hmsdms', alwayssign=False) #will need to edit by hand to change hms dms to : until I find a better way

'''
========================================================================
TAGS
'''
    
CFHT_JD = np.zeros_like(CFHT_RA) + 2000.00 #coordinate frame reference
CFHT_filter_tag1 = ['I' for x in range(0, len(CFHT_RA))]
CFHT_filter_tag2 = ['G' for x in range(0, len(CFHT_RA))]
CFHT_selection_tag = np.zeros_like(CFHT_RA)


'''
========================================================================================================================
HST data
'''

HST_hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/fall18_refine.fits')
HST_data = HST_hdu[1].data 
HST_RA = HST_data['RA']
HST_Dec = HST_data['Dec']
HST_F814W = HST_data['F814W_VEGA']
HST_F475W = HST_data['F475W_VEGA']
HST_F275W = HST_data['F275W_VEGA']
HST_F336W = HST_data['F336W_VEGA']
HST_F110W = HST_data['F110W_VEGA']
HST_F160W = HST_data['F160W_VEGA']
HST_F814W_crowd = HST_data['F814CROWDMAG']
HST_F475W_crowd = HST_data['F475CROWDMAG']
HST_ID = HST_data['TARGTYPE']
HST_FeH = HST_data['FEH']
HST_ind = np.arange(0, len(HST_RA)) #original index 

input_file = h5py.File('/Users/amandaquirk/Documents/M33/Data/HST_anil_strict_julia_isolated.hdf5', 'r')
HST_isolated_tag = input_file["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated

'''
===============================================================================
ELIMINATE CROWDED STARS AND DIM Stars AND non RGB stars
'''

HST_isolated = (HST_F814W < 24) &  (HST_isolated_tag == 0) #(HST_F475W < HST_F475W_crowd)

HST_RA = HST_RA[HST_isolated]
HST_Dec = HST_Dec[HST_isolated]
HST_F814W = HST_F814W[HST_isolated]
HST_F814W_crowd = HST_F814W_crowd[HST_isolated]
HST_F475W = HST_F475W[HST_isolated]
HST_F275W = HST_F275W[HST_isolated]
HST_F336W = HST_F336W[HST_isolated]
HST_F110W = HST_F110W[HST_isolated]
HST_F160W = HST_F160W[HST_isolated]
HST_ID = HST_ID[HST_isolated]
HST_FeH = HST_FeH[HST_isolated]
HST_ind = HST_ind[HST_isolated]

#just want RGB stars
RGB = (HST_ID == 'RGB') | ((HST_ID == 'AGB') & ((HST_F475W - HST_F814W > 1.5) & (4.5 > HST_F475W - HST_F814W))) | ((HST_ID.startswith('RHeB')==True) & ((HST_F475W - HST_F814W > 1.5) & (4.5 > HST_F475W - HST_F814W)))

print('Number of isolated potential RGB stars is {}'.format(sum(RGB)))

HST_RA = HST_RA[RGB]
HST_Dec = HST_Dec[RGB]
HST_F814W = HST_F814W[RGB]
HST_F814W_crowd = HST_F814W_crowd[RGB]
HST_F475W = HST_F475W[RGB]
HST_F275W = HST_F275W[RGB]
HST_F336W = HST_F336W[RGB]
HST_F110W = HST_F110W[RGB]
HST_F160W = HST_F160W[RGB]
HST_ID = HST_ID[RGB]
HST_FeH = HST_FeH[RGB]
HST_ind = HST_ind[RGB]

'''
============================================================================
ISOLATION
list 0 will be guide/alignment stars, 1 is for HST stars that passed the strict isolation criteria, 2 is for CFHT stars that passed the strict isolation, 3 is HST data that did NOT pass the strict isolation, and 4 is for CFHT stars that did not pass the strict isolation
'''

HST_list_assignment = np.zeros_like(HST_RA) + 1

# for i in range(len(HST_list_assignment)):
# 	if HST_ID[i] == 'RGB':
# 		HST_list_assignment[i] = 1
# 	else:
# 		HST_list_assignment[i] = 2

# '''
# ========================================================================
# PRIORITY
# will contain the priority of stars: for CFHT based on brightness, for HST is based on CMD space
# (-1): guidestar
# (-2): alignment star
# '''

HST_priority = np.zeros_like(HST_RA)

for i in range(len(HST_RA)): #prioritize bright RGB
	if (HST_F814W[i] > 20.6) & (HST_F814W[i] < 21.1):
		HST_priority[i] = 600
	elif (HST_F814W[i] > 21.1) & (HST_F814W[i] < 21.6):
		HST_priority[i] = 250
	elif (HST_F814W[i] > 20.3) & (HST_F814W[i] < 20.6):
		HST_priority[i] = 150
	elif (HST_F814W[i] > 21.6) & (HST_F814W[i] < 22):
		HST_priority[i] = 150
	else:
		HST_priority[i] = 70

	if HST_ID[i] == 'RGB':
		HST_priority[i] = HST_priority[i] + 70

'''
========================================================================
TAGS
'''

HST_JD = np.zeros_like(HST_RA) + 2000.00 #coordinate frame reference
HST_filter_tag1 = ['F814W' for x in range(0, len(HST_RA))]
HST_filter_tag2 = ['F475W' for x in range(0, len(HST_RA))]
HST_selection_tag = np.zeros_like(HST_RA) #want these on the mask
HST_labels = []
for i in range(len(HST_RA)):
	HST_labels.append('{}_{}'.format(HST_ID[i], HST_ind[i]))

'''
========================================================================
COORDINATES
adjust and format
'''
#reference data below made with HST_coordinate_adjust.py
reference_RA, reference_Dec, delta_RA, delta_Dec = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/matched_HST_adjustments.txt', unpack=True)
adjusted_HST_RA = np.zeros_like(HST_RA) #deg
adjusted_HST_Dec = np.zeros_like(HST_Dec) #deg
HST_c = SkyCoord(ra=HST_RA, dec=HST_Dec, unit=(u.deg, u.deg))
reference_c = SkyCoord(ra=reference_RA, dec=reference_Dec, unit=(u.deg, u.deg))
for i in range(len(HST_RA)):
	c1 = HST_c[i]
	sep = c1.separation(reference_c)
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

HST_coords = SkyCoord(adjusted_HST_RA, adjusted_HST_Dec, unit=(u.deg, u.deg))
HST_formated_coords = HST_coords.to_string('hmsdms', alwayssign=False) #will need to edit by hand to change hms dms to : until I find a better way

'''
========================================================================================================================
COMBINING AND SAVING DATA
'''

#mags for CHFT
CFHT_mag_extra = np.ones(len(CFHT_ID)) * np.nan

all_IDs = np.concatenate((HST_labels, CFHT_ID), axis = None)
all_coords = np.concatenate((HST_formated_coords, CFHT_formated_coords), axis = None)
all_coord_frame = np.concatenate((HST_JD, CFHT_JD), axis = None)
all_mags1 = np.concatenate((HST_F814W, CFHT_i_mag), axis = None)
all_mag_ref1 = np.concatenate((HST_filter_tag1, CFHT_filter_tag1), axis = None)
all_mags2 = np.concatenate((HST_F475W, CFHT_g_mag), axis = None)
all_mag_ref2 = np.concatenate((HST_filter_tag2, CFHT_filter_tag2), axis = None)
all_priorities = np.concatenate((HST_priority, CFHT_priority), axis = None)
all_list_assignments = np.concatenate((HST_list_assignment, CFHT_list_assignment), axis = None)
all_selection_flag = np.concatenate((HST_selection_tag, CFHT_selection_tag), axis = None)
all_mags5 = np.concatenate((HST_F275W, CFHT_mag_extra), axis = None)
all_mags6 = np.concatenate((HST_F336W, CFHT_mag_extra), axis = None)
all_mags3 = np.concatenate((HST_F110W, CFHT_mag_extra), axis = None)
all_mags4 = np.concatenate((HST_F160W, CFHT_mag_extra), axis = None)


np.savetxt('/Users/amandaquirk/Desktop/target_list_RGB_2019_expanded.in', np.c_[all_IDs, all_coords, all_coord_frame, all_mags1, all_mag_ref1, all_priorities, all_list_assignments, all_selection_flag], fmt="%-s", delimiter='\t', header='ID, coordinates, coordinate reference frame, magnitude1, filter1, priority, list assignment, selection flag, HST_F110W, HST_F160W, HST_F275W, HST_F336W') 
np.savetxt('//Users/amandaquirk/Desktop/target_list_RGB_2019_expanded_6filt.in', np.c_[all_IDs, all_coords, all_coord_frame, all_mags1, all_mag_ref1, all_mags2, all_mag_ref2, all_priorities, all_list_assignments, all_selection_flag, all_mags3, all_mags4, all_mags5, all_mags6], fmt="%-s", delimiter='\t', header='ID, coordinates, coordinate reference frame, magnitude1, filter1, magnitude2, filter2, priority, list assignment, selection flag, HST_F110W, HST_F160W, HST_F275W, HST_F336W') 

# import matplotlib.pyplot as plt 
# plt.scatter(np.array(CFHT_g_mag) - np.array(CFHT_i_mag), CFHT_i_mag, alpha=.2, c='grey')
# plt.ylabel('i')
# plt.xlabel('g - i')
# plt.xlim(-4, 6)
# plt.ylim(25, 12)
# #plt.show()
# plt.savefig('/Users/amandaquirk/Dropbox/CouchWorking/CFHT_CMD_target.png')
# plt.close()

# plt.scatter(np.array(HST_F475W) - np.array(HST_F814W), HST_F814W, alpha=.2, c='grey')
# plt.ylabel('F814W')
# plt.xlabel('F475W - F814W')
# plt.xlim(1, 5.2)
# plt.ylim(22.5, 18)
# #plt.show()
# plt.savefig('/Users/amandaquirk/Dropbox/CouchWorking/HST_CMD_target.png')
# plt.close()



