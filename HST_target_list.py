from astropy.io import fits
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
import h5py 

'''
purpose: sort HST into target lists

TO DO: add guide and alignment stars; what is our mag goal?; from HST for CFHT?
'''

#will need to reformat for Anil's file-- more variables needed?
hdu = fits.open('datafile', memmap=True) 
data = hdu[1].data 
RA = data['ra'] #degrees
Dec = data['dec'] #degrees
F814W = data['f814w_vega']
F475W = data['f475w_vega']
label = *******
A_priority = ******

input_file1 = h5py.File('/Users/amandaquirk/Documents/M33/Data/HST_julia_isolated.hdf5', 'r') #use 1) saving_HST.py and 2) find_neigh.jl to create this file
isolated_tag1 = input_file1["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated 
input_file2 = h5py.File('/Users/amandaquirk/Documents/M33/Data/HST_brick23_julia_isolated.hdf5', 'r') 
isolated_tag2 = input_file2["isolation_tag"][...] 

#combine isolation tags
isolation_tag = np.concatenate((isolated_tag1, isolated_tag2), axis=None)
print('Data loaded')

'''
========================================================================
ELIMINATE DIM STARS
only want stars brighter than F814W 22 or F475W 24
'''

bright = (F814W < 22) | (F475W < 24)
RA = RA[bright]
Dec = Dec[bright]
F814W = F814W[bright]
F475W = F475W[bright]
label = label[bright]
isolated_tag = isolated_tag[bright]
A_priority = A_priority[bright]
print('Dim stars eliminated')

'''
========================================================================
ELIMINATE DIM MWFG Stars
-as tagged by Anil
-we want to keep bright ones for guidestars
'''

M33 = (label.startswith('MWFG') == False) | ((label.startswith('MWFG') == True) & (((F475W[i] > 15) & (F475W[i] < 18)) | ((F814W[i] > 15) & (F814W[i] < 18)))) *********
RA = RA[M33]
Dec = Dec[M33]
F814W = F814W[M33]
F475W = F475W[M33]
label = label[M33]
isolated_tag = isolated_tag[M33]
A_priority = A_priority[M33]
print('Dim MWFG stars eliminated')

'''
========================================================================
ISOLATION
list 0 will be guide/alignment stars, 1 is for HST stars that based the isolation criteria, 2 is for CFHT stars that passed isolation, 3 is HST data that did NOT pass isolation, and 4 is for CFHT stars that did not pass isolation
'''

list_assignment = np.zeros_like(RA) 

for i in range(len(list_assignment)):
	if isolated_tag[i] == 1:
		list_assignment[i] = 3
	else: #guide/alignment stars can be the bright MWFG stars saved above OR bright HST stars
		if (label.startswith('MWFG') == True) | ((label.startswith('MWFG') == False) & (((F475W[i] > 15) & (F475W[i] < 18)) | ((F814W[i] > 15) & (F814W[i] < 18)))):
			list_assignment[i] = 0
		else:
			list_assignment[i] = 1
print('List assignment done')

'''
========================================================================
PRIORITY
will contain the priority of stars: for CFHT based on brightness, for HST is based on CMD space
HST: (see how this compares to Anil's priorities)
-10 for specifically targeted stars (like the HMXB)
-8, 6, 4, 2 for rarest stars to RGB stars respectively 
-0 for MWFG stars
CFHT:
-6: 19 < i < 20.5 or 20 < g < 22.5 
-4: 20.5 < i < 21. or 22.5 < g < 23 
-2: 21 < i < 21.5 or 23 < g > 23.5 
-0: i > 22 or g > 24
(-1): guidestar
(-2): alignment star
'''
priority = np.zeros_like(RA) *********
#is there a faster way to do this?? can i just use Anil's priorities? are there more labels? try to just use Anil's labels so as to avoid these loops; will need to keep the designation of guide/alignment stars so try to eliminated the second big set of loops
for i in range(len(RA)):
	if list_assignment[i] == 0:
		if (label.startswith('MWFG') == True):
			priority[i] = -1
		else:
			priority[i] = -2
	else:
		if label[i].startswith('MS'):
			if (F814W[i] > 21) | (F475W[i] > 23):
				priority[i] = 2
			else:
				priority[i] = 6
		elif label[i].startswith('AGB'):
				priority[i] = 8
		elif label[i].startswith('RGB'):
			if F814W[i] > 21.5:
				priority[i] = 2
			else:
				priority = 4
print('Priorities set')

'''
========================================================================
COORDINATES
need to have the coordinates in HA format for dsimulator
'''
coords = SkyCoord(RA, Dec, unit=(u.deg, u.deg))
formated_coords = coords.to_string('hmsdms', alwayssign=False) #will need to edit by hand to change hms dms to : until I find a better way
print('Coordinates created')

'''
========================================================================
SAVE TO FILE
'''
ID = np.linspace(0, len(RA), len(RA) + 1)
JD = np.zeros_like(RA) + 2000.00 #coordinate frame reference
filter_tag = ['F814W' for x in range(0, len(RA))]

#is this going to be too big as a text file?
np.savetxt('/Users/amandaquirk/Desktop/HST_target_list.in', np.c_[ID[:-1], formated_coords, JD, F814W, filter_tag, priority, list_assignment], fmt='%s', delimiter='  ', header='original index, coordinates, coordinate reference frame, magnitude, filter, priority, list assignment') 

