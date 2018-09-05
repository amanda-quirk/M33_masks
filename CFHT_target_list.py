from astropy.io import fits
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u

'''
purpose: sort CFHT into target lists
'''

#CFHT data-- using the old i filter; this has already been dust corrected
hdu = fits.open('/Users/amandaquirk/Documents/M33/data/M33.GI.matchcat.extcorr.fits')
data = hdu[1].data
RA = data['RA2'] #degrees
Dec = data['DEC2'] #degrees
i_mag = data['MAG2_AUTO_0']
g_mag = data['MAG1_AUTO_0']
isolated_tag = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/CFHT_julia_isolated.txt') #1 is NOT isolated

ID = np.linspace(0, len(RA), len(RA) + 1) 

'''
========================================================================
ELIMINATE DIM STARS
only want stars brighter than i 22 or g 24

do I also want a brightness cut?
'''

bright = (i_mag < 22) | (g_mag < 24)
RA = RA[bright]
Dec = Dec[bright]
i_mag = i_mag[bright]
g_mag = g_mag[bright]
isolated_tag = isolated_tag[bright]
ID = ID[:-1][bright]

'''
========================================================================
ISOLATION
list 0 will be guide/alignment stars, 1 is for HST stars that based the isolation criteria, 2 is for CFHT stars that passed isolation, 3 is HST data that did NOT pass isolation, and 4 is for CFHT stars that did not pass isolation
'''

list_assignment = np.zeros_like(RA) 

for i in range(len(list_assignment)):
	if isolated_tag[i] == 1:
		list_assignment[i] = 4
	else:
		if ((i_mag[i] > 15) & (i_mag[i] < 18)) | ((g_mag[i] > 15) & (g_mag[i] < 18)): #alignment/guide stars
			list_assignment[i] = 0
		else:
			list_assignment[i] = 2

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
priority = np.zeros_like(RA) 

for i in range(len(list_assignment)): #check order of operations
	if list_assignment[i] == 0: #alignment/guide stars 
		priority[i] = -2
	else:
		if ((i_mag[i] > 19) & (i_mag[i] < 20.5)) | ((g_mag[i] > 20) & (g_mag[i] < 22.5)):
			priority[i] = 6
		elif ((i_mag[i] > 20.5) & (i_mag[i] < 21)) | ((g_mag[i] > 22.5) & (g_mag[i] < 23)):
			priority[i] = 4
		elif ((i_mag[i] > 21) & (i_mag[i] < 21.5)) | ((g_mag[i] > 23) & (g_mag[i] < 23.5)):
			priority[i] = 2
		else:
			priority[i] = 0

'''
========================================================================
COORDINATES
need to have the coordinates in HA format for dsimulator

Need to apply any offsets??
'''
coords = SkyCoord(RA, Dec, unit=(u.deg, u.deg))
formated_coords = coords.to_string('hmsdms', alwayssign=False) #will need to edit by hand to change hms dms to : until I find a better way

'''
========================================================================
SAVE TO FILE
'''
    
JD = np.zeros_like(RA) + 2000.00 #coordinate frame reference
filter_tag = ['I' for x in range(0, len(RA))]

np.savetxt('/Users/amandaquirk/Desktop/CFHT_target_list.in', np.c_[ID, formated_coords, JD, i_mag, filter_tag, priority, list_assignment], fmt='%s', delimiter='  ', header='original index, coordinates, coordinate reference frame, magnitude, filter, priority, list assignment') 

