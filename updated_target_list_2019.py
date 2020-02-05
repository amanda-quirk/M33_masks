import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u

#read in all data
#read in all data
ref_data = np.genfromtxt('/Users/amandaquirk/Desktop/target_list_RGB_2019_2mag.in', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, magnitude2, filter2, priority, list_assignment, selection_flag')
ref_ID = ref_data['ID']
ref_ID = [a.decode("utf-8") for a in ref_ID]
mag1 = ref_data['magnitude1']
mag2 = ref_data['magnitude2']
RA = ref_data['ras']
RA = [a.decode("utf-8") for a in RA]
Dec = ref_data['decs']
Dec = [a.decode("utf-8") for a in Dec]
frame = ref_data['frame']
filter1 = ref_data['filter1']
filter1 = [a.decode("utf-8") for a in filter1]
priority_0 = ref_data['priority']
list_assignment_0 = ref_data['list_assignment']

#read in selected IDs
sel_data = np.genfromtxt('/Users/amandaquirk/Desktop/selected_set2.txt', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, priority, list_assignment, selection_flag')
sel_ID = sel_data['ID']
sel_ID = [a.decode("utf-8") for a in sel_ID]

#move observed stars to lower lists
for i in range(len(sel_ID)):
	if sel_ID[i] in ref_ID:
		N = ref_ID.index(sel_ID[i])
		if sel_ID[i].startswith('CFHT'):
			priority_0[N] = 1
		else:
			list_assignment_0[N] = 2 

#move stars not observed to correct list
for i in range(len(ref_ID)):
	if ref_ID[i] not in sel_ID:
		#if (ref_ID[i] == 'AGB') | (ref_ID[i].startswith('RHeB')):
		#	list_assignment_0[i] = 3
		if ref_ID[i].startswith('CFHT'):
			list_assignment_0[i] = 3

#try to prioritize some bluer TRGB stars
# color = mag2 - mag1
# for i in range(len(ref_ID)):
# 	if (ref_ID[i].startswith('RHeB')) | (ref_ID[i] == 'AGB') | (ref_ID[i] == 'RGB'):
# 		if (1.5 < color[i]) & (color[i] < 2) & (20 < mag1[i]) & (mag1[i] < 21):
# 			priority_0[i] = 2100

coords = SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))
coords = coords.to_string('hmsdms', alwayssign=False)
selection_flag = np.zeros(len(priority_0))

np.savetxt('/Users/amandaquirk/Desktop/target_list_RGB_2019_updated.in', np.c_[ref_ID, coords, frame, mag1, filter1, priority_0, list_assignment_0, selection_flag], fmt="%-s", delimiter='\t', header='ID, coordinates, coordinate reference frame, magnitude1, filter1, priority, list assignment, selection flag') 