import numpy as np 

'''
purpose: update the target list by removing bright stars that have already been observed (keep guide/alignment stars and keep dim stars)
'''

mask_output_file = 'mask1.out' #file from dsimulator
bright = 0 #magnitude cutoff between bright and dim
new_input_file = 'updated_target_list_mask1'

#read in the data
data = np.genfromtxt('/Users/amandaquirk/Desktop/{}'.format(mask_output_file), dtype=None, names='ID, ras, decs, frame, mag, filter, priority, list, selection')

ID = data['ID'] 
RA = data['ras'] #does this work or nah? break up into ra and dec?
Dec = data['decs']
coordinate_reference_frame = data['frame'] 
magnitude = data['mag'] 
filter_tag = data['filter']
priority = data['priority']
list_assignment = data['list'] 
selection_flag = data['selection']

all_IDs = []
all_ras = []
all_decs = []
all_coord_frame = []
all_mags = []
all_mag_ref = []
all_priorities = []
all_list_assignments = []
all_selection_flag = []
for i in range(len(ID)):
	if (selection_flag[i] == 0.0) | ((selection_flag[i] == 1.0) & (magnitude[i] > bright)) | ((selection_flag[i] == 1.0) & (list_assignment[i] == 0)): #want stars that have not been selected or have been selected but are dim or are our guide/alignment stars
		all_IDs.append(ID[i])
		all_ras.append(RA[i])
		all_decs.append(Dec[i])
		all_coord_frame.append(coordinate_reference_frame[i])
		all_mags.append(magnitude[i])
		all_mag_ref.append(filter_tag[i])
		all_priorities.append(priority[i])
		all_list_assignments.append(list_assignment[i])
		all_selection_flag.append(selection_flag[i])

np.savetxt('/Users/amandaquirk/Desktop/{}.in'.format(new_input_file), np.c_[all_IDs, all_ras, all_decs, all_coord_frame, all_mags, all_mag_ref, all_priorities, all_list_assignments, all_selection_flag], fmt="%-s", delimiter='\t', header='ID, coordinates, coordinate reference frame, magnitude, filter, priority, list assignment, selection flag')  #formats with b'str' -- need to remove


