import numpy as np 
import matplotlib.pyplot as plt 

#read in all data
ref_data = np.genfromtxt('/Users/amandaquirk/Desktop/target_list_RGB_2019_2mag.in', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, magnitude2, filter2, priority, list_assignment, selection_flag')
ref_ID = ref_data['ID']
ref_ID = [a.decode("utf-8") for a in ref_ID]
mag1 = ref_data['magnitude1']
mag2 = ref_data['magnitude2']

#read in selected IDs
sel_data = np.genfromtxt('/Users/amandaquirk/Desktop/selected.txt', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, priority, list_assignment, selection_flag')
sel_ID = sel_data['ID']
sel_ID = [a.decode("utf-8") for a in sel_ID]

#read in not selected IDs
# not_sel_data = np.genfromtxt('/Users/amandaquirk/Desktop/not_selected.txt', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, priority, list_assignment, selection_flag')
# not_sel_ID = not_sel_data['ID']
# not_sel_ID = [a.decode("utf-8") for a in not_sel_ID]

sel_i = []
sel_g = []
sel_F475W = []
sel_F814W = []
for name in sel_ID:
	if name.startswith('CFHT'):
		N = ref_ID.index(name)
		# sel_i.append(mag1[N])
		# sel_g.append(mag2[N])
	else:
		N = ref_ID.index(name)
		sel_F814W.append(mag1[N])
		sel_F475W.append(mag2[N])

# not_sel_i = []
# not_sel_g = []
# not_sel_F475W = []
# not_sel_F814W = []
# for name in not_sel_ID:
# 	if name.startswith('CFHT'):
# 		N = ref_ID.index(name)
# 		# not_sel_i.append(mag1[N])
# 		# not_sel_g.append(mag2[N])
# 	else:
# 		N = ref_ID.index(name)
# 		not_sel_F814W.append(mag1[N])
# 		not_sel_F475W.append(mag2[N])

# plt.scatter(np.array(not_sel_g) - np.array(not_sel_i), not_sel_i, alpha=.2, c='grey')
# plt.scatter(np.array(sel_g) - np.array(sel_i), sel_i, alpha=.6, c='r')
# plt.ylabel('i')
# plt.xlabel('g - i')
# plt.xlim(-4, 6)
# plt.ylim(24.2, 12)
# #plt.show()
# plt.savefig('/Users/amandaquirk/Desktop/CFHT_CMD.png')
# plt.close()

#plt.scatter(np.array(not_sel_F475W) - np.array(not_sel_F814W), not_sel_F814W, alpha=.2, c='grey')
plt.scatter(np.array(sel_F475W) - np.array(sel_F814W), sel_F814W, c='r', alpha=.3)
plt.ylabel('F814W')
plt.xlabel('F475W - F814W')
plt.xlim(1, 5)
plt.ylim(22.5, 18)
#plt.show()
plt.savefig('/Users/amandaquirk/Desktop/HST_CMD.png', transparent=True)
plt.close()
