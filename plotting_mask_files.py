import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u 

def plot_data(mask_number, color):

	data = np.genfromtxt('/Users/amandaquirk/Desktop/{}M33P.out'.format(mask_number), dtype=None, names='ID, ras, decs, frame, mag, filter, priority, list, selection')
 
	RA =  data['ras'] 
	Dec= data['decs']
	selection_flag = data['selection']

	selected = (selection_flag == 1)
	RA = RA[selected]
	Dec = Dec[selected]

	RA = [a.decode("utf-8") for a in RA]
	Dec = [a.decode("utf-8") for a in Dec]

	coords = SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))

	plt.scatter(coords.ra.deg, coords.dec.deg, c=color, alpha=.4)

	plt.xlim(23.3, 23.7)
	plt.ylim(30.4, 30.9)	

	return 

plot_data('A1', 'r')
plot_data('A2', 'r')
plot_data('B1', 'b')
plot_data('B2', 'b')
plot_data('C1', 'g')
plot_data('C2', 'g')
plot_data('D1', 'm')
plot_data('D2', 'm')
plot_data('E1', 'gray')
plot_data('E2', 'gray')
plot_data('K1', 'k')
plt.show()

#KG data
# data = np.genfromtxt('/Users/amandaquirk/Desktop/Kristen_deimos_target_list.txt', dtype=None, names='ID, filter1, filter2, ras, decs')

# ID = data['ID']
# RA = data['ras'] #does this work or nah? break up into ra and dec?
# Dec = data['decs']

# # coords = SkyCoord(RA, Dec[x], unit=(u.deg, u.deg))
	
# # 	plt.scatter(coords.ra.deg, coords.dec.deg)
# # 	plt.xlim(23.3, 23.7)
# # 	plt.ylim(30.4, 30.9)
# # 	plt.show()

# for i in range(len(ID)):
# 	if (Dec[i] > 30.5) & (Dec[i] < 30.65):
# 		print(ID[i])




