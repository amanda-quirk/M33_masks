import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
import sfdmap 
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import matplotlib.pyplot as plt 
from matplotlib import rc

'''
applies extinction corrections to all of the data filters used to sort into a CMD
resorts stars into age bins
'''

#read in catalogue ===========================================================================================================
data = np.genfromtxt('/Volumes/Titan/M33/Data/M33_2018b_phot_spec.txt', dtype=None, names='ID, ra, dec, F275W, F336W, F475W, F814W, F110W, F160W, z, vel, vel_aband, err, zqual, aband, time, mask, age, HI, CO, Ha')

ID = data['ID']
ID = np.array([a.decode("utf-8") for a in ID])
ra = data['ra']
ra = np.array([a.decode("utf-8") for a in ra])
dec = data['dec']
dec = np.array([a.decode("utf-8") for a in dec])
F275W = data['F275W']
F336W = data['F336W']
F475W = data['F475W']
F814W = data['F814W']
F110W = data['F110W']
F160W = data['F160W']
z = data['z']
vel = data['vel']
vel_aband = data['vel_aband']
err = data['err']
zqual = data['zqual']
aband = data['aband']
time = data['time']
mask = data['mask']
mask = [a.decode("utf-8") for a in mask]
age = data['age']
age = np.array([a.decode("utf-8") for a in age])
HI = data['HI']
CO = data['CO']
Ha = data['Ha']


#read in the extinction map ==================================================================================================
emap = sfdmap.SFDMap('../Data/sfddata-master/', scaling=1.0)
coords = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
extinct_map = emap.ebv(coords)

#filter info
acs_475 = 3.268
acs_606 = 2.471
acs_814 = 1.526
sdss_g = 3.303
sdss_i = 1.698

#apply the extinction corrections
num = [len(x) for x in ID] #length of the ID which will help indicate which catalogue it came from
CFHT = (age == 'CFHT') | (np.array(num) < 6) #F475W ~ SDSS g and F814W ~ SDSS i 

mag1_ext = np.zeros(len(F475W))
mag2_ext = np.zeros(len(F475W))
mag1_name = np.zeros_like(ID)
mag2_name = np.zeros_like(ID)
for i in range(len(CFHT)):
	if CFHT[i] == True:
		mag1_ext[i] = F475W[i] - extinct_map[i] * sdss_g
		mag2_ext[i] = F814W[i] - extinct_map[i] * sdss_i
		mag1_name[i] = 'g'
		mag2_name[i] = 'i'
	elif ID[i].isnumeric() == True:
		if ((int(ID[i]) < 8000000) | (int(ID[i]) >= 11000000)): #2016 HST and blue filter is F606W
			mag1_ext[i] = F475W[i] - extinct_map[i] * acs_606
			mag2_ext[i] = F814W[i] - extinct_map[i] * acs_814
			mag1_name[i] = 'F606W'
			mag2_name[i] = 'F814W'
		else: #2016 HST and blue filter is F475W
			mag1_ext[i] = F475W[i] - extinct_map[i] * acs_475
			mag2_ext[i] = F814W[i] - extinct_map[i] * acs_814
			mag1_name[i] = 'F475W'
			mag2_name[i] = 'F814W'
	else: #2018/19 HST and blue filter is F475W
		mag1_ext[i] = F475W[i] - extinct_map[i] * acs_475
		mag2_ext[i] = F814W[i] - extinct_map[i] * acs_814
		mag1_name[i] = 'F475W'
		mag2_name[i] = 'F814W'

#make a CMD to add age tags to the untagged data =============================================================================
color = mag1_ext - mag2_ext

# rc('font', family = 'serif')
# fig, ax=plt.subplots(1)
# for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(1)
# ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
# ax.tick_params(axis='x',which='both',top='on', direction='in')
# ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# ax.tick_params(axis='y',which='both',right='on', direction='in')
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# plt.tick_params(labelsize=12) 
# plt.minorticks_on()

# plt.scatter(color[CFHT], mag2_ext[CFHT], s=2, alpha=.5)
# plt.gca().invert_yaxis()
# plt.xlabel('g - i')
# plt.ylabel('i')
# plt.savefig('/Users/amandaquirk/Desktop/CFHT_CMD.png')


#CMD sorting =================================================================================================================
'''
currently plan: 
	-keep all of Anil's labels the same for the 18-19 HST
	-for the CFHT 16-19 data, use the CMDs that I worked with with Raja for the young stars and AGB; use Karrie's notebook for a rough RGB polygon
	-for the 16 HST, use the Karrie's polygons and try sorting the others based on their extinction corrected CMDs
'''
age_tag = np.zeros_like(age)

for i in len(age):
	#going dataset by dataset or I will get confused
	if ID[i].isnumeric() == False: 	#2018/2019 data
		if ID[i].startswith('CFHT') == False: #HST data -- keeping Anil's label
			age_tag[i] = age[i]
		else: #CFHT data -- using Karrie's RGB selection box and CMD 
			point = Point(color[i], mag2_ext[i])
			RGB_poly = Polygon[(1.4, 21), (2.4, 21.25), (3.55, 21.9), (3, 21.89), (2.6, 21.95), (1.9, 22.5), (1.8, 22.8), (1.5, 23.4), (1.2, 24.9), (.43, 23.25), (1.1, 22)]
			HeB_poly = Polygon[(.8, 22), (1.1, 19.5), (1.4, 18), (1.8, 19), (1.6, 20), (1.1, 22)]
			AGB_poly = Polygon[(1.4, 21), (2.4, 21.25), (3.55, 21.9), (3.55, 20), (1.8, 20)]
			if RGB_poly.contains(point) == True:
				age_tag[i] = 'RGB'
			elif HeB_poly.contains(point) == True:
				age_tag[i] = 'HeB'
			elif (color[i] < -0.5) & (mag2_ext[i] > 20):
				age_tag[i] = 'MS'
			elif AGB_poly.contains(point) == True:
				age_tag[i] = 'AGB'
			else:
				age_tag[i] = age[i]
	elif num[i] < 6: #2016 CHFT data -- same as above
		point = Point(color[i], mag2_ext[i])
		RGB_poly = Polygon[(1.4, 21), (2.4, 21.25), (3.55, 21.9), (3, 21.89), (2.6, 21.95), (1.9, 22.5), (1.8, 22.8), (1.5, 23.4), (1.2, 24.9), (.43, 23.25), (1.1, 22)]
		HeB_poly = Polygon[(.8, 22), (1.1, 19.5), (1.4, 18), (1.8, 19), (1.6, 20), (1.1, 22)]
		AGB_poly = Polygon[(1.4, 21), (2.4, 21.25), (3.55, 21.9), (3.55, 20), (1.8, 20)]
		if RGB_poly.contains(point) == True:
			age_tag[i] = 'RGB'
		elif HeB_poly.contains(point) == True:
			age_tag[i] = 'HeB'
		elif (color[i] < -0.5) & (mag2_ext[i] > 20):
			age_tag[i] = 'MS'
		elif AGB_poly.contains(point) == True:
			age_tag[i] = 'AGB'
		else:
			age_tag[i] = age[i]
	elif ((int(ID[i]) < 8000000) | (int(ID[i]) >= 11000000)): #2016 HST and blue filter is F606W
		point = Point(color[i], mag2_ext[i])
		RGB_poly = Polygon[(.8, 23.2), (.9, 22), (1, 20.5), (2.5, 21), (1.8, 21.2), (1.25, 21.8), (1.1, 22.2), (1, 23.25)]
		if RGB_poly.contains(point) == True:
			age_tag[i] = 'RGB'
		elif color[i] < 0.5:
 			age_tag[i] = 'MS'
		else:
			age_tag[i] = age[i]
	else: #2016 HST and blue filter is F475W -- use Anil's color coded CMD
		point = Point(color[i], mag2_ext[i])
		RGB_poly = Polygon[]
		if RGB_poly.contains(point) == True:
			age_tag[i] = 'RGB'
		elif color[i] < 1.19:
 			age_tag[i] = 'MS'
 		else:
 			age_tag[i] = age[i]

#save data ===================================================================================================================
# np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_phot_spec.txt', np.c_[ID, ra, dec, F275W, F336W, mag1_ext, mag1_name, mag2_ext, mag2_name, F110W, F160W, z, vel, vel_aband, err, zqual, aband, time, mask, age_tag, HI, CO, Ha], fmt='%s', delimiter='\t', header='ID, RA, Dec, F275W, F336W, mag1, mag1 name, mag2, mag2 name, F110W, F160W, redshift, heliocorrected vel (km/s), helio+aband corrected vel (km/s), velocity error (km/s), zquality, A band, MJD, mask name, age tag, HI (km/s), CO (km/s), Halpha (km/s)') 




