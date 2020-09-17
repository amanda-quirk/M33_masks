import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
import sfdmap 

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

#CMD sorting =================================================================================================================
#currently using a combo of Anil's CMD and what we used in the RGB selection for the 2019 target lists
color = mag1_ext - mag2_ext
age_tag = np.zeros_like(age)

for i in len(age):
	if mag1_name[i] == 'F475W': #can sort all age bins
		if color[i] < 1.19:
			age_tag[i] = 'MS'
	# elif mag1_name[i] == 'F606W': #can sort RGB only (for now?)
	# 	if RGB selction:
	# 		age_tag[i] = 'RGB'
	# 	else:
	# 		age_tag = 'unknown'
	# elif mag1_name[i] == 'g': #can sort RGB only (for now?)
	# 	if RGB selction:
	# 		age_tag[i] = 'RGB'
	# 	else:
	# 		age_tag = 'unknown'


#save data ===================================================================================================================
# np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_phot_spec.txt', np.c_[ID, ra, dec, F275W, F336W, mag1_ext, mag1_name, mag2_ext, mag2_name, F110W, F160W, z, vel, vel_aband, err, zqual, aband, time, mask, age_tag, HI, CO, Ha], fmt='%s', delimiter='\t', header='ID, RA, Dec, F275W, F336W, mag1, mag1 name, mag2, mag2 name, F110W, F160W, redshift, heliocorrected vel (km/s), helio+aband corrected vel (km/s), velocity error (km/s), zquality, A band, MJD, mask name, age tag, HI (km/s), CO (km/s), Halpha (km/s)') 




