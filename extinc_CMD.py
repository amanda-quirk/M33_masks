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
# mag1_cfht = Column(tbdata['G_AUTO_CFHT'] - extinct_map*sdss_g, name='G_AUTO_CFHT_0')
# mag2_cfht = Column(tbdata['I_AUTO_CFHT'] - extinct_map*sdss_i, name='I_AUTO_CFHT_0') 

# mag1a_hst16b = Column(tbdata['F475W_ACS'] - extinct_map*acs_475, name='F475W_ACS_0')
# mag1b_hst16b = Column(tbdata['F606W_ACS'] - extinct_map*acs_606, name='F606W_ACS_0')
# mag2_hst16b = Column(tbdata['F814W_ACS'] - extinct_map*acs_814, name='F814W_ACS_0')

# # current photometry file doesn't have error
# mag1_hst18b = Column(tbdata['F475W_VEGA'] - extinct_map*acs_475, name='F475W_VEGA_0')
# mag2_hst18b = Column(tbdata['F814W_VEGA'] - extinct_map*acs_814, name='F814W_VEGA_0')



