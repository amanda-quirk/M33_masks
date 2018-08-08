from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u
from shapely.geometry import Polygon
from shapely.geometry import Point 

'''
purpose: sort CFHT into target lists
first part of scripts goes through 3 criteria
'''

#CFHT data-- using the old i filter; this has already been dust corrected
CFHT_hdu = fits.open('/Users/amandaquirk/Documents/M33/data/M33.GI.matchcat.extcorr.fits')
CFHT_data = CFHT_hdu[1].data
CFHT_RA = CFHT_data['RA2'] #degrees
CFHT_Dec = CFHT_data['DEC2'] #degrees
CFHT_i = CFHT_data['MAG2_AUTO_0']
CFHT_g = CFHT_data['MAG1_AUTO_0']

#use below just to look for bugs in code
# idx = [0, 1, 2, 3, 4, 5]
# CFHT_RA = CFHT_RA[idx]
# CFHT_Dec = CFHT_Dec[idx]
# CFHT_i = CFHT_i[idx]
# CFHT_g = CFHT_g[idx]

'''
step 1) see if stars are isolated or not =================================================================
-reject a star if it has at least one neighbor that satisfies: I_neighbor < I_target - (d/0.8")^1.5 +2
'''
CFHT_coords = SkyCoord(CFHT_RA, CFHT_Dec, unit=(u.deg, u.deg))

isolated = np.zeros((len(CFHT_RA))) #will contain a 0 if the star doesn't pass the isolation criteria and a 1 if it does
for i in range(len(CFHT_coords)):
	star1 = CFHT_coords[i] #go through star by star
	I_tgt = CFHT_i[i]
	sep = star1.separation(CFHT_coords) #get distances from that star to all other stars in region
	close = (sep.arcsecond < 10) & (sep.arcsecond > 0) #narrow our neighbor search for closer stars to save time; don't want to include the target star
	close_coords = CFHT_coords[close]
	close_mag = CFHT_i[close] 
	close_distances = sep[close]
	reject = I_tgt > close_mag + (close_distances.arcsecond / 0.8)**1.5 - 2 #rejection criteria
	reject_neighbors = close_coords[reject] 
	if len(reject_neighbors) == 0: #don't want to keep any stars that have a neighbor that satisfies the above criteria
		isolated[i] = 1

plt.scatter(CFHT_RA, CFHT_Dec, c=isolated)
plt.xlim(23, 24)
plt.ylim(30, 31)
plt.colorbar()
plt.savefig('/Users/amandaquirk/Desktop/checking_isolation.png')
plt.close()
print('done with isolation!')

'''
step 2) see if stars are in the masks and if so, which mask =================================================================
'''
mask1_center = SkyCoord('1h34m02.7303s', '+30d44m11.000s')
mask2_center = SkyCoord('1h34m07.0646s', '+30d48m07.179s')
mask3_center = SkyCoord('1h33m52.6102s', '+30d32m15.893s')
mask4_center = SkyCoord('1h33m57.4856s', '+30d40m12.844s')
mask5_center = SkyCoord('1h33m55.2257s', '+30d36m14.442s')

half_length = 984 / 2 / 3600 #degrees
half_height = 240 / 2 / 3600 #degrees

RA_min1 = mask1_center.ra.degree - half_length
RA_max1 = mask1_center.ra.degree + half_length
Dec_min1 = mask1_center.dec.degree - half_height
Dec_max1 = mask1_center.dec.degree + half_height

RA_min2 = mask2_center.ra.degree - half_length
RA_max2 = mask2_center.ra.degree + half_length
Dec_min2 = mask2_center.dec.degree - half_height
Dec_max2 = mask2_center.dec.degree + half_height

RA_min3 = mask3_center.ra.degree - half_length
RA_max3 = mask3_center.ra.degree + half_length
Dec_min3 = mask3_center.dec.degree - half_height
Dec_max3 = mask3_center.dec.degree + half_height

RA_min4 = mask4_center.ra.degree - half_length
RA_max4 = mask4_center.ra.degree + half_length
Dec_min4 = mask4_center.dec.degree - half_height
Dec_max4 = mask4_center.dec.degree + half_height

RA_min5 = mask5_center.ra.degree - half_length
RA_max5 = mask5_center.ra.degree + half_length
Dec_min5 = mask5_center.dec.degree - half_height
Dec_max5 = mask5_center.dec.degree + half_height

mask = np.zeros(len(CFHT_RA)) #will contain the mask number (as noted above) that the star is in, 0 if outside of masks
for i in range(len(CFHT_RA)):
	if (CFHT_RA[i] > RA_min1) & (CFHT_RA[i] < RA_max1) & (CFHT_Dec[i] > Dec_min1) & (CFHT_Dec[i] < Dec_max1):
		mask[i] = 1
	elif (CFHT_RA[i] > RA_min2) & (CFHT_RA[i] < RA_max2) & (CFHT_Dec[i] > Dec_min2) & (CFHT_Dec[i] < Dec_max2):
		mask[i] = 2
	elif (CFHT_RA[i] > RA_min3) & (CFHT_RA[i] < RA_max3) & (CFHT_Dec[i] > Dec_min3) & (CFHT_Dec[i] < Dec_max3):
		mask[i] = 3
	elif (CFHT_RA[i] > RA_min4) & (CFHT_RA[i] < RA_max4) & (CFHT_Dec[i] > Dec_min4) & (CFHT_Dec[i] < Dec_max4):
		mask[i] = 4
	elif (CFHT_RA[i] > RA_min5) & (CFHT_RA[i] < RA_max5) & (CFHT_Dec[i] > Dec_min5) & (CFHT_Dec[i] < Dec_max5):
		mask[i] = 5

plt.scatter(CFHT_RA, CFHT_Dec, c=mask)
plt.xlim(23, 24)
plt.ylim(30, 31)
plt.colorbar()
plt.savefig('/Users/amandaquirk/Desktop/checking_mask_placement.png')
plt.close()
print('done with mask assignment!')

'''
step 3) see if stars are bright
'''
mag_cut = np.zeros(len(CFHT_RA)) #will be 0 if the star doesn't pass the magnitude cut and 1 if it does
for i in range(len(CFHT_RA)):
	if (CFHT_i[i] < 22) | (CFHT_g[i] < 24):
		mag_cut[i] = 1

print('done with brightness check!')

'''
step 4) see if stars are in the HST region
'''
#coordinates of HST mask corners
box1_coord1 = SkyCoord('1h34m25.6731s', '+30d42m36.450s')
box1_coord2 = SkyCoord('1h33m31.2569s', '+30d44m38.990s')
box1_coord3 = SkyCoord('1h33m25.7641s', '+30d38m03.325s')
box1_coord4 = SkyCoord('1h34m20.2174s', '+30d35m58.989s')

box2_coord1 = SkyCoord('1h34m38.6055s', '+30d48m51.161s')
box2_coord2 = SkyCoord('1h33m44.3862s', '+30d50m52.929s')
box2_coord3 = SkyCoord('1h33m38.7691s', '+30d44m20.574s')
box2_coord4 = SkyCoord('1h34m33.2784s', '+30d42m17.443s')

box3_coord1 = SkyCoord('1h34m22.1601s', '+30d35m57.480s')
box3_coord2 = SkyCoord('1h33m27.6702s', '+30d38m01.780s')
box3_coord3 = SkyCoord('1h33m22.3027s', '+30d31m23.962s')
box3_coord4 = SkyCoord('1h34m16.6213s', '+30d29m20.769s')

HST_box1 = Polygon(((box1_coord1.ra.degree, box1_coord1.dec.degree), (box1_coord2.ra.degree, box1_coord2.dec.degree), (box1_coord3.ra.degree, box1_coord3.dec.degree), (box1_coord4.ra.degree, box1_coord4.dec.degree)))
HST_box2 = Polygon(((box2_coord1.ra.degree, box2_coord1.dec.degree), (box2_coord2.ra.degree, box2_coord2.dec.degree), (box2_coord3.ra.degree, box2_coord3.dec.degree), (box2_coord4.ra.degree, box2_coord4.dec.degree)))
HST_box3 = Polygon(((box3_coord1.ra.degree, box3_coord1.dec.degree), (box3_coord2.ra.degree, box3_coord2.dec.degree), (box3_coord3.ra.degree, box3_coord3.dec.degree), (box3_coord4.ra.degree, box3_coord4.dec.degree)))

#CFHT points
x_coord = CFHT_coords.ra.degree 
y_coord = CFHT_coords.dec.degree 

from descartes import PolygonPatch
fig, ax = plt.subplots(1)
HST_box1 = PolygonPatch(HST_box1)
HST_box2 = PolygonPatch(HST_box2)
HST_box3 = PolygonPatch(HST_box3)
ax.add_patch(HST_box1)
ax.add_patch(HST_box2)
ax.add_patch(HST_box3)


#check if point is in the HST mask -- this is wrong
HST_region = np.zeros(len(CFHT_RA)) #0 if not in region, 1 if in region
for i in range(len(CFHT_RA)):
	point = Point(x_coord[i], y_coord[i])
	if (HST_box1.contains(point) == True) | (HST_box2.contains(point) == True) | (HST_box3.contains(point) == True):
		HST_region[i] = 1

plt.scatter(x_coord, y_coord, alpha=0.3, c=HST_region)
plt.xlim(23, 24)
plt.ylim(30, 31)
plt.colorbar()
plt.savefig('/Users/amandaquirk/Desktop/checking_HST_region.png')
plt.close()
print('done with HST region!')

'''
save the data =================================================================
'''
np.savetxt('/Users/amandaquirk/Desktop/CFHT_target_criteria.txt', np.c_[CFHT_RA, CFHT_Dec, CFHT_i, CFHT_g, isolated, mask, mag_cut, HST_region], fmt='%1.16f', delimiter=' ', header='RA (deg), Dec (deg), i (mag), g (mag), isolation criteria, mask, brightness criteria, in the HST region?')


