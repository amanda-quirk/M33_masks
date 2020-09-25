import numpy as np 
import matplotlib.pyplot as plt 
import sfdmap 
from astropy.coordinates import SkyCoord
from astropy import units as u
#from deprojecting_rotation_funcs import deprojection_geo
from matplotlib import rc
import h5py
# from shapely.geometry import Point
# from shapely.geometry.polygon import Polygon

#read in the PAndAS catalogue =======================================================================================================
ra, dec, g_mag, dg, flag1, i_mag, di, flag2 = np.loadtxt('../Data/PAndAS_.8deg.tsv', unpack = True)
isolation_file = h5py.File('/Volumes/Titan/M33/Data/PAndAS_isolation_tag.hdf5', 'r') 
isolation_tag = isolation_file["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated 
print('Read in all the data')

#apply extinction corrections =======================================================================================================
sdss_g = 3.303
sdss_i = 1.698

emap = sfdmap.SFDMap('../Data/sfddata-master/', scaling=1.0)
coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
extinct_map = emap.ebv(coords)

g_ext = g_mag - extinct_map * sdss_g
i_ext = i_mag - extinct_map * sdss_i
print('Did the extinction corrections')

#calcualte deprojected geo ==========================================================================================================
star = (flag1 == -1) & (flag2 == -1) #let's look at things that are most likely to be stars
isolated = isolation_tag == 0
#xi, eta, alpha, beta, dist, PA, theta, assigned_PA, assigned_i = deprojection_geo(ra, dec, 'M33', unit='deg')

#calculate color and minimum brightnesses
color = g_ext - i_ext
bright = i_mag < 22 #what we'll target for spectroscopy
minimum_bright = i_mag < 24

# #calculate bins and make CMD =======================================================================================================
# MS = (color < -0.5) & (i_ext > 20)
# HeB = (color > 0.75) & (color < 1.3) & (i_ext < 21.5)
# AGB = (i_ext < 20.75) & (i_ext > 20) & (color > 1.6)
# RGB = (color > 0.5) & (i_ext > 21) 
# # RGB_poly = Polygon([(1.4, 21), (2.4, 21.25), (3.55, 21.9), (3, 21.89), (2.6, 21.95), (1.9, 22.5), (1.8, 22.8), (1.5, 23.4), (1.2, 24.9), (.43, 23.25), (1.1, 22)])
# # HeB_poly = Polygon([(.8, 22), (1.1, 19.5), (1.4, 18), (1.8, 19), (1.6, 20), (1.1, 22)])
# # AGB_poly = Polygon([(1.4, 21), (2.4, 21.25), (3.55, 21.9), (3.55, 20), (1.8, 20)])
# # RGB_flag = np.zeros(len(i_ext))
# # HeB_flag = np.zeros(len(i_ext))
# # AGB_flag = np.zeros(len(i_ext))
# # for i in range(len(i_ext)):
# # 	point = Point(color[i], i_ext[i])
# # 	if RGB_poly.contains(point) == True:
# # 		RGB_flag[i] = 1
# # 	elif HeB_poly.contains(point) == True:
# # 		HeB_flag[i] = 1
# # 	elif AGB_poly.contains(point) == True:
# # 		AGB_flag[i] = 1

# # MS = (color < -0.5) & (i_ext > 20)
# # HeB = HeB_flag == 1
# # AGB = AGB_flag == 1
# # RGB = RGB_flag == 1

# def plot_CMD(distance_range, label, name):
# 	#population statistics
# 	N = sum(minimum_bright * star * distance_range * isolated) #total
# 	N_bright = sum(bright * star * distance_range * isolated) #total
# 	frac_MS =  sum(minimum_bright * star * distance_range * isolated * MS) / N 
# 	frac_HeB = sum(minimum_bright * star * distance_range * isolated * HeB) / N 
# 	frac_AGB = sum(minimum_bright * star * distance_range * isolated * AGB) / N 
# 	frac_RGB = sum(minimum_bright * star * distance_range * isolated * RGB) / N 

# 	frac_MS_bright =  sum(bright * star * distance_range * isolated * MS) / N_bright 
# 	frac_HeB_bright = sum(bright * star * distance_range * isolated * HeB) / N_bright 
# 	frac_AGB_bright = sum(bright * star * distance_range * isolated * AGB) / N_bright 
# 	frac_RGB_bright = sum(bright * star * distance_range * isolated * RGB) / N_bright 

# 	#plotting
# 	rc('font', family = 'serif')
# 	fig, ax=plt.subplots(1)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	plt.tick_params(which='both', width=1)
# 	plt.tick_params(which='major', length=7)
# 	plt.tick_params(which='minor', length=4)
# 	plt.tick_params(labelsize=12) 
# 	plt.minorticks_on()

# 	plt.scatter(color[star * distance_range], i_ext[star * distance_range], alpha = 0.2, s=.4)
# 	plt.xlabel('g - i')
# 	plt.ylabel('i')
# 	plt.plot([-2, 5], [22, 22], c='grey', linestyle='--')
# 	plt.plot([-2, 5], [24, 24], c='black', linestyle='--')
# 	ax.annotate('i < 24, i < 22', xy=(3.8, 25.5), horizontalalignment='right')
# 	ax.annotate('MS = {}, {}'.format(round(frac_MS, 2), round(frac_MS_bright, 2)),    xy=(3.8, 26), horizontalalignment='right')
# 	ax.annotate('HeB = {}, {}'.format(round(frac_HeB, 2), round(frac_HeB_bright, 2)), xy=(3.8, 26.5), horizontalalignment='right')
# 	ax.annotate('AGB = {}, {}'.format(round(frac_AGB, 2), round(frac_AGB_bright, 2)), xy=(3.8, 27), horizontalalignment='right')
# 	ax.annotate('RGB = {}, {}'.format(round(frac_RGB, 2), round(frac_RGB_bright, 2)), xy=(3.8, 27.5), horizontalalignment='right')
# 	plt.xlim(-1, 4)
# 	plt.ylim(28, 18)
# 	plt.title('{}'.format(label))
# 	plt.savefig('/Users/amandaquirk/Desktop/PAndAS_{}_isolated.png'.format(name))
# 	plt.close()

# distance_range1 = (dist > 7) & (dist < 8)
# distance_range2 = (dist > 8) & (dist < 9)
# distance_range3 = (dist > 9) & (dist < 10)
# distance_range4 = (dist > 12) & (dist < 14)
# far = (dist > 14)

# def calc_num_density(rout, rin, distance_range):
# 	#population statistics
# 	N = sum(minimum_bright * star * distance_range * isolated) #total
# 	N_bright = sum(bright * star * distance_range * isolated) #total
# 	num_MS =  sum(minimum_bright * star * distance_range * isolated * MS) 
# 	num_HeB = sum(minimum_bright * star * distance_range * isolated * HeB)
# 	num_AGB = sum(minimum_bright * star * distance_range * isolated * AGB)
# 	num_RGB = sum(minimum_bright * star * distance_range * isolated * RGB)

# 	num_MS_bright =  sum(bright * star * distance_range * isolated * MS) 
# 	num_HeB_bright = sum(bright * star * distance_range * isolated * HeB) 
# 	num_AGB_bright = sum(bright * star * distance_range * isolated * AGB) 
# 	num_RGB_bright = sum(bright * star * distance_range * isolated * RGB) 

# 	area = np.pi * (rout**2 - rin**2)

# 	return num_MS / area, num_HeB / area, num_AGB / area, num_RGB / area, num_MS_bright / area, num_HeB_bright / area, num_AGB_bright / area, num_RGB_bright / area

# plot_CMD(distance_range1, '7 to 8 kpc', '7-8')
# plot_CMD(distance_range2, '8 to 9 kpc', '8-9')
# plot_CMD(distance_range3, '9 to 10 kpc', '9-10')
# plot_CMD(distance_range4, '12 to 14 kpc', '12-14')
# plot_CMD(far, '14 to 20 kpc', '14-20')

#doing a comparison between the OG CFHT catalogue and this real PAndAS catalogue ====================================================
# from astropy.io import fits 

# h = fits. open('../Data/M33.GI.matchcat.fits')
# data = h[1].data 
# ra0 = data['ra1']
# dec0 = data['dec1']
# mag1 = data['MAG1_AUTO'] #g 
# mag2 = data['MAG2_AUTO'] #i 

# #apply extinction corrections
# coords0 = SkyCoord(ra=ra0, dec=dec0, unit=(u.deg, u.deg))
# extinct_map0 = emap.ebv(coords0)
# mag1_ext = mag1 - extinct_map0 * sdss_g
# mag2_ext = mag2 - extinct_map0 * sdss_i
# color0 = mag1_ext - mag2_ext

# #roughly cutting the area that we're looking at
# minor_axis = (ra < 23.5) & (ra > 23.15) & (dec > 30.2) & (dec < 30.5)
# disk_16 = (ra < 23.6) & (ra > 23.3) & (dec > 30.4) & (dec < 30.6)
# disk_18 = (ra < 23.6) & (ra > 23.3) & (dec > 30.6) & (dec < 30.9)
# minor_axis0 = (ra0 < 23.5) & (ra0 > 23.15) & (dec0 > 30.2) & (dec0 < 30.5)
# disk_16_0 = (ra0 < 23.6) & (ra0 > 23.3) & (dec0 > 30.4) & (dec0 < 30.6)
# disk_18_0 = (ra0 < 23.6) & (ra0 > 23.3) & (dec0 > 30.6) & (dec0 < 30.9)

# #examing the CMDs and position maps
# def plot_comp_CMD_map(mask, mask0, name):
# 	rc('font', family = 'serif')
# 	fig, ax=plt.subplots(1, 2, sharey=True, sharex=True)
# 	ax[0].scatter(color[mask], i_ext[mask], alpha = 0.4, s=.8, label='PAndAS', c='b')
# 	ax[1].scatter(color0[mask0], mag2_ext[mask0], alpha=0.2, s=.8, label='original catalogue', c='r')
# 	plt.xlabel('g - i')
# 	plt.ylabel('i')
# 	plt.xlim(-1, 4)
# 	plt.ylim(28, 18)
# 	plt.legend()
# 	plt.savefig('/Users/amandaquirk/Desktop/CFHT_comp_CMD_{}_all.png'.format(name))
# 	plt.close()

# 	rc('font', family = 'serif')
# 	fig, ax=plt.subplots(1, 2, sharey=True, sharex=True)
# 	ax[0].scatter(ra[mask], dec[mask], alpha = 0.4, s=.8, label='PAndAS', c='b')
# 	ax[1].scatter(ra0[mask0], dec0[mask0], alpha=0.2, s=.8, label='original catalogue', c='r')
# 	plt.gca().invert_xaxis()
# 	plt.xlabel('RA (deg)')
# 	plt.ylabel('Dec (deg)')
# 	plt.legend()
# 	plt.savefig('/Users/amandaquirk/Desktop/CFHT_comp_map_{}_all.png'.format(name))
# 	plt.close()

# plot_comp_CMD_map(minor_axis, minor_axis0, 'minor_axis')
# plot_comp_CMD_map(disk_16, disk_16_0, '2016_disk')
# plot_comp_CMD_map(disk_18, disk_18_0, '2018_disk')

#actually start the target list stuff!!! ============================================================================================
required = minimum_bright * isolated * star #want things to be bright and isolated and likely point sources, since the catalogue is big enough
coords = coords[required]
i_mag = i_mag[required] #not using extincton corrected values for alignment stars
g_mag = g_mag[required] #not using extincton corrected values for alignment stars
i_ext = i_ext[required]
g_ext = g_ext[required]
color = color[required]
list_num = np.ones(len(i_mag)) 
priority = np.zeros(len(i_mag)) 

#list assignments -- saving list 2 for duplicates if we have repeat positions
for i in range(len(i_mag)):
	if ((i_mag[i] > 15) & (i_mag[i] < 18.5)) | ((g_mag[i] > 15) & (g_mag[i] < 18.5)): #alignment/guide stars
		list_num[i] = 0
	elif color[i] < 0.6: #want to avoid blue stars
		list_num[i] = 3
	else:
		list_num[i] = 1

#priorities
for i in range(len(i_ext)):
	if list_num[i] == 0: #guide/alignment stars
		priority[i] = -2
	elif list_num[i] == 3:
		priority[i] = 0 #don't want any blue stars on our masks
	elif (i_ext[i] > 21) & (i_ext[i] < 22): #really what we want
		priority[i] = 100
	elif (i_ext[i] > 20.5) & (i_ext[i] < 21):
		priority[i] = 55
	elif (i_ext[i] > 22) & (i_ext[i] < 22.5):
		priority[i] = 30
	else:
		priority[i] = 1

print('Sorted by lists and priorities')
#exta info for dsimulator <3 
ID = np.arange(0, len(di))[required] #index in orginal file!!
formated_coords = coords.to_string('hmsdms', alwayssign=False) #will need to edit by hand to change hms dms to : until I find a better way
JD = np.zeros(len(i_ext)) + 2000
filter_tag = ['I' for x in range(0, len(i_ext))]
selection_tag = np.zeros(len(i_ext))

print('Saving the data')
#save the data
np.savetxt('/Users/amandaquirk/Desktop/2020B_targestlist.in', np.c_[ID, formated_coords, JD, i_mag, filter_tag, priority, list_num, selection_tag], fmt="%-s", delimiter='\t', header='ID, coordinates, coordinate reference frame, magnitude, filter, priority, list assignment, selection flag')
np.savetxt('/Users/amandaquirk/Desktop/2020B_targestlist_check.in', np.c_[ID, formated_coords, i_ext, g_ext, priority, list_num], fmt="%-s", delimiter='\t', header='ID, coordinates, i, g, priority, list assignment')

#testing the target list ============================================================================================================
data = np.genfromtxt('/Users/amandaquirk/Desktop/2020B_targestlist_check.in', dtype=None, names='ID, ra, dec, i, g, priority, list')
i_mag = data['i']
g_mag = data['g']
priority = data['priority']
color = g_mag - i_mag

rc('font', family = 'serif')
fig, ax=plt.subplots(1)
for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1)
ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
ax.tick_params(axis='x',which='both',top='on', direction='in')
ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
ax.tick_params(axis='y',which='both',right='on', direction='in')
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
plt.tick_params(labelsize=12) 
plt.minorticks_on()

plt.scatter(color, i_mag, alpha = 0.2, s=.4, c=priority, vmin=-5, vmax=100)
plt.colorbar()
plt.xlabel('g - i')
plt.ylabel('i')
plt.xlim(-1, 4)
plt.ylim(24.2, 12)
plt.savefig('/Users/amandaquirk/Desktop/CMD_target_list.png')
plt.close()




