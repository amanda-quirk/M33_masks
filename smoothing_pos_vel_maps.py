'''
reads in all of the M33 data and outputs smoothed LOS velocities, velocity dispersion, and position velocity maps
'''

import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt 
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator
from matplotlib import patches

#read in catalogue ===================================================================================================================================
data = np.genfromtxt('/Users/amandaquirk/Documents/M33/Data/M33_2018b_phot_spec.txt', dtype=None, names='ID, ra, dec, F275W, F336W, F475W, F814W, F110W, F160W, z, vel, err, zqual, aband, time, mask, age, HI, CO, Ha')

ID = data['ID']
ID = np.array([a.decode("utf-8") for a in ID])
ra = data['ra']
ra = np.array([a.decode("utf-8") for a in ra])
dec = data['dec']
dec = np.array([a.decode("utf-8") for a in dec])
vel = data['vel']
err = data['err']
age = data['age']
age = np.array([a.decode("utf-8") for a in age])
zqual = data['zqual']
HI = data['HI']
CO = data['CO']
Ha = data['Ha']
mask = data['mask']
mask = [a.decode("utf-8") for a in mask]

#only look at things that have high enough zqual =====================================================================================================
good_qual = (zqual == 1) | (zqual > 2)

ID = ID[good_qual]
ra = ra[good_qual]
dec = dec[good_qual]
vel = vel[good_qual]
err = err[good_qual]
age = age[good_qual]
zqual = zqual[good_qual]
HI = HI[good_qual]
CO = CO[good_qual]
Ha = Ha[good_qual]

#separate into age bins ==============================================================================================================================
MS = age == 'MS'
AGB = age == 'AGB'
RGB = (age == 'RGB') | ('R' in mask == True)
HeB = (age == 'RHeB') | (age == 'BHeB') | (age == 'HeB') #want to separate this out or lump into diff groups?

#smoothing ===========================================================================================================================================
#function to calculate the weights
def calc_weights(err):
        return 1 / (err**2)

def normed_weight(w):
        sum_weights=sum(w)
        return w / sum_weights

#function does the weighted meean
def weighted_mean(data,norm_w):
	return sum(data * norm_w)

#function does the weighted RMSE
def weighted_rmse(norm_w, data, mean):
	diff_sq = (data - mean)**2
	return np.sqrt(sum(diff_sq * norm_w))

def smoothing(ids, zqual, ras, decs, errs, HI, CO, Ha, velocities, circleSize):
	smoothed_v = []
	dispersion = []
	#below these values are not actually smoothed, just saving the ones for good centers
	ra_goodcenter = []
	dec_goodcenter = []
	smoothed_err = []
	HI_goodcenter = []
	CO_goodcenter = []
	Ha_goodcenter = []
	ID_goodcenter = []
	zqual_goodcenter = []

	#remove stars that have unreliable velocities
	reliable = (abs(velocities) < 1000) & (errs > 0) & (errs < 999) #km/s
	ras = ras[reliable]
	decs = decs[reliable]
	errs = errs[reliable]
	HI = HI[reliable]
	CO = CO[reliable]
	Ha = Ha[reliable]
	velocities = velocities[reliable]
	ID = ids[reliable]
	zqual = zqual[reliable]

	weight = calc_weights(errs) #error is already adjusted so can just calculate the weights
	sc = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle,u.deg))
	for i in range(len(ras)):
		c1 = SkyCoord(ras[i], decs[i], unit=(u.hourangle,u.deg)) #go through all coordinates one at a time
		sep = c1.separation(sc)
		area = sep.arcsecond < circleSize #put stars into smoothing circle of this size
		velocities_circ = velocities[area]
		weight_circ = weight[area]
		if len(velocities_circ) > 15: #only want circles with at least 15 stars
			normed_weights = normed_weight(weight_circ)
			smoothed_v.append(weighted_mean(velocities_circ, normed_weights)) #average the velocites
			dispersion.append(weighted_rmse(normed_weights, velocities_circ, weighted_mean(velocities_circ, normed_weights)))
			ra_goodcenter.append(ras[i]) #ha
			dec_goodcenter.append(decs[i]) #deg
			smoothed_err.append(errs[i])
			HI_goodcenter.append(HI[i])
			CO_goodcenter.append(CO[i])
			Ha_goodcenter.append(Ha[i])
			ID_goodcenter.append(ids[i])
			zqual_goodcenter.append(zqual[i])
	return ra_goodcenter, dec_goodcenter, smoothed_v, smoothed_err, dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter

MS_smoothed_data = smoothing(ID[MS], zqual[MS], ra[MS], dec[MS], err[MS], HI[MS], CO[MS], Ha[MS], vel[MS], 300)
print('done with MS smoothing')
AGB_smoothed_data = smoothing(ID[AGB], zqual[AGB], ra[AGB], dec[AGB], err[AGB], HI[AGB], CO[AGB], Ha[AGB], vel[AGB], 300)
print('done with AG smoothing')
HeB_smoothed_data = smoothing(ID[HeB], zqual[HeB], ra[HeB], dec[HeB], err[HeB], HI[HeB], CO[HeB], Ha[HeB], vel[HeB], 300)
print('done with HeB smoothing')
RGB_smoothed_data = smoothing(ID[RGB], zqual[RGB], ra[RGB], dec[RGB], err[RGB], HI[RGB], CO[RGB], Ha[RGB], vel[RGB], 300)
print('done with RGB smoothing')

#velocity position maps ==============================================================================================================================
def position_map(ra, dec, ind_vel, ra_smoothed, dec_smoothed, vel_smoothed, dispersion, circle_size, scale_factor, age): #individual by age for now
	#convert to xi and eta centered on M33
	m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
	sc_ind = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
	c_inm33_ind = sc_ind.transform_to(m33.skyoffset_frame())
	xi_ind, eta_ind = c_inm33_ind.lon, c_inm33_ind.lat
	xi_ind = xi_ind.degree
	eta_ind = eta_ind.degree
	sc_sm = SkyCoord(ra=ra_smoothed, dec=dec_smoothed, unit=(u.hourangle,u.deg))
	c_inm33_sm = sc_sm.transform_to(m33.skyoffset_frame())
	xi_sm, eta_sm = c_inm33_sm.lon, c_inm33_sm.lat
	xi_sm = xi_sm.degree
	eta_sm = eta_sm.degree

	#set up the diagram of smoothing circle
	centerx = 10
	centery = 0
	radius = circle_size / 60 / 60 * scale_factor
	c0 = plt.Circle((centerx, centery), radius, color='k', fill=False)
	c1 = plt.Circle((centerx, centery), radius, color='k', fill=False)

	#set up ellipse for reference
	ylength = 22 / 3 #PA
	xlength = ylength * np.cos(54 * np.pi / 180) #inclination
	e0 = patches.Ellipse((0, 0), xlength, ylength, angle=180 - 22, linewidth=.5, fill=False, zorder=0)
	e1 = patches.Ellipse((0, 0), xlength, ylength, angle=180 - 22, linewidth=.5, fill=False, zorder=0)
	e2 = patches.Ellipse((0, 0), xlength, ylength, angle=180 - 22, linewidth=.5, fill=False, zorder=0)

	#plot the stuff -- first column = individual v_LOS; second = smoothed v_LOS; third = velocity dispersion
	f, axes = plt.subplots(1,3, sharey=True, sharex=False, figsize=(11.1,10.26 / 4))

	axes[0].add_patch(e0)
	axes[1].add_patch(e1)
	axes[2].add_patch(e2)

	axes[1].add_artist(c0)
	axes[2].add_artist(c1)

	im0 = axes[0].scatter(xi_ind * scale_factor, eta_ind * scale_factor, c=ind_vel, cmap='plasma', s=6, vmin=-260,vmax=-100) 
	im1 = axes[1].scatter(xi_sm * scale_factor, eta_sm * scale_factor, c=vel_smoothed, cmap='plasma', s=6, vmin=-260,vmax=-100) 
	im2 = axes[2].scatter(xi_sm * scale_factor, eta_sm * scale_factor, c=dispersion, cmap='viridis', s=6, vmin=0,vmax=90) 

	for ax in axes:
		ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
		#ax.set_ylim(-2.5,15.5)
		ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
		ax.tick_params(axis='y',which='both',right='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		nbins = 7
		ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)

	axes[0].set_ylabel(r'$\eta\ (kpc)$', fontsize=13)
	f.subplots_adjust(right=0.885)
	cbar_ax1 = f.add_axes([0.889, 0, 0.01, 1])
	cbar_ax2 = f.add_axes([0.990, 0, 0.01, 1])
	clb1 = f.colorbar(im0, cax=cbar_ax1)
	clb2 = f.colorbar(im2, cax=cbar_ax2)
	clb1.set_label(r'$\rm Individual,\ Mean\ LOS\ velocity:\ v, \ \overline{v}\ (km\ s^{-1})$', fontsize=11)
	clb2.set_label(r'$\rm Velocity\ Dispersion: \sigma\ (km\ s^{-1})$', fontsize=11)
	axes[0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1].set_title('{}'.format(age), fontsize=13)
	axes[0].annotate(r'$v_{LOS}$', xy=(-1.7, 3), horizontalalignment='right', fontsize=12)
	axes[1].annotate(r'$\overline{v}_{LOS}$', xy=(-1.7, 3), horizontalalignment='right', fontsize=12)
	axes[2].annotate(r'$\sigma$', xy=(-1.7, 3), horizontalalignment='right', fontsize=12)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig('/Users/amandaquirk/Desktop/M33_maps_{}.png'.format(age), bbox_inches='tight')

position_map(ra[MS], dec[MS], vel[MS], MS_smoothed_data[0], MS_smoothed_data[1], MS_smoothed_data[2], MS_smoothed_data[4], 300, 14.12, 'MS')
position_map(ra[AGB], dec[AGB], vel[AGB], AGB_smoothed_data[0], AGB_smoothed_data[1], AGB_smoothed_data[2], AGB_smoothed_data[4], 300, 14.12, 'AGB')
position_map(ra[HeB], dec[HeB], vel[HeB], HeB_smoothed_data[0], HeB_smoothed_data[1], HeB_smoothed_data[2], HeB_smoothed_data[4], 300, 14.12, 'HeB')
position_map(ra[RGB], dec[RGB], vel[RGB], RGB_smoothed_data[0], RGB_smoothed_data[1], RGB_smoothed_data[2], RGB_smoothed_data[4], 300, 14.12, 'RGB')

#save the data into catalogues  ======================================================================================================================
#I am being lazy and saving as separate files. might be easier in the future to save as one file but who knows
np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_smoothed_kinematics_MS.txt', np.c_[MS_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter')
np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_smoothed_kinematics_AGB.txt', np.c_[AGB_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter')
np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_smoothed_kinematics_HeB.txt', np.c_[HeB_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter')
np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_smoothed_kinematics_RGB.txt', np.c_[RGB_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter')

