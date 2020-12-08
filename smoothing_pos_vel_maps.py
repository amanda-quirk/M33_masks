'''
reads in all of the M33 data and outputs smoothed LOS velocities, velocity dispersion, and position velocity maps
'''

import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt 
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator
from matplotlib import patches

max_circle = 300
#read in catalogue ===========================================================================================================
data = np.genfromtxt('/Volumes/Titan/M33/Data/M33_2018b_phot_spec_CMD_sorted.txt', dtype=None, names='ID, ra, dec, F275W, F336W, mag1_ext, mag1_name, mag2_ext, mag2_name, F110W, F160W, z, vel, vel_aband, err, zqual, aband, time, mask, age_tag, HI, CO, Ha')

ID = data['ID']
ID = np.array([a.decode("utf-8") for a in ID])
ra = data['ra']
ra = np.array([a.decode("utf-8") for a in ra])
dec = data['dec']
dec = np.array([a.decode("utf-8") for a in dec])
vel = data['vel']
vel_aband = data['vel_aband']
err = data['err']
age = data['age_tag']
age = np.array([a.decode("utf-8") for a in age])
zqual = data['zqual']
HI = data['HI']
CO = data['CO']
Ha = data['Ha']
mask = data['mask']
mask = [a.decode("utf-8") for a in mask]
aband = data['aband']


fg = np.genfromtxt('/Volumes/Titan/M33/Data/M33_FG_IDs.txt', dtype=None)
fg_id = np.array([a.decode("utf-8") for a in fg])

#only look at things that have good quality =============================================================================
qual = ((zqual == 1) | (zqual > 2)) & (vel < 500) & (vel > -500) & (aband * const.c.to(u.km/u.s).value < 80) & (aband * const.c.to(u.km/u.s).value > -80)
in_m33 = np.invert(np.in1d(ID, fg_id))
good_qual = qual * in_m33
#Karrie also has some SN cuts
print('All data =', len(ID), '; Data that passed cut =', sum(good_qual))

#separate into age bins ======================================================================================================
MS = age == 'MS'
AGB = age == 'AGB' 
RGB = age == 'RGB'
HeB_all = (age == 'RHeB') | (age == 'BHeB') | (age == 'HeB')
young = MS | HeB_all #add weak CN to this later -- will pull from AGB stars


# #histograms of LOS v =========================================================================================================
plt.hist(vel[MS * good_qual], bins=range(-300, 100, 15), label=r'$\mu, \sigma$' + r'$={}\ , {}$'.format(round(np.median(vel[MS * good_qual])), round(np.std(vel[MS * good_qual]))) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, stacked=True, fill=False, color='b')
plt.legend()
plt.savefig('/Volumes/Titan/M33/Plots/MS_vel_hist.png')
plt.close()

plt.hist(vel[AGB * good_qual], bins=range(-300, 100, 15), label=r'$\mu, \sigma$' + r'$={}\ , {}$'.format(round(np.median(vel[AGB * good_qual])), round(np.std(vel[AGB * good_qual]))) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, stacked=True, fill=False, color='green')
plt.legend()
plt.savefig('/Volumes/Titan/M33/Plots/AGB_vel_hist.png')
plt.close()

plt.hist(vel[RGB * good_qual], bins=range(-300, 100, 15), label=r'$\mu, \sigma$' + r'$={}\ , {}$'.format(round(np.median(vel[RGB * good_qual])), round(np.std(vel[RGB * good_qual]))) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, stacked=True, fill=False, color='red')
plt.legend()
plt.savefig('/Volumes/Titan/M33/Plots/RGB_vel_hist.png')
plt.close()

plt.hist(vel[HeB_all * good_qual], bins=range(-300, 100, 15), label=r'$all\ \mu, \sigma$' + r'$={}\ , {}$'.format(round(np.median(vel[HeB_all * good_qual])), round(np.std(vel[HeB_all * good_qual]))) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, stacked=True, fill=False, color='black')
plt.legend()
plt.savefig('/Volumes/Titan/M33/Plots/HeB_vel_hist.png')
plt.close()

# #smoothing ===================================================================================================================
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

def smoothing(ids, zqual, ras, decs, errs, HI, CO, Ha, velocities, RGB=False):
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
	radius_goodcenter = []

	#remove stars that have unreliable velocities
	weight = calc_weights(errs) #error is already adjusted so can just calculate the weights
	sc = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle,u.deg))
	circle_radius = np.zeros(len(zqual))
	if RGB == False:
		for i in range(len(ras)):
			c1 = SkyCoord(ras[i], decs[i], unit=(u.hourangle,u.deg)) #go through all coordinates one at a time
			sep = c1.separation(sc)
			circleSize = 50 #arcseconds, starting size
			area = sep.arcsecond < circleSize #put stars into smoothing circle of this size
			while sum(area) < 15 and circleSize < max_circle: #let the circle grow until there are 15 stars in it but cut its growth
				circleSize += 5
				area = sep.arcsecond < circleSize
			circle_radius[i] = circleSize
			velocities_circ = velocities[area]
			weight_circ = weight[area]
			if sum(area) >= 15 and circleSize <= max_circle: #only want circles that aren't too big and have low statistics
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
				radius_goodcenter.append(circle_radius[i])
	else: #going to increase the number of neighbors and only look at stars that are likely disk stars
		for i in range(len(ras)):
			c1 = SkyCoord(ras[i], decs[i], unit=(u.hourangle,u.deg)) #go through all coordinates one at a time
			sep = c1.separation(sc)
			circleSize = 50 #arcseconds, starting size
			area = sep.arcsecond < circleSize #put stars into smoothing circle of this size
			while sum(area) < 30 and circleSize < 150: #let the circle grow until there are 15 stars in it but cut its growth
				circleSize += 5
				area = sep.arcsecond < circleSize
			circle_radius[i] = circleSize
			velocities_circ = velocities[area]
			weight_circ = weight[area]
			#add the disk and halo gaussin fitting here
			if sum(area) >= 30 and circleSize <= 150: #only want circles that aren't too big and have low statistics
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
				radius_goodcenter.append(circle_radius[i])
	return ra_goodcenter, dec_goodcenter, smoothed_v, smoothed_err, dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter, radius_goodcenter

#I am currently excluding the aband correction -- change to vel_aband if want to include it
MS_smoothed_data = smoothing(ID[MS * good_qual], zqual[MS * good_qual], ra[MS * good_qual], dec[MS * good_qual], err[MS * good_qual], HI[MS * good_qual], CO[MS * good_qual], Ha[MS * good_qual], vel[MS * good_qual])
print('done with MS smoothing: {1} / {0}'.format(sum(MS), len(MS_smoothed_data[0])))
AGB_smoothed_data = smoothing(ID[AGB * good_qual], zqual[AGB * good_qual], ra[AGB * good_qual], dec[AGB * good_qual], err[AGB * good_qual], HI[AGB * good_qual], CO[AGB * good_qual], Ha[AGB * good_qual], vel[AGB * good_qual])
print('done with AGB smoothing: {1} / {0}'.format(sum(AGB), len(AGB_smoothed_data[0])))
HeB_all_smoothed_data = smoothing(ID[HeB_all * good_qual], zqual[HeB_all * good_qual], ra[HeB_all * good_qual], dec[HeB_all * good_qual], err[HeB_all * good_qual], HI[HeB_all * good_qual], CO[HeB_all * good_qual], Ha[HeB_all * good_qual], vel[HeB_all * good_qual])
print('done with HeB_all smoothing: {1} / {0}'.format(sum(HeB_all), len(HeB_all_smoothed_data[0])))
RGB_smoothed_data = smoothing(ID[RGB * good_qual], zqual[RGB * good_qual], ra[RGB * good_qual], dec[RGB * good_qual], err[RGB * good_qual], HI[RGB * good_qual], CO[RGB * good_qual], Ha[RGB * good_qual], vel[RGB * good_qual])
print('done with RGB smoothing: {1} / {0}'.format(sum(RGB), len(RGB_smoothed_data[0])))
young_smoothed_data = smoothing(ID[young * good_qual], zqual[young * good_qual], ra[young * good_qual], dec[young * good_qual], err[young * good_qual], HI[young * good_qual], CO[young * good_qual], Ha[young * good_qual], vel[young * good_qual])
print('done with young smoothing: {1} / {0}'.format(sum(young), len(young_smoothed_data[0])))

#velocity position maps ======================================================================================================
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
	min_radius = min(circle_size) / 60 / 60 * scale_factor
	max_radius = max(circle_size) / 60 / 60 * scale_factor
	c0 = plt.Circle((-4, -3), min_radius, color='k', fill=False)
	c1 = plt.Circle((-4, -3), min_radius, color='k', fill=False)
	c2 = plt.Circle((-3, 2.5), max_radius, color='k', fill=False)
	c3 = plt.Circle((-3, 2.5), max_radius, color='k', fill=False)

	#set up ellipse for reference
	ylength = 22 / 3 #PA
	xlength = ylength * np.cos(54 * np.pi / 180) #inclination
	e0 = patches.Ellipse((0, 0), xlength, ylength, angle=180 - 22, linewidth=.5, fill=False, zorder=0)
	e1 = patches.Ellipse((0, 0), xlength, ylength, angle=180 - 22, linewidth=.5, fill=False, zorder=0)
	e2 = patches.Ellipse((0, 0), xlength, ylength, angle=180 - 22, linewidth=.5, fill=False, zorder=0)

	#plot the stuff -- first column = individual v_LOS; second = smoothed v_LOS; third = velocity dispersion
	f, axes = plt.subplots(1,3, sharey=True, sharex=True, figsize=(11.1, 14 / 4))

	axes[0].add_patch(e0)
	axes[1].add_patch(e1)
	axes[2].add_patch(e2)

	axes[1].add_artist(c0)
	axes[2].add_artist(c1)
	axes[1].add_artist(c2)
	axes[2].add_artist(c3)

	im0 = axes[0].scatter(xi_ind * scale_factor, eta_ind * scale_factor, c=ind_vel, cmap='plasma', s=6, vmin=-260,vmax=-100) 
	im1 = axes[1].scatter(xi_sm * scale_factor, eta_sm * scale_factor, c=vel_smoothed, cmap='plasma', s=6, vmin=-260,vmax=-100) 
	im2 = axes[2].scatter(xi_sm * scale_factor, eta_sm * scale_factor, c=dispersion, cmap='viridis', s=6, vmin=0,vmax=75) 

	for ax in axes:
		ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
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
		ax.axis('equal')

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
	axes[0].set_ylim(-6.8, 3.8)
	axes[1].set_ylim(-6.8, 3.8)
	axes[2].set_ylim(-6.8, 3.8)
	axes[0].set_xlim(6.5, -5)
	axes[1].set_xlim(6.5, -5)
	axes[2].set_xlim(6.5, -5)
	axes[1].set_title('{}'.format(age), fontsize=13)
	axes[0].annotate(r'$v_{LOS}$', xy=(3, 3), horizontalalignment='right', fontsize=12)
	axes[1].annotate(r'$\overline{v}_{LOS}$', xy=(3, 3), horizontalalignment='right', fontsize=12)
	axes[2].annotate(r'$\sigma$', xy=(3, 3), horizontalalignment='right', fontsize=12)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig('/Volumes/Titan/M33/Plots/M33_maps_{}.png'.format(age), bbox_inches='tight')
	plt.close()

position_map(ra[MS * good_qual], dec[MS * good_qual], vel[MS * good_qual], MS_smoothed_data[0], MS_smoothed_data[1], MS_smoothed_data[2], MS_smoothed_data[4], MS_smoothed_data[-1], 14.12, 'MS')
position_map(ra[AGB * good_qual], dec[AGB * good_qual], vel[AGB * good_qual], AGB_smoothed_data[0], AGB_smoothed_data[1], AGB_smoothed_data[2], AGB_smoothed_data[4], AGB_smoothed_data[-1], 14.12, 'AGB')
position_map(ra[HeB_all * good_qual], dec[HeB_all * good_qual], vel[HeB_all * good_qual], HeB_all_smoothed_data[0], HeB_all_smoothed_data[1], HeB_all_smoothed_data[2], HeB_all_smoothed_data[4], HeB_all_smoothed_data[-1], 14.12, 'HeB_all')
position_map(ra[RGB * good_qual], dec[RGB * good_qual], vel[RGB * good_qual], RGB_smoothed_data[0], RGB_smoothed_data[1], RGB_smoothed_data[2], RGB_smoothed_data[4], RGB_smoothed_data[-1], 14.12, 'RGB')
position_map(ra[young * good_qual], dec[young * good_qual], vel[young * good_qual], young_smoothed_data[0], young_smoothed_data[1], young_smoothed_data[2], young_smoothed_data[4], young_smoothed_data[-1], 14.12, 'Younger Populations')


# # #save the data into catalogues  ==============================================================================================
#I am being lazy and saving as separate files. might be easier in the future to save as one file but who knows
np.savetxt('/Volumes/Titan/M33/Data/M33_2018b_smoothed_kinematics_MS.txt', np.c_[MS_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter, smoothing circle radius')
np.savetxt('/Volumes/Titan/M33/Data/M33_2018b_smoothed_kinematics_AGB.txt', np.c_[AGB_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter, smoothing circle radius')
np.savetxt('/Volumes/Titan/M33/Data/M33_2018b_smoothed_kinematics_HeB_all.txt', np.c_[HeB_all_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter, smoothing circle radius')
np.savetxt('/Volumes/Titan/M33/Data/M33_2018b_smoothed_kinematics_RGB.txt', np.c_[RGB_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter, smoothing circle radius')
np.savetxt('/Volumes/Titan/M33/Data/M33_2018b_smoothed_kinematics_young.txt', np.c_[young_smoothed_data], fmt='%s', delimiter='\t', header='ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter, smoothing circle radius')

