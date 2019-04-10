from astropy.io import fits
from astropy.table import Table
from astropy import units as u 
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import rc 

def import_data(maskname):
	hdu = fits.open('/Users/amandaquirk/Desktop/old_masks/zspec.{}.fits'.format(maskname))
	data = Table(hdu[1].data)
	ra = data['RA']
	dec = data['DEC']
	z = data['Z']
	error = data['Z_ERR']
	zqual = data['ZQUALITY']
	i_mag = data['IMAG']
	r_mag = data['RMAG']
	aband = data['ABAND']
	ID = data['OBJNAME']

	c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
	ra = c.ra.value 
	dec = c.dec.value 

	v = z * 299792.458 #km/s

	color = i_mag - r_mag

	good_data = ((zqual > 2) | (zqual == 1)) & (error > 0) & (error < 100) & (ra > 23) & (ra < 24) & (dec > 30) & (dec < 31)

	return ra[good_data], dec[good_data], v[good_data], error[good_data], zqual[good_data], i_mag[good_data], r_mag[good_data], aband[good_data], ID[good_data], color[good_data] 

data1 = import_data('M33D2A')
data2 = import_data('M33D2B')
data3 = import_data('M33D3A')
data4 = import_data('M33D3B')
data5 = import_data('M33D3D')
data6 = import_data('M33D4A')
data7 = import_data('M33D4B')

rc('font', family = 'serif')
fig, ax=plt.subplots(1)
for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
ax.tick_params(axis='x',which='both',top='on', direction='in')
ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
ax.tick_params(axis='y',which='both',right='on', direction='in')
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
plt.tick_params(labelsize=12) 
plt.minorticks_on()
# plt.scatter(data1[0], data1[1], c=data1[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data2[0], data2[1], c=data2[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data3[0], data3[1], c=data3[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data4[0], data4[1], c=data4[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data5[0], data5[1], c=data5[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data6[0], data6[1], c=data6[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data7[0], data7[1], c=data7[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# clb = plt.colorbar()
# clb.set_label(r'$\rm Velocity\ (km\ s^{-1})$', fontsize = 12)
# plt.xlabel(r'$\rm RA\ (deg)$', fontsize=12)
# plt.ylabel(r'$\rm Dec\ (deg)$', fontsize=12)
# plt.gca().invert_xaxis()
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/velocity_maps.pdf', bbox_inches='tight')

all_vs = data1[2].tolist() + data2[2].tolist() + data3[2].tolist() + data4[2].tolist() + data5[2].tolist() + data6[2].tolist() + data7[2].tolist()

plt.hist(all_vs, bins=range(-300, 0, 8), stacked=True, color='b')
plt.xlabel(r'$\rm Velocity\ (km\ s^{-1})$', fontsize=12)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/velocity_hist.pdf', bbox_inches='tight')

# plt.scatter(data1[5], data1[9], c=data1[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data2[5], data2[9], c=data2[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data3[5], data3[9], c=data3[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data4[5], data4[9], c=data4[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data5[5], data5[9], c=data5[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data6[5], data6[9], c=data6[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# plt.scatter(data7[5], data7[9], c=data7[2], cmap='plasma', s=15, vmin=-250,vmax=0)#, alpha=0.7)
# clb = plt.colorbar()
# clb.set_label(r'$\rm Velocity\ (km\ s^{-1}$', fontsize = 12)
# plt.xlabel(r'$i$', fontsize=12)
# plt.ylabel(r'$i-r$', fontsize=12)
# plt.gca().invert_xaxis()
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/CMD.pdf', bbox_inches='tight')



