from astropy.io import fits
import numpy as np 
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt 
from matplotlib import rc 

def import_data(maskname): #import spectroscopy data from zspec files
	hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	data = hdu[1].data
	zqual = data['ZQUALITY'] #eventually only want to use 1, 3, or 4
	ID = data['OBJNAME']
	return zqual, ID   

z_A, ID_A = import_data('A1M33P')
z_B, ID_B = import_data('B1M33P')
z_C, ID_C = import_data('C1M33P')
z_D, ID_D = import_data('D1M33P')
z_E, ID_E = import_data('E1M33P')

ref_data = np.genfromtxt('/Users/amandaquirk/Documents/M33/Data/all_target_list.in', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, magnitude2, filter2, priority, list assignment, selection flag, HST F110W, HST F160W, HST F275W, HST F336W')

ref_ID = ref_data['ID']
ref_ID = [a.decode("utf-8") for a in ref_ID]
mag1 = ref_data['magnitude1']
mag2 = ref_data['magnitude2']

zspec_IDs = list(ID_A) + list(ID_B) + list(ID_C) + list(ID_D) + list(ID_E)
zspec_zquals = list(z_A) + list(z_B) + list(z_C) + list(z_D) + list(z_E)

F475W = []
F814W = []
zquals = []
for i in range(len(zspec_IDs)):
	if zspec_IDs[i].startswith('CFHT') == False and zspec_IDs[i].startswith('serendip') == False:
		N = ref_ID.index(zspec_IDs[i])
		F475W.append(mag2[N])
		F814W.append(mag1[N])
		zquals.append(zspec_zquals[i])

def single_plot():
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
	return

color = np.array(F475W) - np.array(F814W)
zquals = np.array(zquals)
F814W = np.array(F814W)
good_data = (zquals == 1) | (zquals > 2)
bad_data = (zquals == 2) | (zquals < 0)

single_plot()
#cmaplist = ['dodgerblue', 'r', 'green', 'm']
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist)
#norm = matplotlib.colors.BoundaryNorm([0.5,1.5,2.5,3.5,4.5], cmap.N)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
plt.scatter(color[bad_data], F814W[bad_data], c='deeppink', alpha=.7, s=10)#c=zquals, cmap=cmap, norm=norm, edgecolor='none', vmin=1, s=25, alpha=.9)
plt.scatter(color[good_data], F814W[good_data], c='orange', edgecolors='peru', alpha=.6, s=10)
plt.xlim(-2 ,8)
plt.ylim(24.6, 13.5)
plt.xlabel('F475W-F814W')
plt.ylabel('F814W')
#clb = plt.colorbar(ticks=np.linspace(1,4,4))
#clb.set_label('Z Quality')
plt.savefig('/Users/amandaquirk/Desktop/CMD_talks.pdf', transparent=True)
plt.close()


