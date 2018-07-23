import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits

'''
compiling data files for Laurent
'''

#data
#xi (kpc), eta (kpc), average v(km/s), v err, var, n, HI main, HI close, ID, orginal index
MS_xi, MS_eta, MS_v, MS_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)
AGy_xi, AGy_eta, AGy_v, AGy_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)
AGo_xi, AGo_eta, AGo_v, AGo_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)
RG_xi, RG_eta, RG_v, RG_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,2,4,), unpack=True)

AGy_xi_additional, AGy_eta_additional, AGy_v_additional, AGy_var_additional=np.loadtxt('/Users/amandaquirk/Documents/Ellipsoid/Data/AGy_additional_data.txt', unpack=True)
AGo_xi_additional, AGo_eta_additional, AGo_v_additional, AGo_var_additional=np.loadtxt('/Users/amandaquirk/Documents/Ellipsoid/Data/AGo_additional_data.txt', unpack=True)
RG_xi_additional, RG_eta_additional, RG_v_additional, RG_var_additional=np.loadtxt('/Users/amandaquirk/Documents/Ellipsoid/Data/RG_additional_data.txt', unpack=True)


#r (kpc)
MS_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_master_vrot.txt', usecols=(0,), unpack=True)
AGy_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_master_vrot.txt', usecols=(0,), unpack=True)
AGo_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_master_vrot.txt', usecols=(0,), unpack=True)
RG_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_master_vrot.txt', usecols=(0,), unpack=True)

#converts xi (kpc) and eta (kpc) to x in the shifted coordinates frame
def x_kpc(xi, eta):
	xi_deg=[a/13.67 for a in xi]
	eta_deg=[a/13.67 for a in eta]
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	x=[(a*cosine)-(b*sine) for a,b in zip(xi_deg, eta_deg)]
	return x 

#converts xi (kpc) and eta (kpc) to y in the shifted coordinates frame
def y_kpc(xi, eta):
	xi_deg=[a/13.67 for a in xi]
	eta_deg=[a/13.67 for a in eta]
	sine=np.sin(float(37*np.pi) / 180)
	cosine=np.cos(float(37*np.pi) / 180)
	y=[(a*cosine)+(b*sine) for a,b in zip(eta_deg, xi_deg)]
	return y

def distance(x, y):
	inclination_factor=np.cos(float(77*np.pi) / 180)**2
	ang_distance_sq=[(a**2)+(float(b**2)/inclination_factor) for a,b in zip(y,x)]
	ang_dist=[np.sqrt(a) for a in ang_distance_sq]
	dist=[a * 13.67 for a in ang_dist]
	return dist

AGy_x_additional = x_kpc(AGy_xi_additional, AGy_eta_additional)
AGy_y_additional = y_kpc(AGy_xi_additional, AGy_eta_additional)
AGo_x_additional = x_kpc(AGo_xi_additional, AGo_eta_additional)
AGo_y_additional = y_kpc(AGo_xi_additional, AGo_eta_additional)
RG_x_additional = x_kpc(RG_xi_additional, RG_eta_additional)
RG_y_additional = y_kpc(RG_xi_additional, RG_eta_additional)

AGy_r_additional = distance(AGy_x_additional, AGy_y_additional)
AGo_r_additional = distance(AGo_x_additional, AGo_y_additional)
RG_r_additional = distance(RG_x_additional, RG_y_additional)

#adding four age groups into two
young_xi = np.concatenate((np.concatenate((MS_xi, AGy_xi)), AGy_xi_additional))
young_eta = np.concatenate((np.concatenate((MS_eta, AGy_eta)), AGy_eta_additional))
young_v = np.concatenate((np.concatenate((MS_v, AGy_v)), AGy_v_additional))
young_var = np.concatenate((np.concatenate((MS_var, AGy_var)), AGy_var_additional))
young_r = np.concatenate((np.concatenate((MS_r, AGy_r)), AGy_r_additional))

old_xi = np.concatenate((np.concatenate((np.concatenate((RG_xi, AGo_xi)), AGo_xi_additional)), RG_xi_additional))
old_eta = np.concatenate((np.concatenate((np.concatenate((RG_eta, AGo_eta)), AGo_eta_additional)), RG_eta_additional))
old_v = np.concatenate((np.concatenate((np.concatenate((RG_v, AGo_v)), AGo_v_additional)), RG_v_additional))
old_var = np.concatenate((np.concatenate((np.concatenate((RG_var, AGo_var)), AGo_var_additional)), RG_var_additional))
old_r = np.concatenate((np.concatenate((np.concatenate((RG_r, AGo_r)), AGo_r_additional)), RG_r_additional))


#creating arrays for the age tag
young_tag = np.zeros(len(young_r))
old_tag = np.ones(len(old_r))

#getting PA
young_x=x_kpc(young_xi, young_eta)
old_x=x_kpc(old_xi, old_eta)
young_y=y_kpc(young_xi, young_eta)
old_y=y_kpc(old_xi, old_eta)

young_PA=[]
for i in range(len(young_x)):
	if young_y[i] > 0 and young_x[i]>0:
		rad=np.arctan(float(young_y[i])/young_x[i])
		deg=90-(float(rad*180)/np.pi)
	elif young_y[i] > 0 and young_x[i] < 0:
		rad=np.arctan(float(young_y[i])/young_x[i])
		deg=270-(float(rad*180)/np.pi)
	elif young_y[i] < 0 and young_x[i] < 0:
		rad=np.arctan(float(young_y[i])/young_x[i])
		deg=270-(float(rad*180)/np.pi)
	elif young_y[i] < 0 and young_x[i] > 0:
		rad=np.arctan(float(young_y[i])/young_x[i])
		deg=90-(float(rad*180)/np.pi)
	young_PA.append(deg)

old_PA=[]
for i in range(len(old_x)):
	if old_y[i] > 0 and old_x[i]>0:
		rad=np.arctan(float(old_y[i])/old_x[i])
		deg=90-(float(rad*180)/np.pi)
	elif old_y[i] > 0 and old_x[i] <0:
		rad=np.arctan(float(old_y[i])/old_x[i])
		deg=270-(float(rad*180)/np.pi)
	elif old_y[i] < 0 and old_x[i] < 0:
		rad=np.arctan(float(old_y[i])/old_x[i])
		deg=270-(float(rad*180)/np.pi)
	elif old_y[i] < 0 and old_x[i] > 0:
		rad=np.arctan(float(old_y[i])/old_x[i])
		deg=90-(float(rad*180)/np.pi)
	old_PA.append(deg)	

young_PA=np.array(young_PA)
old_PA = np.array(old_PA)

#plotting for sanity check
# print(len(young_PA), len(young_x))
# print(len(old_PA), len(old_x))

# plt.scatter(young_x, young_y, c=young_PA)
# plt.scatter(old_x, old_y, c=old_PA)
# plt.xlim(1.5,-1.5)
# plt.ylim(-1.5,1.5)
# plt.plot([0,0],[-11,16], c='k')
# plt.plot([-16, 6], [0,0], c='k')
# plt.colorbar()
# plt.show()

#only keeping data in the radial bin 0-30

keep_young = young_r > 0 
keep_old = old_r > 0

young_r = young_r[keep_young]
old_r = old_r[keep_old]
young_var = young_var[keep_young]
old_var = old_var[keep_old]
young_v = young_v[keep_young]
old_v = old_v[keep_old]
young_PA = young_PA[keep_young]
old_PA = old_PA[keep_old]
young_tag = young_tag[keep_young]
old_tag = old_tag[keep_old]

keep_young = young_r <= 30
keep_old = old_r <= 30

young_r = young_r[keep_young]
old_r = old_r[keep_old]
young_var = young_var[keep_young]
old_var = old_var[keep_old]
young_v = young_v[keep_young]
old_v = old_v[keep_old]
young_PA = young_PA[keep_young]
old_PA = old_PA[keep_old]
young_tag = young_tag[keep_young]
old_tag = old_tag[keep_old]

#figuring out which radial bin each star belongs to 
#parameters
Rmin = 0 #kpc
Rmax = 31
delta_r= 0.9 #kpc

#making radial bins
r_bins=np.linspace(Rmin, Rmax, (Rmax - Rmin) / delta_r + 1)

#will contain the number of the bin a star belongs to 
young_bin=np.zeros(len(young_r)) #number of the bin it's in
old_bin=np.zeros(len(old_r))
young_rmid=np.zeros(len(young_r)) #middle radius of bin
old_rmid=np.zeros(len(old_r))

#dividing each age groups into radial bin
for i in range(len(young_r)):
	for j in range(len(r_bins)-1):
		if young_r[i]>=r_bins[j] and young_r[i]<r_bins[j+1]:
			young_bin[i]=j
			young_rmid[i] = r_bins[j] + r_bins[j+1] / 2

for i in range(len(old_r)):
	for j in range(len(r_bins)-1):
		if old_r[i]>=r_bins[j] and old_r[i]<r_bins[j+1]:
			old_bin[i]=j
			old_rmid[i] = 0.5 * (r_bins[j] + r_bins[j+1])

#averaging dispersions in each bin
#will contain number of stars in each bin
young_bin_disp=np.zeros(len(r_bins))
old_bin_disp=np.zeros(len(r_bins))

#dividing each age groups into radial bin
for i in range(len(r_bins)-1):
	young_stars=[b for a,b in zip(young_r, young_var) if a>=r_bins[i] and a<r_bins[i+1]]
	young_bin_disp[i]=np.mean(young_stars)
	old_stars=[b for a,b in zip(old_r, old_var) if a>=r_bins[i] and a<r_bins[i+1]]
	old_bin_disp[i]=np.mean(old_stars)

#adding data together
Rbin=np.concatenate((young_bin, old_bin))
Rbin_mid=np.concatenate((young_rmid, old_rmid))
Radius=np.concatenate((young_r, old_r))
v_LOS=np.concatenate((young_v, old_v))
dispersion=np.concatenate((young_var, old_var))
PA=np.concatenate((young_PA, old_PA))
age=np.concatenate((young_tag, old_tag))

#for averaged file
avg_young_tag= np.zeros(len(r_bins)-1)
avg_old_tag= np.ones(len(r_bins)-1)
avg_age=np.concatenate((avg_young_tag, avg_old_tag))
Rbin_disp= np.concatenate((young_bin_disp[:-1], old_bin_disp[:-1]))
lower_bins= r_bins[:-1]

#writing data to files
c1=fits.Column(name='Rbin_number', array=Rbin, format='K')
c2=fits.Column(name='mid_rbin_kpc', array=Rbin_mid, format='D')
c3=fits.Column(name='Radius_kpc', array=Radius, format='D')
c4=fits.Column(name='v_LOS_kms', array=v_LOS, format='D')
c5=fits.Column(name='dispersion_kms', array=dispersion, format='D')
c6=fits.Column(name='PA_deg', array=PA, format='D')
c7=fits.Column(name='Age', array=age, format='K')
t= fits.BinTableHDU.from_columns([c1, c2, c3,c4,c5,c6,c7])
t.writeto('/Users/amandaquirk/Desktop/smoothing_circle_masterfile_quirk.fits')

# c1=fits.Column(name='lower_rbin', array=np.concatenate((lower_bins, lower_bins)), format='D')
# c2=fits.Column(name='avg_dispersion_kms', array=Rbin_disp, format='D')
# c3=fits.Column(name='Age', array=avg_age, format='K')
# t= fits.BinTableHDU.from_columns([c1, c2, c3])
# t.writeto('/Users/amandaquirk/Desktop/smoothing_circle_rbin_quirk.fits')

