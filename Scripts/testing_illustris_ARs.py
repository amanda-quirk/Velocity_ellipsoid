import numpy as np 
import matplotlib.pyplot as plt 

#checking ARs
# star1_med_disp_R, star2_med_disp_R, star3_med_disp_R, star4_med_disp_R, star1_med_disp_Z, star2_med_disp_Z, star3_med_disp_Z, star4_med_disp_Z, star1_med_disp_phi, star2_med_disp_phi, star3_med_disp_phi, star4_med_disp_phi = np.loadtxt('/Volumes/Titan/analogs/data/star_particle_dispersion.txt', usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), unpack=True)

# star1_AR1 = star1_med_disp_phi / star1_med_disp_R
# star2_AR1 = star2_med_disp_phi / star2_med_disp_R
# star3_AR1 = star3_med_disp_phi / star3_med_disp_R
# star4_AR1 = star4_med_disp_phi / star4_med_disp_R

# star1_AR2 = star1_med_disp_Z / star1_med_disp_R
# star2_AR2 = star2_med_disp_Z / star2_med_disp_R
# star3_AR2 = star3_med_disp_Z / star3_med_disp_R
# star4_AR2 = star4_med_disp_Z / star4_med_disp_R

# plt.hist(star1_AR1[~np.isnan(star1_AR1)], bins=np.linspace(0.4, 1, 25), alpha=0.5, color='b', label='avg={}, std={}'.format(round(np.median(star1_AR1[~np.isnan(star1_AR1)]), 3), round(np.std(star1_AR1[~np.isnan(star1_AR1)]), 3))) 
# plt.hist(star2_AR1[~np.isnan(star2_AR1)], bins=np.linspace(0.4, 1, 25), alpha=0.5, color='m', label='avg={}, std={}'.format(round(np.median(star2_AR1[~np.isnan(star2_AR1)]), 3), round(np.std(star2_AR1[~np.isnan(star2_AR1)]), 3)))
# plt.hist(star3_AR1[~np.isnan(star3_AR1)], bins=np.linspace(0.4, 1, 25), alpha=0.5, color='g', label='avg={}, std={}'.format(round(np.median(star3_AR1[~np.isnan(star3_AR1)]), 3), round(np.std(star3_AR1[~np.isnan(star3_AR1)]), 3)))
# plt.hist(star4_AR1[~np.isnan(star4_AR1)], bins=np.linspace(0.4, 1, 25), alpha=0.5, color='r', label='avg={}, std={}'.format(round(np.median(star4_AR1[~np.isnan(star4_AR1)]), 3), round(np.std(star4_AR1[~np.isnan(star4_AR1)]), 3)))
# plt.xlabel(r'$\sigma{\phi} / \sigma{R}$')
# plt.xlim(0.39, 1.11)
# plt.legend()
# plt.savefig('/Users/amandaquirk/Desktop/AR1_M33.png')
# plt.close()

# plt.hist(star1_AR2[~np.isnan(star1_AR2)], bins=np.linspace(0.6, 1.5, 25), alpha=0.5, color='b', label='avg={}, std={}'.format(round(np.median(star1_AR2[~np.isnan(star1_AR2)]), 3), round(np.std(star1_AR2[~np.isnan(star1_AR2)]), 3)))
# plt.hist(star2_AR2[~np.isnan(star2_AR2)], bins=np.linspace(0.6, 1.5, 25), alpha=0.5, color='m', label='avg={}, std={}'.format(round(np.median(star2_AR2[~np.isnan(star2_AR2)]), 3), round(np.std(star2_AR2[~np.isnan(star2_AR2)]), 3)))
# plt.hist(star3_AR2[~np.isnan(star3_AR2)], bins=np.linspace(0.6, 1.5, 25), alpha=0.5, color='g', label='avg={}, std={}'.format(round(np.median(star3_AR2[~np.isnan(star3_AR2)]), 3), round(np.std(star3_AR2[~np.isnan(star3_AR2)]), 3)))
# plt.hist(star4_AR2[~np.isnan(star4_AR2)], bins=np.linspace(0.6, 1.5, 25), alpha=0.5, color='r', label='avg={}, std={}'.format(round(np.median(star4_AR2[~np.isnan(star4_AR2)]), 3), round(np.std(star4_AR2[~np.isnan(star4_AR2)]), 3)))
# plt.xlabel(r'$\sigma{Z} / \sigma{R}$')
# plt.xlim(0.45, 1.65)
# plt.legend()
# plt.savefig('/Users/amandaquirk/Desktop/AR2_M33.png')
# plt.close()

#checking velocity ellipsoid for individual halo ==================================================================================
AR2_value = 0.97

def sine(angle): #deg
	return np.sin(np.deg2rad(angle))

def cosine(angle): #deg
	return np.cos(np.deg2rad(angle)) 

def sigma_LOS(oR, oP, i, PA): #km/s, deg
	oLOS2 = (oR**2 * sine(PA)**2 + oP**2 * cosine(PA)**2) * sine(i)**2 + AR2_value**2 * oR**2 * cosine(i)**2
	return np.sqrt(oLOS2) #km/s

#data from one Illustris analog
x, y, oR_real, oLOS_real, oPhi_real = np.loadtxt('/Volumes/Titan/analogs/data/419510_star2_particle_dispersion.txt', unpack=True)

#caluclate position arguments
inc = np.zeros(len(x)) #since we rotated them to face on

#PA
def PA_calc(x, y):
	angle = np.arctan(y / x ) * 180 / np.pi #deg
	PA = [] #east of north
	for i in range(len(angle)):
		if y[i] < 0:
			PA.append(angle[i] + 180)
		else:
			PA.append(angle[i])
	return(PA)

PA = PA_calc(x, y)

oLOS_calc = sigma_LOS(oR_real, oPhi_real, inc, PA)

diff = oLOS_calc - oLOS_real

plt.hist(oLOS_calc, label='calc')
plt.hist(oLOS_real, alpha=0.4, label='real')
#plt.hist(diff)
plt.xlabel('sigma_LOS')
plt.legend()
plt.savefig('/Users/amandaquirk/Desktop/sigma_LOS_star2.png')
plt.close()

plt.hist(diff)
plt.xlabel('sigma_calc - sigma_LOS')
plt.savefig('/Users/amandaquirk/Desktop/sigma_LOS_diff_star2.png')
plt.close()

plt.hist(oR_real)
plt.xlabel('sigma_R')
plt.savefig('/Users/amandaquirk/Desktop/sigma_R_star2.png')
plt.close()

plt.hist(oPhi_real)
plt.xlabel('sigma_phi')
plt.savefig('/Users/amandaquirk/Desktop/sigma_phi_star2.png')
plt.close()




