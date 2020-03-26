import numpy as np 
import matplotlib.pyplot as plt 

def sine(angle):
	return np.sin(np.deg2rad(angle))

def cosine(angle):
	return np.cos(np.deg2rad(angle))

def sigma_LOS(oR, oZ, oP, i, PA): #i and PA come from tilted ring table (any dispersion values?)
	oLOS2 = (oR**2 * sine(PA)**2 + oP**2 * cosine(PA)**2) * sine(i)**2 + oZ**2 * cosine(i)**2
	return np.sqrt(oLOS2)

def va2_vavc(oR, oP, oZ, R, vc): #vc and R from titled ring table
	Rd = 5.76 #kpc
	terms = np.zeros(len(oR))
	for j in range(len(oR)):
		if R[j] < 10:
			k = 2
		else:
			k = 1
		terms[j] = -oR[j]**2 * ((oP[j]**2 / oR[j]**2) - 1.5 + (k * R[j] / Rd) + (oZ[j]**2 / (2 * oR[j]**2)))
	return terms

def va(oR, oP, oZ, R, vc):
	C = va2_vavc(oR, oP, oZ, R, vc)
	pos = vc + np.sqrt(vc**2 + C)
	neg = vc - np.sqrt(vc**2 + C)
	return pos, neg #want to use the negative root

# AR1 = .8
# AR2 = np.sqrt(.2)
# num_data_points = 100

# #going to make num_data_points random oR terms and get oZ and oP from that
# fake_oR = np.random.normal(80, 30, num_data_points)
# fake_oR = [abs(a) for a in fake_oR] #don't want any negative dispersions

# #make factor so that the axis ratio is obsurced
# AR1_wiggle = np.random.normal(.5, .5, num_data_points) * 2 
# AR1_wiggle = [abs(a) for a in AR1_wiggle]
# fake_oP = fake_oR * np.array(AR1_wiggle) * AR1

# AR2_wiggle = np.random.normal(.4, .4, num_data_points) * 3 
# AR2_wiggle = [abs(a) for a in AR2_wiggle]
# fake_oZ = fake_oR * np.array(AR2_wiggle) * AR2

# # plt.clf()
# # plt.scatter(fake_oR, fake_oZ, c='k')
# # plt.plot(np.linspace(min(fake_oR), max(fake_oR)), AR2 * np.linspace(min(fake_oR), max(fake_oR)), c='b')
# # plt.show() 

# # #assign each fake star a R, vc, and i based on HI tilted ring model 
# # all_rs, all_is, all_vcs = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', usecols=(0, 2, 3), unpack=True)
# # N = np.random.uniform(2, len(all_vcs) - 1, num_data_points)
# # N = [int(a) for a in N]

# # radii = all_rs[N]
# # inclinations = all_is[N]
# # v_circs = all_vcs[N]

# # #assign each fake star a PA at random
# # pos_angs = np.random.uniform(0, 360, num_data_points)

# # #calculate the oLOS and va term using above functions
# fake_oR = np.array(fake_oR)
# # o_LOS = sigma_LOS(fake_oR, fake_oZ, fake_oP, inclinations, pos_angs)
# # #va_term = va2_vavc(fake_oR, fake_oP, fake_oZ, radii, v_circs)
# # va_pos, va_neg = va(fake_oR, fake_oP, fake_oZ, radii, v_circs)

# #checking that there is a roughly linear relationship between oP**2 and va
# # plt.scatter(fake_oP**2, va_neg)
# # plt.xlim(0, 20000)
# # plt.ylim(-20, 140)
# # plt.show()

# #save data file -- has axis ratio and jean's
# #np.savetxt('/Users/amandaquirk/Documents/Ellipsoid/Data/fake_data.txt', np.c_[radii, inclinations, pos_angs, v_circs, fake_oR, fake_oP, fake_oZ, o_LOS, va_neg], header='r (kpc), inclination (deg), PA (deg), circular vel (km/s), sigmaR (km/s), sigmaPhi (km/s), sigmaZ (km/s), sigmaLOS (km/s), AD (km/s)')

# #creating fake data that doesn't rely on jean's but has axis ratios
# # PAs = np.random.normal(180, 180, num_data_points)
# # incs = np.random.normal(77, 10, num_data_points)
# # o_LOS = sigma_LOS(fake_oR, fake_oZ, fake_oP, incs, PAs)
# #np.savetxt('/Users/amandaquirk/Desktop/fake_data_nojeans_nowiggle.txt', np.c_[fake_oR, fake_oP, fake_oZ, o_LOS, incs, PAs], header='sigmaR (km/s), sigmaPhi (km/s), sigmaZ (km/s), sigmaLOS (km/s), inclination (deg), PA (deg)')

# #has no axis ratios
# fake_oR = np.random.normal(80, 30, num_data_points)
# fake_oZ = np.random.normal(110, 40, num_data_points)
# fake_oP = np.random.normal(50, 10, num_data_points)
# PAs = np.random.normal(180, 180, num_data_points)
# incs = np.random.normal(77, 10, num_data_points)
# o_LOS = sigma_LOS(fake_oR, fake_oZ, fake_oP, incs, PAs)
# np.savetxt('/Users/amandaquirk/Desktop/fake_data_nojeans_noARs.txt', np.c_[fake_oR, fake_oP, fake_oZ, o_LOS, incs, PAs], header='sigmaR (km/s), sigmaPhi (km/s), sigmaZ (km/s), sigmaLOS (km/s), inclination (deg), PA (deg)')

#making a fake data set based off of a combo of data from Illustris and M31 data =================================================
import random
xi, eta, r, oLOS = np.loadtxt('../Data/M31_RG.txt', usecols=(0,1,2,5), unpack=True)

#randomly grab 250 points so that my fake data set is no super large
size = 250
indices = np.linspace(0, len(xi) - 1, len(xi))
ind = random.sample(list(indices), size)
ind = [int(a) for a in ind]

xi = xi[ind]
eta = eta[ind]
r = r[ind]
oLOS = oLOS[ind]

oR, oPhi = np.loadtxt('/Volumes/Titan/analogs/data/419510_star2_particle_dispersion.txt', usecols=(2,4,), unpack=True)
indices_I = np.linspace(0, len(oR) - 1, len(oR))
ind_I = random.sample(list(indices_I), size)
ind_I = [int(a) for a in ind_I]

oR = oR[ind_I]
oPhi = oPhi[ind_I]

plt.hist(oLOS)
plt.xlabel('sigma_LOS')
plt.savefig('/Users/amandaquirk/Desktop/sigma_LOS_fakedata.png')
plt.close()

plt.hist(oR)
plt.xlabel('sigma_R')
plt.savefig('/Users/amandaquirk/Desktop/sigma_R_fakedata.png')
plt.close()

plt.hist(oPhi)
plt.xlabel('sigma_phi')
plt.savefig('/Users/amandaquirk/Desktop/sigma_phi_fakedata.png')
plt.close()

AR1 = oPhi / oR
plt.hist(AR1, label='{}'.format(round(np.median(AR1), 3)))
plt.xlabel('AR1')
plt.legend()
plt.savefig('/Users/amandaquirk/Desktop/AR1_fakedata.png')
plt.close()

np.savetxt('../Data/fakedata_M31_Illustris.txt', np.c_[xi, eta, r, oLOS, oR, oPhi], header='xi (kpc), eta (kpc), r (kpc), oLOS (km/s), oR (km/s), oPhi (km/s)')
