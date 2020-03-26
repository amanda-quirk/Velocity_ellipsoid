import emcee 
import numpy as np 
import matplotlib.pyplot as plt 
#from scipy.optimize import fsolve
import corner

'''
updating orginial_dispersion_MCMC.py to remove dependence on the jean's equation

This will be the MCMC procedure to (hopefully) decompose the velocity ellipsoid
'''
AR2_value = 1.02
galaxy_inc = 77
galaxy_PA = 37
galaxy_name = 'M31'

#data models ===================================================================================================================
def sine(angle): #deg
	return np.sin(np.deg2rad(angle))

def cosine(angle): #deg
	return np.cos(np.deg2rad(angle))

#can't constrain AR2 so not fitting for it
def sigma_LOS(oR, oP, i, PA): #km/s, deg
	oLOS2 = (oR**2 * sine(PA)**2 + oP**2 * cosine(PA)**2) * sine(i)**2 + AR2_value**2 * oR**2 * cosine(i)**2
	return np.sqrt(oLOS2) #km/s

#emcee functions ===============================================================================================================
#logliklihood function -- right now has no errors incorporated
def ll(theta, oLOS_real, i, PA):
	oR_guess, oPhi_guess = theta
	oLOS_calc = sigma_LOS(oR_guess, oPhi_guess, i, PA)
	oLOS_comp = -1 / 2 * (oLOS_real - oLOS_calc)**2
	return oLOS_comp 

#priors -- assuming normal distribution but play around with this and check it later
def lprior(theta):
	#pulling priors from theta
	oR_guess, oPhi_guess = theta 
	#mean and standard deviations of components -- will eventually want to change these for each age group
	if (oR_guess > 0 and oR_guess < 150) and (oPhi_guess > 0 and oPhi_guess < 150):
		mu_oR = 76
		mu_oPhi = 55
		std = 10
		lprior_R = -1 / 2 * (oR_guess - mu_oR)**2 / std**2
		lprior_Phi = -1 / 2 * (oPhi_guess - mu_oPhi)**2 / std**2
		prior = lprior_R + lprior_Phi
	else: #want to avoid unphysical parameters
		prior = -1e100
	return prior 

#total probablitiy
def lnprob(theta, oLOS_real, i, PA):
	lp = lprior(theta)
	# if np.isfinite(lp):
	# 	return -np.inf
	return ll(theta, oLOS_real, i, PA) + lp 

#read in fake dataset -- replace with real data eventually (won't have the oR or oP) =============================================
xi, eta, r, oLOS_real, oR_real, oPhi_real = np.loadtxt('../Data/fakedata_M31_Illustris.txt', unpack=True)
oLOS_real = oLOS_real - 40
#xi, eta, r, oLOS_real = np.loadtxt('../Data/M31_MS.txt', usecols=(0,1,2,5), unpack=True)

#prepare the fake data
#go from xi and eta to PA,deproj
def PA_calc(xi, eta): #kpc, kpc
	sine = np.sin(galaxy_PA * np.pi / 180) #sine of M33's PA
	cosine = np.cos(galaxy_PA * np.pi / 180) #cosine of M33's PA
	x = xi * cosine - eta * sine #deg
	y = eta * cosine + xi * sine #deg
	#angle stuff below
	projected = np.tan(x / y) #rad; x over y because of shifted quadrants
	deprojected = projected / np.cos(np.deg2rad(galaxy_inc)) #inclincation of M31
	return np.arctan(deprojected) * 180 / np.pi #deg; this will not return the correct quadrant exactly BUT we only care about the cosine and sine SQUARED values
PA = PA_calc(xi, eta)

#grab i from the tilted ring model
def assign_TR_params(dist, galaxy): #stellar distance to center in kpc
	#read in the titled ring parameters
	if galaxy == 'M31':
		r_kpc, i_tr = np.loadtxt('../Data/HI_PA_i_vrot.txt', usecols=(0, 2,), unpack=True) 
		print('M31')
	elif galaxy == 'M33':
		r_kpc, i_tr = np.loadtxt('../Data/M33_HI_tilted_ring.txt', usecols=(1, 4,), unpack=True)
		print('M33')

	#find the ring that each star fits in
	assigned_i = np.zeros(len(dist))
	for i in range(len(dist)):
		difference = r_kpc - dist[i] #difference between the star radii and all the rings
		pos_diff = np.array([a for a in difference if a > 0]) #eliminate rings inner to the star
		smallest_diff = min(pos_diff) #find the smallest positive difference
		ring = list(difference).index(smallest_diff) #find where in the original list the right ring is
		assigned_i[i] = i_tr[ring]
	return assigned_i
inc = assign_TR_params(r, galaxy_name)

#set the emcee run =============================================================================================================
ndim = 2 #number of parameters 
nwalkers = 300 #number of walkers exploring the parameter space simultaneously
nsteps = 300 #how long the walkers explore

#initial parameter guesses; this and priors will need to fix/play with 
oR_0 = 76
oPhi_0 = 55
oR_int = np.random.normal(loc=oR_0, scale=0.3 * oR_0, size=nwalkers)
oR_int[oR_int < 0] = 70
oPhi_int = np.random.normal(loc=oPhi_0, scale=0.1 * oPhi_0, size=nwalkers)
oPhi_int[oPhi_int < 0] = 50
pos = np.vstack((oR_int, oPhi_int)).T

#run MCMC
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(oLOS_real, inc, PA))
sampler.run_mcmc(pos, nsteps, skip_initial_state_check=True)#, progress=True)

#checking results
plt.clf()
#check mixing
for j in range(nwalkers):
	plt.plot(range(nsteps), sampler.chain[j, :,0], c='k', alpha=.2)
plt.plot([0, nsteps], [np.median(oR_real), np.median(oR_real)], c='b')
plt.xlabel('step number')
plt.ylabel('oR')
#plt.xlim(0, 10)
plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/combodata/oR_mix.png')
plt.close()
for j in range(nwalkers):
	plt.plot(range(nsteps), sampler.chain[j, :,1], c='k', alpha=.2)
plt.plot([0, nsteps], [np.median(oPhi_real), np.median(oPhi_real)], c='b')
plt.xlabel('step number')
plt.ylabel('oPhi')
#plt.xlim(0, 10)
plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/combodata/oPhi_mix.png')
plt.close()
#check for burn in
for j in range(nwalkers):
	plt.plot(range(nsteps), sampler.lnprobability[j,:], c='k', alpha=.2)
plt.title('Mean acceptance fraction: {0:.3f}'.format(np.mean(sampler.acceptance_fraction)))
plt.xlabel('step number')
plt.ylabel('Log(Probability)')
plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/combodata/burnin.png')
plt.close()
#making corner plots
samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=['oR', 'oPhi'], show_titles=True, truths=[np.median(oR_real), np.median(oPhi_real)])
fig.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/combodata/corner.png')
plt.close()





