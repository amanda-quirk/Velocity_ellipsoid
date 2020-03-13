import emcee 
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve

'''
updating orginial_dispersion_MCMC.py to remove dependence on the jean's equation

This will be the MCMC procedure to (hopefully) decompose the velocity ellipsoid
'''

AR1_value = .8
AR2_value = np.sqrt(.2)

#data models ===================================================================================================================
def sine(angle): #deg
	return np.sin(np.deg2rad(angle))

def cosine(angle): #deg
	return np.cos(np.deg2rad(angle))

#sigma R is the individual star parameter while AR1 and AR2 are global age bin parameters
def sigma_LOS(oR, oZ, oP, i, PA): #km/s, deg
	oLOS2 = (oR**2 * sine(PA)**2 + oP**2 * cosine(PA)**2) * sine(i)**2 + oZ**2 * cosine(i)**2
	return np.sqrt(oLOS2) #km/s

#emcee functions ===============================================================================================================
#logliklihood function -- right now has no errors incorporated
def ll(theta, oLOS_real, i, PA):
	oR_guess, oZ_guess, oPhi_guess = theta
	oLOS_calc = sigma_LOS(oR_guess, oZ_guess, oPhi_guess, i, PA)
	oLOS_comp = -1 / 2 * (oLOS_real - oLOS_calc)**2
	return oLOS_comp 

#priors -- assuming normal distribution but play around with this and check it later
def lprior(theta):
	#pulling priors from theta
	oR_guess, oZ_guess, oPhi_guess = theta 
	#mean and standard deviations of components -- will eventually want to change these for each age group
	if (oR_guess > 0 and oR_guess < 500) and (oZ_guess > 0 and oZ_guess < 500) and (oPhi_guess > 0 and oPhi_guess < 500):
		mu_oR = 80
		mu_oZ = AR2_value * mu_oR
		mu_oPhi = AR1_value * mu_oR
		std = 30
		#want to do for each parameter (R, z, phi components)
		lprior_R = -1 / 2 * (oR_guess - mu_oR)**2 / std**2
		lprior_Z = -1 / 2 * (oZ_guess - mu_oZ)**2 / std**2
		lprior_Phi = -1 / 2 * (oPhi_guess - mu_oPhi)**2 / std**2
		prior = lprior_R + lprior_Z + lprior_Phi
	else: #want to avoid unphysical parameters
		prior = -1e100
	return prior 

#total probablitiy
def lnprob(theta, oLOS_real, i, PA):
	lp = lprior(theta)
	# if np.isfinite(lp):
	# 	return -np.inf
	return ll(theta, oLOS_real, i, PA) + lp 

#read in fake dataset -- replace with real data eventually =====================================================================
#sigmaR (km/s), sigmaPhi (km/s), sigmaZ (km/s), sigmaLOS (km/s), inclination (deg), PA (deg)
oR_real, oPhi_real, oZ_real, oLOS_real, inc, PA = np.loadtxt('/Volumes/Titan/Velocity_ellipsoid/Data/fake_data_nojeans.txt', unpack=True)

#set the emcee run =============================================================================================================
ndim = 3 #number of parameters 
nwalkers = 500 #number of walkers exploring the parameter space simultaneously
nsteps = 170 #how long the walkers explore

#initial parameter guess -- will want to change for each age bin =================================================================
#use scipy root finder to get a good initial guess
def func(sigma_R, *data): #km/s, deg, deg 
	sigma_LOS, inc, PA, AR1, AR2 = data
	return -sigma_LOS**2 + sigma_R**2 * sine(inc)**2 * (sine(PA)**2 + AR1**2 * cosine(PA)**2) + AR2**2 * sigma_R**2 * cosine(inc)**2

#go star by star and do MCMC, saving some diagnostic output along the way=========================================================
import corner
for i in range(4):#len(oLOS_real)):
	print(i)

	#set initial positions for walkers based on data from this point
	oR_guess = np.median(oLOS_real) / np.sqrt(3) #initialize one guess -- kind of random
	data_for_root = ((oLOS_real[i], inc[i], PA[i], AR1_value, AR2_value))
	oR_scipy = fsolve(func, oR_guess, args=data_for_root)
	print(oR_scipy, oR_real[i])
	oR_guesses = abs(np.random.normal(loc=oR_scipy, scale= 1 * oR_guess, size=nwalkers)) #will want to vary this scale factor later!
	oZ_guesses = oR_guesses * AR2_value  
	oPhi_guesses = oR_guesses * AR1_value 
	pos = np.vstack((oR_guesses, oZ_guesses, oPhi_guesses)).T
	#run MCMC
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(oLOS_real, inc, PA))
	sampler.run_mcmc(pos, nsteps, skip_initial_state_check=True)

	plt.clf()
	#check mixing
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.chain[j, :,0], c='k', alpha=.2)
	plt.plot([0, nsteps], [oR_real[i], oR_real[i]], c='b')
	plt.xlabel('step number')
	plt.ylabel('oR')
	#plt.xlim(0, 10)
	plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/fakedata/mixing/{}_oR_mix.png'.format(i))
	plt.close()
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.chain[j, :,1], c='k', alpha=.2)
	plt.plot([0, nsteps], [oZ_real[i], oZ_real[i]], c='b')
	plt.xlabel('step number')
	plt.ylabel('oZ')
	#plt.xlim(0, 10)
	plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/fakedata/mixing/{}_oZ_mix.png'.format(i))
	plt.close()
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.chain[j, :,2], c='k', alpha=.2)
	plt.plot([0, nsteps], [oPhi_real[i], oPhi_real[i]], c='b')
	plt.xlabel('step number')
	plt.ylabel('oPhi')
	#plt.xlim(0, 10)
	plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/fakedata/mixing/{}_oPhi_mix.png'.format(i))
	plt.close()
	#check for burn in
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.lnprobability[j,:], c='k', alpha=.2)
	plt.title('Mean acceptance fraction: {0:.3f}'.format(np.mean(sampler.acceptance_fraction)))
	plt.xlabel('step number')
	plt.ylabel('Log(Probability)')
	plt.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/fakedata/burnin/{}_burnin.png'.format(i))
	plt.close()
	#making corner plots
	samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
	fig = corner.corner(samples, labels=['oR', 'oZ', 'oPhi'], truths=[oR_real[i], oZ_real[i], oPhi_real[i]])
	fig.savefig('/Volumes/Titan/Velocity_ellipsoid/Plots/fakedata/corner/{}_corner.png'.format(i))
	plt.close()





