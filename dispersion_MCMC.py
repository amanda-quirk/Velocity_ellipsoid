import emcee 
import numpy as np 
import matplotlib.pyplot as plt 

'''
purpose: do MCMC (not hierarchical) to first familarize self with emcee and to later compare to hierarchical MCMC
'''

#data models =============================================================================================================================
def sine(angle):
	return np.sin(np.deg2rad(angle))

def cosine(angle):
	return np.cos(np.deg2rad(angle))

def sigma_LOS(oR, oZ, oP, i, PA): #i and PA come from tilted ring table (any dispersion values?)
	oLOS2 = (oR**2 * sine(PA)**2 + oP**2 * cosine(PA)**2) * sine(i)**2 + oZ**2 * cosine(i)**2
	return np.sqrt(oLOS2)

# def va2_vavc(oR, oP, oZ, R, vc): #vc and R from titled ring table
# 	Rd = 5.76 #kpc
# 	terms = np.zeros(len(R))
# 	for j in range(len(R)):
# 		if R[j] < 10:
# 			k = 2
# 		else:
# 			k = 1
# 		terms[j] = -oR**2 * ((oP**2 / oR**2) - 1.5 + (k * R[j] / Rd) + (oZ**2 / (2 * oR**2)))
# 	return terms

# def va(oR, oP, oZ, R, vc):
# 	C = va2_vavc(oR, oP, oZ, R, vc)
# 	pos = vc + np.sqrt(vc**2 + C)
# 	neg = vc - np.sqrt(vc**2 + C)
# 	return pos, neg #want to use the negative root

#emcee functions =======================================================================================================================
#define log likelihood function -- no errors right now; can I add in errors later?
def ll(theta, oLOS_real, i, PA):#ll(theta, oLOS_real, va_real, i, PA, R, vc):
	oR_guess, oZ_guess, oPhi_guess = theta
	oLOS_calc = sigma_LOS(oR_guess, oZ_guess, oPhi_guess, i, PA)
	oLOS_comp = -1 / 2 * (oLOS_real - oLOS_calc)**2
	#va_calc = va(oR_guess, oZ_guess, oPhi_guess, R, vc)
	#va_comp = -1 / 2 * (va_real - va_calc)**2
	return oLOS_comp #+ va_comp

#assuming normal distribution -- check this though
def lprior(theta):
	#pulling priors from theta
	oR_guess, oZ_guess, oPhi_guess = theta
	#negative ratio isn't physical 
	if oR_guess < 0 or oZ_guess < 0 or oPhi_guess < 0:
		prior = -np.inf 
	#mean and standard deviations of components -- will eventually want to change these for each age group
	else:
		mu_oR = 70
		mu_oZ = np.sqrt(2) * mu_oR
		mu_oPhi = 0.8 * mu_oR
		std = 50
		#want to do for each parameter (R, z, phi components)
		lprior_R = -1 / 2 * (oR_guess - mu_oR)**2 / std**2
		lprior_Z = -1 / 2 * (oZ_guess - mu_oZ)**2 / std**2
		lprior_Phi = -1 / 2 * (oPhi_guess - mu_oPhi)**2 / std**2
		prior = lprior_R + lprior_Z + lprior_Phi
	return prior 

#total probability
def lnprob(theta, oLOS_real, i, PA):# lnprob(theta, oLOS_real, va_real, i, PA, R, vc):
	return ll(theta, oLOS_real, i, PA) + lprior(theta)#ll(theta, oLOS_real, va_real, i, PA, R, vc) + lprior(theta)

#read in fake dataset ==================================================================================================================
#header='r (kpc), inclination (deg), PA (deg), circular vel (km/s), sigmaR (km/s), sigmaPhi (km/s), sigmaZ (km/s), sigmaLOS (km/s), AD (km/s)'
r, incs, pos_angs, vcs, sigmaR, sigmaPhi, sigmaZ, sigmaLOS, AD = np.loadtxt('../Data/fake_data.txt', unpack=True)

#set the emcee run =====================================================================================================================
#fit for the components of the LOS dispersion -- again, for this exercise, we are not assuming that these components follow axis ratios
ndim = 3
nwalkers = 100
initial_guess = np.array([80, 50, 60]) #initial guess for theta
#do I want to do this in a different way -- formally make a normal distribution?
pos = [initial_guess + 40 * np.random.randn(ndim) for i in range(nwalkers)] 
pos = [abs(a) for a in pos] #start all walkers from a positive value
nsteps = 500
#checking the distribution of initial guesses
plt.hist([row[0] for row in pos])
plt.xlabel('oR guess')
plt.savefig('/Users/amandaquirk/Desktop/oR_guesses.png')
plt.close()
plt.hist([row[1] for row in pos])
plt.xlabel('oZ guess')
plt.savefig('/Users/amandaquirk/Desktop/oZ_guesses.png')
plt.close()
plt.hist([row[2] for row in pos])
plt.xlabel('oPhi guess')
plt.savefig('/Users/amandaquirk/Desktop/oPhi_guesses.png')
plt.close()

#go star by star
import corner
for i in range(len(sigmaLOS)):
	plt.clf()
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(sigmaLOS[i], incs[i], pos_angs[i]))#args=(sigmaLOS, AD, incs, pos_angs, r, vcs))
	sampler.run_mcmc(pos, nsteps)
	#check mixing
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.chain[j, :,0], c='k', alpha=.2)
	plt.plot([0, nsteps], [sigmaR[i], sigmaR[i]], c='b')
	plt.xlabel('step number')
	plt.ylabel('oR')
	plt.savefig('/Users/amandaquirk/Desktop/fake_data/mixing/{}_oR_mix.png'.format(i))
	plt.close()
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.chain[j, :,1], c='k', alpha=.2)
	plt.plot([0, nsteps], [sigmaZ[i], sigmaZ[i]], c='b')
	plt.xlabel('step number')
	plt.ylabel('oZ')
	plt.savefig('/Users/amandaquirk/Desktop/fake_data/mixing/{}_oZ_mix.png'.format(i))
	plt.close()
	for j in range(nwalkers):
		plt.plot(range(nsteps), sampler.chain[j, :,2], c='k', alpha=.2)
	plt.plot([0, nsteps], [sigmaPhi[i], sigmaPhi[i]], c='b')
	plt.xlabel('step number')
	plt.ylabel('oPhi')
	plt.savefig('/Users/amandaquirk/Desktop/fake_data/mixing/{}_oPhi_mix.png'.format(i))
	plt.close()

	#making corner plots
	samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
	fig = corner.corner(samples, labels=['oR', 'oZ', 'oPhi'], truths=[sigmaR[i], sigmaZ[i], sigmaPhi[i]])
	fig.savefig('/Users/amandaquirk/Desktop/fake_data/corner/{}_corner.png'.format(i))
	plt.close()











