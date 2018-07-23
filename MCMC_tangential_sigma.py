import numpy as np 
import matplotlib.pyplot as plt

#convert degrees to radians
# def cosine(degree):
# 	return np.cos(degree*np.pi/180.)

# def sine(degree):
# 	return np.sin(degree*np.pi/180.)

# #define chi squared
# def chisq(sig_obs, model, error):
#     #chisq=[]
#     #for i in range(len(sig_obs)):
#         #chisq.append(((sig_obs[i]-model[i])/error[i])**2.)
#     chisq=(((sig_obs-model)/error)**2.)
#     #return sum(chisq)
#     return chisq

# #log likelihood equation
# def log_likelihood(theta, PA, sig_obs, errors):
# 	A, B, C = theta
# 	model=predicted_sigma(theta, PA)
# 	#errors_log=[]
# 	#for i in range(len(sig_obs)):
# 	#	errors_log.append(np.log(errors[i]))
# 	errors_log = np.log(errors)
# 	#ll = (-0.5*len(sig_obs)*np.log(2.*np.pi))- np.sum(errors_log)-0.5*chisq(sig_obs, model,errors)
# 	ll = (-0.5*np.log(2.*np.pi)) - np.sum(errors_log)-0.5*chisq(sig_obs, model,errors)
# 	return ll

# def claire_ll(theta, PA, sig_obs, errors):
# 	A, B, C = theta
# 	C=1
# 	model=predicted_sigma(theta, PA)
# 	stat = (sig_obs - model)**2 / (2 * errors**2) + 0.5 * np.log10(2 * np.pi * A * A)
# 	return -stat 

# #log priors -- right now these are uniform priors, might be worth doing log uniform priors later (Dark Matter Final Project)
# def log_prior(theta, bounds):
# 	A, B, C = theta
# 	Amin, Amax, Bmin, Bmax, Cmin, Cmax=bounds
# 	if Amin < A < Amax and Bmin < B < Bmax and Cmin < C < Cmax:
# 		return 0.0
# 	return -np.inf #guess doesn't fit into prior range

# #probability of parameters being correct
# def log_prob(theta, PA, sig_obs, errors):#data, model, errors, theta, bounds):
# 	lp= log_prior(theta, bounds)
# 	if not np.isfinite(lp):
# 		return -np.inf 
# 	return lp + claire_ll(theta, PA, sig_obs, errors)

# #model
# def predicted_sigma(theta, PA):
# 	A, B, C = theta
# 	i=77. #check with tilted ring model
# 	sq=(((sine(PA)**2.)+(B*cosine(PA)**2.))*(sine(i)**2.)) + (C*cosine(i)**2.)
# 	return A*np.sqrt(sq)

# #calculate the PA
# def x(xi, eta):
# 	xi_deg=xi/13.67
# 	eta_deg=eta/13.67
# 	sine=np.sin(float(37*np.pi) / 180)
# 	cosine=np.cos(float(37*np.pi) / 180)
# 	x=(xi_deg*cosine)-(eta_deg*sine)
# 	return x 

# def y(xi, eta):
# 	xi_deg=xi/13.67
# 	eta_deg=eta/13.67
# 	sine=np.sin(float(37*np.pi) / 180)
# 	cosine=np.cos(float(37*np.pi) / 180)
# 	y=(eta_deg*cosine)+(xi_deg*sine)
# 	return y

# def PA(x,y): #defined as west of north as positioned in maps (inverted x axis)
# 	if x>0:
# 		rad=np.arctan(float(y)/x)
# 		deg=90-(float(rad*180)/np.pi)
# 	else:
# 		rad=np.arctan(float(y)/x)
# 		deg=270-(float(rad*180)/np.pi)
# 	return deg + 37

# #read in data
# #is this the right verr??
# MS_xi, MS_eta,  MS_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)
# AGy_xi, AGy_eta, AGy_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)
# AGo_xi, AGo_eta, AGo_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)
# RG_xi, RG_eta, RG_var=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_smoothed_chemin.txt', usecols=(0,1,4,), unpack=True)

# MS_verr=np.loadtxt('/Users/amandaquirk/Desktop/MS_dispersion_err.txt')
# AGy_verr=np.loadtxt('/Users/amandaquirk/Desktop/AGy_dispersion_err.txt')
# AGo_verr=np.loadtxt('/Users/amandaquirk/Desktop/AGo_dispersion_err.txt')
# RG_verr=np.loadtxt('/Users/amandaquirk/Desktop/RG_dispersion_err.txt')

# MS_x=x(MS_xi, MS_eta) 
# MS_y=y(MS_xi, MS_eta)
# AGy_x=x(AGy_xi, AGy_eta) 
# AGy_y=y(AGy_xi, AGy_eta)
# AGo_x=x(AGo_xi, AGo_eta) 
# AGo_y=y(AGo_xi, AGo_eta)
# RG_x=x(RG_xi, RG_eta) 
# RG_y=y(RG_xi, RG_eta)

# MS_PA=[]
# for i in range(len(MS_x)):
# 	MS_PA.append(PA(MS_x[i], MS_y[i]))
# AGy_PA=[]
# for i in range(len(AGy_x)):
# 	AGy_PA.append(PA(AGy_x[i], AGy_y[i]))
# AGo_PA=[]
# for i in range(len(AGo_x)):
# 	AGo_PA.append(PA(AGo_x[i], AGo_y[i]))
# RG_PA=[]
# for i in range(len(RG_x)):
# 	RG_PA.append(PA(RG_x[i], RG_y[i]))

# #terms of MCMC
# jmax=1000
# iteration=100

# #MS
# #initial guesses
# sig_R_guess, o_r_guess, z_r_guess=np.array([29, 1.,0.2])
# bounds=[0,150,.5,1.5,0.,1.]

# MS_A=[]
# MS_B=[]
# MS_C=[]
# for i in range(len(MS_xi)):
# #running MCMC
# 	As=[]
# 	Bs=[]
# 	Cs=[]
# 	log_likelihoods=[]
# 	accepted_steps=[]
	
# 	theta_0=np.array([sig_R_guess, o_r_guess, z_r_guess])
# 	PA=MS_PA[i]
# 	sig_obs=MS_var[i]
# 	errors=MS_verr[i]
	
# 	for j in range(jmax):
# 		log_likelihood_0=log_prob(theta_0, PA, sig_obs, errors)
	
# 		A_trial=np.random.normal(theta_0[0], .2*sig_R_guess)
# 		B_trial=np.random.normal(theta_0[1], .2*o_r_guess)
# 		C_trial=np.random.normal(theta_0[2], .2*z_r_guess)
# 		theta_trial=np.array([A_trial,B_trial,C_trial])
	
# 		log_likelihood_trial=log_prob(theta_trial, PA, sig_obs, errors)
	
# 		if log_likelihood_trial - log_likelihood_0 > 0:
# 			As.append(A_trial)
# 			Bs.append(B_trial)
# 			Cs.append(C_trial)
# 			log_likelihoods.append(log_likelihood_trial)
# 			accepted_steps.append(j)
# 			theta_0=np.array([A_trial, B_trial, C_trial])
# 		else: #Metropolis Rule
# 			if np.exp(log_likelihood_trial- log_likelihood_0)> np.random.uniform(0,1):
# 				As.append(A_trial)
# 				Bs.append(B_trial)
# 				Cs.append(C_trial)
# 				log_likelihoods.append(log_likelihood_trial)
# 				accepted_steps.append(j)
# 				theta_0=np.array([A_trial, B_trial, C_trial])

# 	#*

# 	A_burn=[]
# 	B_burn=[]
# 	C_burn=[]
# 	log_likelihood_burn=[]
# 	for i in range(len(accepted_steps)):
# 		if accepted_steps[i]>iteration:
# 			A_burn.append(As[i])
# 			B_burn.append(Bs[i])
# 			C_burn.append(Cs[i])
# 			log_likelihood_burn.append(log_likelihoods[i])
	
# 	n=np.argmax(log_likelihood_burn)
# 	MS_A.append(A_burn[n])
# 	MS_B.append(B_burn[n])
# 	MS_C.append(C_burn[n])

# print(np.mean(MS_A), np.mean(MS_B), np.mean(MS_C), len(MS_A))

# #young AGB
# #initial guesses
# sig_R_guess, o_r_guess, z_r_guess=np.array([79, 0.49,0.2])
# bounds=[0,150,0.,1.,0.,1.]

# AGy_A=[]
# AGy_B=[]
# AGy_C=[]
# for i in range(len(AGy_xi)):
# #running MCMC
# 	As=[]
# 	Bs=[]
# 	Cs=[]
# 	log_likelihoods=[]
# 	accepted_steps=[]
	
# 	theta_0=np.array([sig_R_guess, o_r_guess, z_r_guess])
# 	PA=AGy_PA[i]
# 	sig_obs=AGy_var[i]
# 	errors=AGy_verr[i]
	
# 	for j in range(jmax):
# 		log_likelihood_0=log_prob(theta_0, PA, sig_obs, errors)
	
# 		A_trial=np.random.normal(theta_0[0], .2*sig_R_guess)
# 		B_trial=np.random.normal(theta_0[1], .2*o_r_guess)
# 		C_trial=np.random.normal(theta_0[2], .2*z_r_guess)
# 		theta_trial=np.array([A_trial,B_trial,C_trial])
	
# 		log_likelihood_trial=log_prob(theta_trial, PA, sig_obs, errors)
	
# 		if log_likelihood_trial - log_likelihood_0 > 0:
# 			As.append(A_trial)
# 			Bs.append(B_trial)
# 			Cs.append(C_trial)
# 			log_likelihoods.append(log_likelihood_trial)
# 			accepted_steps.append(j)
# 			theta_0=np.array([A_trial, B_trial, C_trial])
# 		else: #Metropolis Rule
# 			if np.exp(log_likelihood_trial- log_likelihood_0)> np.random.uniform(0,1):
# 				As.append(A_trial)
# 				Bs.append(B_trial)
# 				Cs.append(C_trial)
# 				log_likelihoods.append(log_likelihood_trial)
# 				accepted_steps.append(j)
# 				theta_0=np.array([A_trial, B_trial, C_trial])

# 	#*

# 	A_burn=[]
# 	B_burn=[]
# 	C_burn=[]
# 	log_likelihood_burn=[]
# 	for i in range(len(accepted_steps)):
# 		if accepted_steps[i]>iteration:
# 			A_burn.append(As[i])
# 			B_burn.append(Bs[i])
# 			C_burn.append(Cs[i])
# 			log_likelihood_burn.append(log_likelihoods[i])
	
# 	n=np.argmax(log_likelihood_burn)
# 	AGy_A.append(A_burn[n])
# 	AGy_B.append(B_burn[n])
# 	AGy_C.append(C_burn[n])

# print(np.mean(AGy_A), np.mean(AGy_B), np.mean(AGy_C), len(AGy_A))

# #old AGB
# #initial guesses
# sig_R_guess, o_r_guess, z_r_guess=np.array([79, 0.49,0.2])
# bounds=[0,150,0.,1.,0.,1.]

# AGo_A=[]
# AGo_B=[]
# AGo_C=[]
# for i in range(len(AGo_xi)):
# #running MCMC
# 	As=[]
# 	Bs=[]
# 	Cs=[]
# 	log_likelihoods=[]
# 	accepted_steps=[]
	
# 	theta_0=np.array([sig_R_guess, o_r_guess, z_r_guess])
# 	PA=AGo_PA[i]
# 	sig_obs=AGo_var[i]
# 	errors=AGo_verr[i]
	
# 	for j in range(jmax):
# 		log_likelihood_0=log_prob(theta_0, PA, sig_obs, errors)
	
# 		A_trial=np.random.normal(theta_0[0], .2*sig_R_guess)
# 		B_trial=np.random.normal(theta_0[1], .2*o_r_guess)
# 		C_trial=np.random.normal(theta_0[2], .2*z_r_guess)
# 		theta_trial=np.array([A_trial,B_trial,C_trial])
	
# 		log_likelihood_trial=log_prob(theta_trial, PA, sig_obs, errors)
	
# 		if log_likelihood_trial - log_likelihood_0 > 0:
# 			As.append(A_trial)
# 			Bs.append(B_trial)
# 			Cs.append(C_trial)
# 			log_likelihoods.append(log_likelihood_trial)
# 			accepted_steps.append(j)
# 			theta_0=np.array([A_trial, B_trial, C_trial])
# 		else: #Metropolis Rule
# 			if np.exp(log_likelihood_trial- log_likelihood_0)> np.random.uniform(0,1):
# 				As.append(A_trial)
# 				Bs.append(B_trial)
# 				Cs.append(C_trial)
# 				log_likelihoods.append(log_likelihood_trial)
# 				accepted_steps.append(j)
# 				theta_0=np.array([A_trial, B_trial, C_trial])

# 	#*

# 	A_burn=[]
# 	B_burn=[]
# 	C_burn=[]
# 	log_likelihood_burn=[]
# 	for i in range(len(accepted_steps)):
# 		if accepted_steps[i]>iteration:
# 			A_burn.append(As[i])
# 			B_burn.append(Bs[i])
# 			C_burn.append(Cs[i])
# 			log_likelihood_burn.append(log_likelihoods[i])
	
# 	n=np.argmax(log_likelihood_burn)
# 	AGo_A.append(A_burn[n])
# 	AGo_B.append(B_burn[n])
# 	AGo_C.append(C_burn[n])

# print(np.mean(AGo_A), np.mean(AGo_B), np.mean(AGo_C), len(AGo_A))

# #RGB
# #initial guesses
# sig_R_guess, o_r_guess, z_r_guess=np.array([110, 0.64,0.2])
# bounds=[0,150,0.,1.,0.,1.]

# RG_A=[]
# RG_B=[]
# RG_C=[]
# for i in range(len(RG_xi)):
# #running MCMC
# 	As=[]
# 	Bs=[]
# 	Cs=[]
# 	log_likelihoods=[]
# 	accepted_steps=[]
	
# 	theta_0=np.array([sig_R_guess, o_r_guess, z_r_guess])
# 	PA=RG_PA[i]
# 	sig_obs=RG_var[i]
# 	errors=RG_verr[i]
	
# 	for j in range(jmax):
# 		log_likelihood_0=log_prob(theta_0, PA, sig_obs, errors)
	
# 		A_trial=np.random.normal(theta_0[0], .2*sig_R_guess)
# 		B_trial=np.random.normal(theta_0[1], .2*o_r_guess)
# 		C_trial=np.random.normal(theta_0[2], .2*z_r_guess)
# 		theta_trial=np.array([A_trial,B_trial,C_trial])
	
# 		log_likelihood_trial=log_prob(theta_trial, PA, sig_obs, errors)
	
# 		if log_likelihood_trial - log_likelihood_0 > 0:
# 			As.append(A_trial)
# 			Bs.append(B_trial)
# 			Cs.append(C_trial)
# 			log_likelihoods.append(log_likelihood_trial)
# 			accepted_steps.append(j)
# 			theta_0=np.array([A_trial, B_trial, C_trial])
# 		else: #Metropolis Rule
# 			if np.exp(log_likelihood_trial- log_likelihood_0)> np.random.uniform(0,1):
# 				As.append(A_trial)
# 				Bs.append(B_trial)
# 				Cs.append(C_trial)
# 				log_likelihoods.append(log_likelihood_trial)
# 				accepted_steps.append(j)
# 				theta_0=np.array([A_trial, B_trial, C_trial])

# 	#*

# 	A_burn=[]
# 	B_burn=[]
# 	C_burn=[]
# 	log_likelihood_burn=[]
# 	for i in range(len(accepted_steps)):
# 		if accepted_steps[i]>iteration:
# 			A_burn.append(As[i])
# 			B_burn.append(Bs[i])
# 			C_burn.append(Cs[i])
# 			log_likelihood_burn.append(log_likelihoods[i])
	
# 	n=np.argmax(log_likelihood_burn)
# 	RG_A.append(A_burn[n])
# 	RG_B.append(B_burn[n])
# 	RG_C.append(C_burn[n])

# print(np.mean(RG_A), np.mean(RG_B), np.mean(RG_C), len(RG_A))

# # import corner
# # data=np.vstack([A_burn, B_burn, C_burn])
# # figure=corner.corner(data.T, labels=['A', 'B', 'C'], quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
# # plt.savefig('/Users/amandaquirk/Desktop/corner.png')
# # plt.close()

# def tangential_component(A, B):
# 	return A * np.sqrt(B)

# MS_tang=tangential_component(MS_A, MS_B)
# AGy_tang=tangential_component(AGy_A, AGy_B)
# AGo_tang=tangential_component(AGo_A, AGo_B)
# RG_tang=tangential_component(RG_A, RG_B)

# file=open('/Users/amandaquirk/Desktop/MS_components.txt', 'w')
# file.write('#sig_tan, sig_R, (sig_tan/sig_R)^2, (sig_z/sig_R)^2 \n')
# for i in range(len(MS_xi)):
# 	file.write('{} {} {} {}\n'.format(MS_tang[i], MS_A[i], MS_B[i], MS_C[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGy_components.txt', 'w')
# file.write('#sig_tan, sig_R, (sig_tan/sig_R)^2, (sig_z/sig_R)^2 \n')
# for i in range(len(AGy_xi)):
# 	file.write('{} {} {} {}\n'.format(AGy_tang[i], AGy_A[i], AGy_B[i], AGy_C[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/AGo_components.txt', 'w')
# file.write('#sig_tan, sig_R, (sig_tan/sig_R)^2, (sig_z/sig_R)^2 \n')
# for i in range(len(AGo_xi)):
# 	file.write('{} {} {} {}\n'.format(AGo_tang[i], AGo_A[i], AGo_B[i], AGo_C[i]))
# file.close()

# file=open('/Users/amandaquirk/Desktop/RG_components.txt', 'w')
# file.write('#sig_tan, sig_R, (sig_tan/sig_R)^2, (sig_z/sig_R)^2 \n')
# for i in range(len(RG_xi)):
# 	file.write('{} {} {} {}\n'.format(RG_tang[i], RG_A[i], RG_B[i], RG_C[i]))
# file.close()

MS_tang=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_components.txt', usecols=(0,), unpack=True)
AGy_tang=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_components.txt', usecols=(0,), unpack=True)
AGo_tang=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_components.txt', usecols=(0,), unpack=True)
RG_tang=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_components.txt', usecols=(0,), unpack=True)

MS_ad=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_ad.txt')
AGy_ad=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_ad.txt')
AGo_ad=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_ad.txt')
RG_ad=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_ad.txt')

from matplotlib import rc 
from matplotlib.ticker import MaxNLocator

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter([a*a for a in MS_tang], MS_ad, s=5, c='b')
axes[0].plot([0,22500], [-10,60], c='w', linewidth=2)
axes[0].plot([0,22500], [-10,60], c='dimgrey')
axes[1].scatter([a*a for a in AGy_tang], AGy_ad, s=5, c='m')
axes[1].plot([0,22500], [10,80], c='dimgrey')
axes[2].scatter([a*a for a in AGo_tang], AGo_ad, s=5, c='k')
axes[2].plot([0,22500], [10,80], c='w', linewidth=3)
axes[2].plot([0,22500], [10,80], c='dimgrey')
axes[3].scatter([a*a for a in RG_tang], RG_ad, s=5, c='r')
axes[3].plot([0,22500], [20,110], c='dimgrey')
axes[0].annotate('MS', xy=(19000,-9), horizontalalignment='right', fontsize=12)
axes[1].annotate('young AGB', xy=(19000,-9), horizontalalignment='right', fontsize=12)
axes[2].annotate('older AGB', xy=(19000,-9), horizontalalignment='right', fontsize=12)
axes[3].annotate('RGB', xy=(19000,-9), horizontalalignment='right', fontsize=12)

for ax in axes:
	ax.set_xlim(0, 20000)
	#ax.set_ylabel(r'$\rm Asymmetric\ Drift\ (km\ s^{-1})$', fontsize=13)
	ax.set_ylim(-20,140)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=2)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
axes[3].set_xlabel(r'$\sigma_{\phi}^{2}$', fontsize=13)
nbins = len(axes[0].get_yticklabels())
axes[0].yaxis.set_major_locator(MaxNLocator(nbins=nbins))#, prune='upper'))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)	        	
f.text(0.008, 0.5, r'$\rm Asymmetric\ Drift:\ v_{\rm a}\ (km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)        
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/lag_tang_sigma.pdf', bbox_inches='tight')

# plt.figure(figsize=(10, 6))
# plt.scatter([a*a for a in MS_tang], MS_ad, s=5, c='k')
# plt.plot([0,22500], [-10,60], c='b')
# plt.xlabel(r'$\sigma_{\phi}^{2}$')
# plt.ylabel('Asymmetric Drift')
# plt.ylim(-20,140)
# plt.xlim(0, 20000)
# plt.savefig('/Users/amandaquirk/Desktop/MS_ad_sig_tan.png')
# plt.close()

# plt.figure(figsize=(10, 6))
# plt.scatter([a*a for a in AGy_tang], AGy_ad, s=5, c='k')
# plt.plot([0,22500], [10,80], c='m')
# plt.xlabel(r'$\sigma_{\phi}^{2}$')
# plt.ylabel('Asymmetric Drift')
# plt.ylim(-20,140)
# plt.xlim(0, 20000)
# plt.savefig('/Users/amandaquirk/Desktop/AGy_ad_sig_tan.png')
# plt.close()

# plt.figure(figsize=(10, 6))
# plt.scatter([a*a for a in AGo_tang], AGo_ad, s=5, c='k')
# plt.plot([0,22500], [10,80], c='m')
# plt.xlabel(r'$\sigma_{\phi}^{2}$')
# plt.ylabel('Asymmetric Drift')
# plt.ylim(-20,140)
# plt.xlim(0, 20000)
# plt.savefig('/Users/amandaquirk/Desktop/AGo_ad_sig_tan.png')
# plt.close()

# plt.figure(figsize=(10, 6))
# plt.scatter([a*a for a in RG_tang], RG_ad, s=5, c='k')
# plt.plot([0,22500], [20,110], c='r')
# plt.xlabel(r'$\sigma_{\phi}^{2}$')
# plt.ylabel('Asymmetric Drift')
# plt.ylim(-20,140)
# plt.xlim(0, 20000)
# plt.savefig('/Users/amandaquirk/Desktop/RG_ad_sig_tan.png')
# plt.close()

#below goes at the star
	#below checks walking and burn in 
	# print(np.std(As))
	# print(np.std(Bs))
	# print(np.std(Cs))
	
	# plt.plot(accepted_steps, log_likelihoods)
	# plt.xlabel('Iteration', fontsize=16)
	# #plt.xlim(0,100)
	# plt.ylabel('Log Likelihood', fontsize=16)
	# plt.title('Searching for the Burn-In')
	# acceptance=(len(accepted_steps)/jmax)*100
	# #print('Accepted Rate={}%'.format(acceptance))
	# plt.savefig('/Users/amandaquirk/Desktop/looking_forburn.png')
	# plt.close()
	
	# plt.plot(accepted_steps, log_likelihoods)
	# plt.xlabel('Iteration', fontsize=16)
	# plt.xlim(0,100)
	# plt.ylabel('Log Likelihood', fontsize=16)
	# plt.title('Searching for the Burn-In')
	# acceptance=(len(accepted_steps)/jmax)*100
	# plt.savefig('/Users/amandaquirk/Desktop/zoom.png')
	# plt.close()
	
	# plt.plot(accepted_steps, As)
	# plt.savefig('/Users/amandaquirk/Desktop/A_walking.png')
	# plt.close()
	
	# plt.plot(accepted_steps, Bs)
	# plt.savefig('/Users/amandaquirk/Desktop/B_walking.png')
	# plt.close()
	
	# plt.plot(accepted_steps, Cs)
	# plt.savefig('/Users/amandaquirk/Desktop/C_walking.png')
	# plt.close()
	
	#burn in safely by iteration 2000





