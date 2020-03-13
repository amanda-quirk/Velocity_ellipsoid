import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve

'''
updating orginial_dispersion_MCMC.py to remove dependence on the jean's equation

This does scipy to solve for sigma_R given axis ratios; will vary the axis ratios eventually 
'''

#data models ======================================================================================================================
def sine(angle): #deg
	return np.sin(np.deg2rad(angle))

def cosine(angle): #deg
	return np.cos(np.deg2rad(angle))

#sigma R is the individual star parameter while AR1 and AR2 are global age bin parameters
def func(sigma_R, *data): #km/s, deg, deg 
	sigma_LOS, inc, PA, AR1, AR2 = data
	return -sigma_LOS**2 + sigma_R**2 * sine(inc)**2 * (sine(PA)**2 + AR1**2 * cosine(PA)**2) + AR2**2 * sigma_R**2 * cosine(inc)**2
	
#read in fake dataset -- replace with real data eventually ========================================================================
#sigmaR (km/s), sigmaPhi (km/s), sigmaZ (km/s), sigmaLOS (km/s), inclination (deg), PA (deg)
oR_real, oPhi_real, oZ_real, oLOS_real, inc, PA = np.loadtxt('/Volumes/Titan/Velocity_ellipsoid/Data/fake_data_nojeans.txt', unpack=True)

#FOR NOW I am assigning values to AR1 and AR2 -- this will be done selectively when I add in the hierarchical part ================
AR1 = np.zeros(len(oR_real)) + .8
AR2 = np.zeros(len(oR_real)) + np.sqrt(.2)

#solving for sigma_R ==============================================================================================================
oR_guess = np.zeros(len(oR_real)) + 55 #keeping it the same for each point -- this is definitely a drawback from MCMC
data = ((oLOS_real, inc, PA, AR1, AR2))
oR_solved = fsolve(func, oR_guess, args=data)


