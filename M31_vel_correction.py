import numpy as np 
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as const
from astropy.time import Time

#original data ================================================================================================================================
hdu = fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/subMasterSPLASH.fits', memmap = True)
data = hdu[1].data
ra = data['RA']
dec = data['Dec']
redshift = data['Z']
time = data['MJD']
aband = data['ABAND']
print('read in subMasterSPLASH')

#read in data from age groups =================================================================================================================
def read_data(agebin):
	#xi, eta, verr, n, HI, index in subMasterSPLASH
	age_data = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/{}_individual_chemin.txt'.format(agebin), usecols=(0, 1, 3, 4, 5, 8), unpack=True)
	return age_data

MS_data = read_data('MS')
AGy_data = read_data('AGy')
AGo_data = read_data('AGo')
RG_data = read_data('RG')
print('read in age data')

#correct the velocities ======================================================================================================================
keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)

def correct_vel(data_array):
	inds = data_array[-1]
	inds = [int(a) for a in inds]

	ras = ra[inds]
	decs = dec[inds]
	sc = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle, u.deg))
	times = time[inds]
	zs = redshift[inds]
	abands = aband[inds]

	heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(times, format='mjd'), location=keck) 
	heliocorr_km_s = heliocorr.to(u.km/u.s) 
	vraw = zs * const.c.to(u.km/u.s)
	vcorr = vraw + heliocorr_km_s - abands * const.c.to(u.km/u.s)

	return vcorr.value #km/s

MS_v = correct_vel(MS_data)
AGy_v = correct_vel(AGy_data)
AGo_v = correct_vel(AGo_data)
RG_v = correct_vel(RG_data)
print('did velocity corrections')

#smoothing data ============================================================================================================================
#function to calculate the weights
def calc_weights(err):
        return 1 / (err**2)

def normed_weight(w):
        sum_weights=sum(w)
        return w / sum_weights

#function does the weighted meean
def weighted_mean(data,norm_w):
	return sum(data * norm_w)

#function does the weighted RMSE
def weighted_rmse(norm_w, data, mean):
	diff_sq = (data - mean)**2
	return np.sqrt(sum(diff_sq * norm_w))

def smoothing(data_array, velocities, circleSize):
	smoothed_v = []
	dispersion = []
	#below these values are not actually smoothed, just saving the ones for good centers
	xi_goodcenter = []
	eta_goodcenter = []
	smoothed_n = []
	smoothed_err = []
	smoothed_HI = []
	smoothed_ind = []

	#remove stars that have unreliable velocities
	reliable = abs(velocities) < 1000
	ras = data_array[0][reliable]
	decs = data_array[1][reliable]
	errs = data_array[2][reliable]
	ns = data_array[3][reliable]
	HI = data_array[4][reliable]
	ind = data_array[5][reliable]
	velocities = velocities[reliable]

	weight = calc_weights(errs) #error is already adjusted so can just calculate the weights
	sc = SkyCoord(ra=ras, dec=decs, unit=(u.deg,u.deg))
	for i in range(len(ras)):
		c1 = SkyCoord(ras[i], decs[i], unit=(u.deg,u.deg)) #go through all coordinates one at a time
		sep = c1.separation(sc)
		good = sep.arcsecond < circleSize #put stars into smoothing circle of this size
		velocities_circ = velocities[good]
		weight_circ = weight[good]
		if len(velocities_circ) > 15: #only want circles with at least 15 stars
			normed_weights = normed_weight(weight_circ)
			smoothed_v.append(weighted_mean(velocities_circ, normed_weights)) #average the velocites
			dispersion.append(weighted_rmse(normed_weights, velocities_circ, weighted_mean(velocities_circ, normed_weights)))
			xi_goodcenter.append(ras[i] * 13.67) #kpc
			eta_goodcenter.append(decs[i]* 13.67) #kpc
			smoothed_n.append(ns[i])
			smoothed_err.append(errs[i])
			smoothed_HI.append(HI[i])
			smoothed_ind.append(ind[i])
	return xi_goodcenter, eta_goodcenter, smoothed_v, smoothed_err, dispersion, smoothed_HI, smoothed_n, smoothed_ind

MS_smoothed_data = smoothing(MS_data, MS_v, 200)
print('done with MS smoothing')
AGy_smoothed_data = smoothing(AGy_data, AGy_v, 275)
print('done with AGy smoothing')
AGo_smoothed_data = smoothing(AGo_data, AGo_v, 275)
print('done with AGo smoothing')
RG_smoothed_data = smoothing(RG_data, RG_v, 200)
print('done with RG smoothing')

# #calculating rotation velocity =================================================================================================================
#deproject coordinates
def x(xi, eta): #xi and eta in kpc
	xi_deg = xi / 13.67
	eta_deg = eta / 13.67
	sine = np.sin(37 * np.pi / 180)
	cosine = np.cos(37 * np.pi / 180)
	x =(xi_deg * cosine) - (eta_deg * sine)
	return x 

def y(xi, eta): #xi and eta in kpc
	xi_deg = xi / 13.67
	eta_deg = eta / 13.67
	sine = np.sin(37 * np.pi / 180)
	cosine = np.cos(37 * np.pi / 180)
	y = (eta_deg * cosine) + (xi_deg * sine)
	return y

#calculate the position angle -- THIS FUNCTION NEEDS TO BE MODIFIED IF Y < 0
def PA(xi, eta): #xi and eta in kpc
	x_coord = x(xi, eta)
	y_coord = y(xi, eta)
	deg = np.zeros_like(x_coord)
	for i in range(len(x_coord)):
		if x_coord[i] > 0:
			rad = np.arctan(y_coord[i] / x_coord[i])
			deg[i] = 90 - rad * 180 / np.pi
		else:
			rad = np.arctan(y_coord[i] / x_coord[i])
			deg[i] = 270 - rad * 180 / np.pi
	return deg + 37 #incorporate tilt of M31

#calculate deprojected radial distance
def distance(xi, eta): #xi and eta in kpc
	x_coord = x(xi, eta)
	y_coord = y(xi, eta)
	inclination_factor = np.cos(77 * np.pi / 180)**2
	ang_dist = np.sqrt(y_coord**2 + x_coord**2 / inclination_factor)
	return ang_dist * 13.67

#bring in the ringed HI data for tilted ring model
HI_r, HI_PA, HI_i, HI_v = np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/HI_PA_i_vrot.txt', unpack=True)
HI_PA = np.array(HI_PA)
HI_i = np.array(HI_i)

def find_nearest_ring(radius):
	indices = np.zeros_like(radius)
	for i in range(len(radius)):
		diff = abs(HI_r - radius[i])
		indices[i] = int(np.argmin(diff))
	return indices #index of the radius closest to the star's radius

#rotation velocity using tilted ring model
def Vrot_tilted_ring(xi, eta, v): #LOS v
	PA_star = PA(xi, eta)
	r = distance(xi, eta)

	#HI ring parameters 
	HI_inds = find_nearest_ring(r)
	PA_ring = np.zeros(len(HI_inds))
	i_ring = np.zeros(len(HI_inds))
	for i in range(len(r)):
		N = int(HI_inds[i])
		PA_ring[i] = HI_PA[N]
		i_ring[i] = HI_i[N]

	deg = np.pi / 180
	vsys = -300 #km/s, as defined in Claire's thesis
	A = (v - vsys) / np.sin(i_ring * deg)
	B = np.tan((PA_ring - PA_star) * deg)**2
	C = np.cos(i_ring * deg)**2
	rotation_velocity = A * np.sqrt(1 + B / C)
	return abs(rotation_velocity)

#calculating rotation velocities
MS_vrot = Vrot_tilted_ring(np.array(MS_smoothed_data[0]), np.array(MS_smoothed_data[1]), np.array(MS_smoothed_data[2]))
MS_HI_vrot = Vrot_tilted_ring(np.array(MS_smoothed_data[0]), np.array(MS_smoothed_data[1]), np.array(MS_smoothed_data[5]))
MS_r = distance(np.array(MS_smoothed_data[0]), np.array(MS_smoothed_data[1]))
AGy_vrot = Vrot_tilted_ring(np.array(AGy_smoothed_data[0]), np.array(AGy_smoothed_data[1]), np.array(AGy_smoothed_data[2]))
AGy_HI_vrot = Vrot_tilted_ring(np.array(AGy_smoothed_data[0]), np.array(AGy_smoothed_data[1]), np.array(AGy_smoothed_data[5]))
AGy_r = distance(np.array(AGy_smoothed_data[0]), np.array(AGy_smoothed_data[1]))
AGo_vrot = Vrot_tilted_ring(np.array(AGo_smoothed_data[0]), np.array(AGo_smoothed_data[1]), np.array(AGo_smoothed_data[2]))
AGo_HI_vrot = Vrot_tilted_ring(np.array(AGo_smoothed_data[0]), np.array(AGo_smoothed_data[1]), np.array(AGo_smoothed_data[5]))
AGo_r = distance(np.array(AGo_smoothed_data[0]), np.array(AGo_smoothed_data[1]))
RG_vrot = Vrot_tilted_ring(np.array(RG_smoothed_data[0]), np.array(RG_smoothed_data[1]), np.array(RG_smoothed_data[2]))
RG_HI_vrot = Vrot_tilted_ring(np.array(RG_smoothed_data[0]), np.array(RG_smoothed_data[1]), np.array(RG_smoothed_data[5]))
RG_r = distance(np.array(RG_smoothed_data[0]), np.array(RG_smoothed_data[1]))

#save the data!!
np.savetxt('/Users/amandaquirk/Documents/Ellipsoid/Data/M31_MS.txt', np.c_[MS_smoothed_data[0], MS_smoothed_data[1], MS_r, MS_smoothed_data[2], MS_smoothed_data[3], MS_smoothed_data[4], MS_vrot, MS_smoothed_data[5], MS_HI_vrot, MS_smoothed_data[6], MS_smoothed_data[7]], delimiter=' ', header='xi (kpc), eta (kpc), r (kpc), LOS_v (km/s), verr (km/s), dispersion (km/s), vrot (km/s), HI v (km/s), HI vrot (km/s), n components, subMasterSPLASH index')
np.savetxt('/Users/amandaquirk/Documents/Ellipsoid/Data/M31_AGy.txt', np.c_[AGy_smoothed_data[0], AGy_smoothed_data[1], AGy_r, AGy_smoothed_data[2], AGy_smoothed_data[3], AGy_smoothed_data[4], AGy_vrot, AGy_smoothed_data[5], AGy_HI_vrot, AGy_smoothed_data[6], AGy_smoothed_data[7]], delimiter=' ', header='xi (kpc), eta (kpc), r (kpc), LOS_v (km/s), verr (km/s), dispersion (km/s), vrot (km/s), HI v (km/s), HI vrot (km/s), n components, subMasterSPLASH index')
np.savetxt('/Users/amandaquirk/Documents/Ellipsoid/Data/M31_AGo.txt', np.c_[AGo_smoothed_data[0], AGo_smoothed_data[1], AGo_r, AGo_smoothed_data[2], AGo_smoothed_data[3], AGo_smoothed_data[4], AGo_vrot, AGo_smoothed_data[5], AGo_HI_vrot, AGo_smoothed_data[6], AGo_smoothed_data[7]], delimiter=' ', header='xi (kpc), eta (kpc), r (kpc), LOS_v (km/s), verr (km/s), dispersion (km/s), vrot (km/s), HI v (km/s), HI vrot (km/s), n components, subMasterSPLASH index')
np.savetxt('/Users/amandaquirk/Documents/Ellipsoid/Data/M31_RG.txt', np.c_[RG_smoothed_data[0], RG_smoothed_data[1], RG_r, RG_smoothed_data[2], RG_smoothed_data[3], RG_smoothed_data[4], RG_vrot, RG_smoothed_data[5], RG_HI_vrot, RG_smoothed_data[6], RG_smoothed_data[7]], delimiter=' ', header='xi (kpc), eta (kpc), r (kpc), LOS_v (km/s), verr (km/s), dispersion (km/s), vrot (km/s), HI v (km/s), HI vrot (km/s), n components, subMasterSPLASH index')
