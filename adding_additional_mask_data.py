import numpy as np 
from astropy.io import fits 
import matplotlib.pyplot as plt 
from astropy import units as u
from astropy.coordinates import SkyCoord

#import data
hdu = fits.open('/Users/amandaquirk/Documents/AsymmetricDrift/Data/submasterSPLASH.fits', memmap=True)
data = hdu[1].data 

mask = data['MASK']
RA = data['RA']
Dec = data['DEC']
z = data['Z']
zqual = data['ZQUAL']
temp = data['TEFF']
likelihood = data['LIKELIHOOD']
ID = data['OBJNAME']

# get data from new region
RA_region = np.array([a for a, b, c, d in zip(RA, mask, zqual, likelihood) if b=='SE7' or b=='SE8' or b=='SE9' or b=='M32' or b=='M32_2' or b=='M32_3' or b=='M32_4' or b=='M32_5' or b=='M32_6' and c>=3 and np.isnan(d)==True or d>0])
Dec_region = np.array([a for a, b, c, d in zip(Dec, mask, zqual, likelihood) if b=='SE7' or b=='SE8' or b=='SE9' or b=='M32' or b=='M32_2' or b=='M32_3' or b=='M32_4' or b=='M32_5' or b=='M32_6' and c>=3 and np.isnan(d)==True or d>0])
v_region = np.array([a * 3 * 10**5 for a, b, c, d in zip(z, mask, zqual, likelihood) if b=='SE7' or b=='SE8' or b=='SE9' or b=='M32' or b=='M32_2' or b=='M32_3' or b=='M32_4' or b=='M32_5' or b=='M32_6' and c>=3 and np.isnan(d)==True or d>0])
temp_region = np.array([a for a, b, c, d in zip(temp, mask, zqual, likelihood) if b=='SE7' or b=='SE8' or b=='SE9' or b=='M32' or b=='M32_2' or b=='M32_3' or b=='M32_4' or b=='M32_5' or b=='M32_6' and c>=3 and np.isnan(d)==True or d>0])
ID_region = np.array([a for a, b, c, d in zip(ID, mask, zqual, likelihood) if b=='SE7' or b=='SE8' or b=='SE9' or b=='M32' or b=='M32_2' or b=='M32_3' or b=='M32_4' or b=='M32_5' or b=='M32_6' and c>=3 and np.isnan(d)==True or d>0])

#get rid of nan values
ok = np.isnan(temp_region) == False

RA_region = RA_region[ok]
Dec_region = Dec_region[ok]
v_region = v_region[ok]
temp_region = temp_region[ok]
ID_region = ID_region[ok]

#eliminate stars too close to M32 -- need to figure out error with below
m32 = SkyCoord(10.6708 * u.deg, 40.8653 * u.deg, frame = 'icrs')
c = SkyCoord(ra = RA_region, dec= Dec_region, unit=(u.hourangle, u.deg))
sep = m32.separation(c)
good = sep.arcsecond > 200

RA_region = RA_region[good]
Dec_region = Dec_region[good]
v_region = v_region[good]
temp_region = temp_region[good]
ID_region = ID_region[good]

#covert to xi and eta
m31 = SkyCoord(10.6847083 * u.deg, 41.26875 * u.deg, frame='icrs')
c = SkyCoord(ra = RA_region, dec= Dec_region, unit=(u.hourangle, u.deg))
c_inm31 = c.transform_to(m31.skyoffset_frame())
xi, eta=c_inm31.lon, c_inm31.lat

#eliminating data from extended area	
xi_area = np.array([a for a, b in zip(xi.value, eta.value) if a> -0.4 and b < 0.05 and b > -0.55 and b > 2 * a - 1.3])
eta_area = np.array([b for a, b in zip(xi.value, eta.value) if a> -0.4 and b < 0.05 and b > -0.55 and b > 2 * a - 1.3])
temp_area = np.array([c for a, b, c in zip(xi.value, eta.value, temp_region) if a> -0.4 and b < 0.05 and b > -0.55 and b > 2 * a - 1.3])
ID_area = np.array([c for a, b, c in zip(xi.value, eta.value, ID_region) if a> -0.4 and b < 0.05 and b > -0.55 and b > 2 * a - 1.3])
v_area = np.array([c for a, b, c in zip(xi.value, eta.value, v_region) if a> -0.4 and b < 0.05 and b > -0.55 and b > 2 * a - 1.3])

#adding in i mag data and getting rid of stars that don't have a good flag or data
i_prime, flag = np.loadtxt('/Users/amandaquirk/Documents/Ellipsoid/Data/iprime_data.txt', usecols=(1,2,), unpack = True)
err = np.loadtxt('/Users/amandaquirk/Documents/Ellipsoid/Data/iprime_errors_data.txt', usecols=(1,), unpack = True)

ok = np.isnan(i_prime) == False #all flags are greater than 0 and are ok here so not adding to statement

xi_area = xi_area[ok]
eta_area = eta_area[ok]
temp_area = temp_area[ok]
ID_area = ID_area[ok]
v_area = v_area[ok]
imag_area = i_prime[ok]
err_area = err[ok]

ok = np.isnan(err_area) == False

xi_area = xi_area[ok]
eta_area = eta_area[ok]
temp_area = temp_area[ok]
ID_area = ID_area[ok]
v_area = v_area[ok]
imag_area = imag_area[ok]
err_area = err_area[ok]

#sorting into AGB and RGB populations based on Teff and i band mag
RG_xi = np.array([a for a, b, c in zip(xi_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])  
RG_eta = np.array([a for a, b, c in zip(eta_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])
RG_v = np.array([a for a, b, c in zip(v_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])
RG_ID = np.array([a for a, b, c in zip(ID_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])
RG_err = np.array([a for a, b, c in zip(err_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])
RG_temp = np.array([b for b, c in zip(temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])
RG_i = np.array([c for b, c in zip(temp_area, imag_area) if b > 1000 and b< 10000 and c < 22 and c > 19 and c > -3/4900 * b + 1111/49])

AGy_xi = np.array([a for a, b, c in zip(xi_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476]) 
AGy_eta = np.array([a for a, b, c in zip(eta_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476])
AGy_v = np.array([a for a, b, c in zip(v_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476])
AGy_ID = np.array([a for a, b, c in zip(ID_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476])
AGy_err = np.array([a for a, b, c in zip(err_area, temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476])
AGy_temp = np.array([b for b, c in zip(temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476])
AGy_i = np.array([c for b, c in zip(temp_area, imag_area) if b > 1000 and b< 10000 and c < 20.3 and c < -0.000597561 * b + 22.4476])

AGo_xi = np.array([a for a, b, c in zip(xi_area, temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])
AGo_eta = np.array([a for a, b, c in zip(eta_area, temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])
AGo_v = np.array([a for a, b, c in zip(v_area, temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])
AGo_ID = np.array([a for a, b, c in zip(ID_area, temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])
AGo_err = np.array([a for a, b, c in zip(err_area, temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])
AGo_temp = np.array([b for b, c in zip(temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])  
AGo_i = np.array([c for b, c in zip(temp_area, imag_area) if b > 1000 and b< 10000 and c > 20.3 and c < -0.000597561 * b + 22.4476])

#plotting to check age bins
# print(len(RG_i), len(AGy_i), len(AGo_i), len(temp_area))
# plt.scatter(RG_temp, RG_i, c='r')
# plt.scatter(AGy_temp, AGy_i, c='m')
# plt.scatter(AGo_temp, AGo_i, c='k')
# plt.xlim(10000, 1000)
# plt.ylim(22, 19)
# plt.savefig('/Users/amandaquirk/Desktop/i_temp.png')
# plt.close()

# plt.scatter(RG_xi, RG_eta, c='r')
# plt.scatter(AGy_xi, AGy_eta, c='m')
# plt.scatter(AGo_xi, AGo_eta, c='k')
# plt.savefig('/Users/amandaquirk/Desktop/map.png')
# plt.close()

#smoothing velocities
#Claire does this
def adjust_err(err):
        return np.sqrt((3.2 * err)**2 + 1.2) 

#function to calculate the weights
def weights(err):
        return 1 / (err**2)

def normed_weight(w):
        sum_weights=sum(w)
        return w / sum_weights

AGy_err = adjust_err(AGy_err)
AGo_err = adjust_err(AGo_err)
RG_err = adjust_err(RG_err)

AGy_weights = weights(AGy_err)
AGo_weights = weights(AGo_err)
RG_weights = weights(RG_err)

#function does the weighted meean
def weighted_mean(data,norm_w):
	return sum([a*b for a,b in zip(data, norm_w)])

#function does the weighted RMSE
def weighted_rmse(norm_w, data, mean):
	differences=[x - mean for x in data]
	diff_sq=[d**2 for d in differences]
	products=[a*b for a,b in zip(diff_sq, norm_w)]
	return np.sqrt(sum(products))

# # m32 = SkyCoord(10.6708 * u.deg, 40.8653 * u.deg, frame = 'icrs')
# # c = SkyCoord(ra = RA_region, dec= Dec_region, unit=(u.hourangle, u.deg))
# # sep = m32.separation(c)
# # good = sep.arcsecond > 200

AGy_smoothed_v = []
AGy_xi_goodcenter = []
AGy_eta_goodcenter = []
AGy_smoothed_err = []
AGy_dispersion = []
c = SkyCoord(ra = AGy_xi, dec = AGy_eta, unit=(u.deg,u.deg))
for i in range(len(AGy_xi)):
	c1 = SkyCoord(AGy_xi[i], AGy_eta[i], unit=(u.deg,u.deg)) #go through all coordinates one at a time
	sep = c1.separation(c)
	good = sep.arcsecond < 275 #put stars into smoothing circle of this size
	velocities = AGy_v[good]
	weights = AGy_weights[good]
	if len(velocities) > 15: #only want circles with at least 15 stars
		normed_weights = normed_weight(weights)
		avg = weighted_mean(velocities, normed_weights)
		AGy_smoothed_v.append(avg) #average the velocites
		AGy_xi_goodcenter.append(AGy_xi[i] * 13.67) #kpc
		AGy_eta_goodcenter.append(AGy_eta[i] * 13.67) #kpc
		AGy_smoothed_err.append(AGy_err[i])
		AGy_dispersion.append(weighted_rmse(normed_weights, velocities, avg))

AGo_smoothed_v = []
AGo_xi_goodcenter = []
AGo_eta_goodcenter = []
AGo_smoothed_err = []
AGo_dispersion = []
c = SkyCoord(ra = AGo_xi, dec = AGo_eta, unit=(u.deg,u.deg))
for i in range(len(AGo_xi)):
	c1 = SkyCoord(AGo_xi[i], AGo_eta[i], unit=(u.deg,u.deg))
	sep = c1.separation(c)
	good = sep.arcsecond < 275
	velocities = AGo_v[good]
	weights = AGo_weights[good]
	if len(velocities) > 15:
		normed_weights = normed_weight(weights)
		avg = weighted_mean(velocities, normed_weights)
		AGo_smoothed_v.append(avg)
		AGo_xi_goodcenter.append(AGo_xi[i] * 13.67) #kpc
		AGo_eta_goodcenter.append(AGo_eta[i] * 13.67) #kpc
		AGo_smoothed_err.append(AGo_err[i])
		AGo_dispersion.append(weighted_rmse(normed_weights, velocities, avg))

RG_smoothed_v = []
RG_xi_goodcenter = []
RG_eta_goodcenter = []
RG_smoothed_err = []
RG_dispersion = []
c = SkyCoord(ra = RG_xi, dec = RG_eta, unit=(u.deg,u.deg))
for i in range(len(RG_xi)):
	c1 = SkyCoord(RG_xi[i], RG_eta[i], unit=(u.deg,u.deg))
	sep = c1.separation(c)
	good = sep.arcsecond < 200
	velocities = RG_v[good]
	weights = RG_weights[good]
	if len(velocities) > 15:
		normed_weights = normed_weight(weights)
		avg = weighted_mean(velocities, normed_weights)
		RG_smoothed_v.append(avg)
		RG_xi_goodcenter.append(RG_xi[i] * 13.67) #kpc
		RG_eta_goodcenter.append(RG_eta[i] * 13.67) #kpc
		RG_smoothed_err.append(RG_err[i])
		RG_dispersion.append(weighted_rmse(normed_weights, velocities, avg))

print(len(AGy_smoothed_v), len(AGo_smoothed_v), len(RG_smoothed_v))

f = open('/Users/amandaquirk/Desktop/AGy_additional_data.txt', 'w')
f.write('#xi (kpc), eta (kpc), avg v (km/s), dispersion (km/s)\n')
for i in range(len(AGy_smoothed_v)):
	f.write('{} {} {} {}\n'.format(AGy_xi_goodcenter[i], AGy_eta_goodcenter[i], AGy_smoothed_v[i], AGy_dispersion[i]))
f.close()

f = open('/Users/amandaquirk/Desktop/AGo_additional_data.txt', 'w')
f.write('#xi (kpc), eta (kpc), avg v (km/s), dispersion (km/s)\n')
for i in range(len(AGo_smoothed_v)):
	f.write('{} {} {} {}\n'.format(AGo_xi_goodcenter[i], AGo_eta_goodcenter[i], AGo_smoothed_v[i], AGo_dispersion[i]))
f.close()

f = open('/Users/amandaquirk/Desktop/RG_additional_data.txt', 'w')
f.write('#xi (kpc), eta (kpc), avg v (km/s), dispersion (km/s)\n')
for i in range(len(RG_smoothed_v)):
	f.write('{} {} {} {}\n'.format(RG_xi_goodcenter[i], RG_eta_goodcenter[i], RG_smoothed_v[i], RG_dispersion[i]))
f.close()