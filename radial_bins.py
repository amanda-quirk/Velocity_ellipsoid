import numpy as np 
import matplotlib.pyplot as plt 

'''
want to find how many stars are in various radial bin
'''

#data
MS_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/MS_pixel_vrot.txt', usecols=(0,), unpack=True)
AGy_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGy_pixel_vrot.txt', usecols=(0,), unpack=True)
AGo_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/AGo_pixel_vrot.txt', usecols=(0,), unpack=True)
RG_r=np.loadtxt('/Users/amandaquirk/Documents/AsymmetricDrift/Data/RG_pixel_vrot.txt', usecols=(0,), unpack=True)

#parameters
Rmin=min(min(MS_r), min(AGy_r), min(AGo_r), min(RG_r)) #kpc
Rmax=max(max(MS_r), max(AGy_r), max(AGo_r), max(RG_r)) + 1
delta_r= 0.9 #kpc

#making radial bins
r_bins=np.linspace(Rmin, Rmax, Rmax/delta_r)

#will contain number of stars in each bin
MS_num=np.zeros(len(r_bins))
AGy_num=np.zeros(len(r_bins))
AGo_num=np.zeros(len(r_bins))
RG_num=np.zeros(len(r_bins))

#dividing each age groups into radial bin
for i in range(len(r_bins)-1):
	MS_stars=[a for a in MS_r if a>=r_bins[i] and a<r_bins[i+1]]
	MS_num[i]=len(MS_stars)
	AGy_stars=[a for a in AGy_r if a>=r_bins[i] and a<r_bins[i+1]]
	AGy_num[i]=len(AGy_stars)
	AGo_stars=[a for a in AGo_r if a>=r_bins[i] and a<r_bins[i+1]]
	AGo_num[i]=len(AGo_stars)
	RG_stars=[a for a in RG_r if a>=r_bins[i] and a<r_bins[i+1]]
	RG_num[i]=len(RG_stars)

#array of all stars
total_num =[a + b + c + d for a, b, c, d in zip(MS_num, AGy_num, AGo_num, RG_num)]
young_num = [a + b + c for a, b, c, in zip(MS_num, AGy_num, AGo_num)]
#plotting
print(sum(MS_num))
print(sum(AGy_num))
print(sum(AGo_num))
print(sum(RG_num))
#plt.scatter(r_bins[:-1], MS_num[:-1], c='b',  s=4, label='MS')
#plt.scatter(r_bins[:-1], AGy_num[:-1], c='m', s=4, label='young AGB')
#plt.scatter(r_bins[:-1], AGo_num[:-1], c='k', s=4, label='older AGB')
plt.scatter(r_bins[:-1], young_num[:-1], c='b',  s=4, label='MS + AGB')
plt.scatter(r_bins[:-1], RG_num[:-1], c='r',  s=4, label='RGB')
#plt.scatter(r_bins[:-1], total_num[:-1], c='darkgrey',  s=4, label='total')
plt.plot([0,30], [100, 100], linestyle='--', c='k')
plt.xlim(0,30)
plt.ylabel('Number of Stars in R Bin')
plt.xlabel('Radius (kpc), bin size= {} kpc'.format(delta_r))
plt.legend(frameon=False)
plt.savefig('/Users/amandaquirk/Desktop/binning_{}.png'.format(delta_r))

#determining how many bins are needed to have 100 stars in each
MS_r_order=sorted(MS_r)
AGy_r_order=sorted(AGy_r)
AGo_r_order=sorted(AGo_r)
RG_r_order=sorted(RG_r)

# total_r = np.concatenate((MS_r, AGy_r, AGo_r, RG_r))
# total_r_order = sorted(total_r)

# i=0
# while i<len(total_r_order):
# 	print(total_r_order[i])

# 	i=i+100






	
