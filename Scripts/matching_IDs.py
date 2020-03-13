import numpy as np 
from astropy.io import fits 

ID_data = np.genfromtxt('IDs.txt')
# ID_imag, imag, flag = np.loadtxt('i_prime.txt', unpack=True)

# ID = []
# i_prime = []
# i_flag = []
# for i in range(len(ID_data)):
# 	if ID_data[i] in ID_imag:
# 		n=np.where(ID_data[i]==ID_imag)
# 		N=n[0]
# 		ID.append(int(ID_data[i]))
# 		i_prime.append(imag[n])
# 		i_flag.append(flag[n])
# 	else:
# 		ID.append(ID_data[i])
# 		i_prime.append(np.nan)
# 		i_flag.append(np.nan)

# ID_imag, imag, flag = np.loadtxt('i_prime.txt', unpack=True)

ID_C, err = np.genfromtxt('claire_errors_data.txt', unpack=True)

ID = []
error = []
for i in range(len(ID_data)):
	if ID_data[i] in ID_C:
		n=np.where(ID_data[i]==ID_C)
		N=n[0]
		ID.append(int(ID_data[i]))
		error.append(err[n])
	else:
		ID.append(ID_data[i])
		error.append(np.nan)

print(len(ID_data), len(ID))

f = open('iprime_errors_data.txt', 'w')
for i in range(len(ID)):
	f.write('{} {}\n'.format(ID[i], error[i]))
f.close()

