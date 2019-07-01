import numpy as np 
import matplotlib.pyplot as plt 

def sine(angle):
	return np.sin(np.deg2rad(angle))

def cosine(angle):
	return np.cos(np.deg2rad(angle))

def sigma_LOS(oR, oZ, oP, i, PA): #i and PA come from tilted ring table (any dispersion values?)
	oLOS2 = (oR**2 * sine(PA)**2 + oP**2 * cosine(PA)**2) * sine(i)**2 + oZ**2 * cosine(i)**2
	return np.sqrt(oLOS2)

def va_term(oR, oP, oZ, R, vc): #vc and R from titled ring table
	Rd = 5.76 #kpc
	if R < 10:
		k = 2
	else:
		k = 1
	return oR**2 / (2 * vc) * (oP**2 / oR**2 - 1.5 + k * R / Rd + oZ**2 / (2 * oR**2))

AR1 = .8
AR2 = np.sqrt(.2)

#going to make 50 random oR terms and get oZ and oP from that
fake_oR = np.random.normal(80, 30, 50)
fake_oR = [abs(a) for a in fake_oR] #don't want any negative dispersions

#make factor so that the axis ratio is obsurced
AR1_wiggle = np.random.normal(3, .5, 50) 
AR1_wiggle = [abs(a) for a in AR1_wiggle]
fake_oP = fake_oR * np.array(AR1_wiggle) * AR1
plt.clf()
plt.scatter(fake_oR, fake_oP, c='k')
plt.plot(np.linspace(min(fake_oR), max(fake_oR)), AR1 * np.linspace(min(fake_oR), max(fake_oR)), c='b')
plt.show() 