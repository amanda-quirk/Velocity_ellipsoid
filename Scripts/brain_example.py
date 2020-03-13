 #posterior for simultaneous fit of many galaxies, called by simulfit
def slogpost(thetasimul, vfs, xs, ys, ivars, proj, gxs, papnsas, shears):
  #AQ thetasimul is the array of parameters to fit
   #0=vmax, 1=inc, 2=pak, 3=hrot, 4=gx, 5 = papnsa
   #lambda value (shear scaling hyperparameter)
   l = thetasimul[-1] #AQ this is Brian's hyperparameter
#uniform prior on l 
   if l < 0 or l > 1:
       return -np.inf

#AQ i think the below section is not needed for me =============================
   #assign predicted shear values, kept consistent
   if type(shears) != bool:
       thetasimul[4::6] = shears

   #parameter for choosing between lambda being a continuous normalization
   #or a discrete model preference
   discrete = True
   if discrete:
       l = int(round(l,0))

       #set shear to 0 if l below .5, otherwise leave it
       if l == 0:
           thetasimul[4::6] = 0
       else:
          pass

   #scale shear by l
   else:
       thetasimul[4::6] = thetasimul[4::6] * l
#AQ===============================================================================

   #pick out all params for each individual galaxy -- AQ: these are the individual parameters
   vmaxi = thetasimul[0::6][:-1]
   inci = thetasimul[1::6]# % (2*np.pi)
   paki = thetasimul[2::6]
   hroti = thetasimul[3::6]
   gxi = thetasimul[4::6]
   papnsai = thetasimul[5::6]

   #cumulative values for all galaxies -- AQ starting them at 0
   lprior = 0
   llike = 0
   good = 0

   #iterate through all galaxies and run prior and like on all of them
   #AQ: where is l in the below loop?
   for i in range(len(vmaxi)): #AQ: Brian says this for loop is what makes it hierarchical
       theta = (vmaxi[i], inci[i], paki[i], hroti[i], papnsai[i], gxi[i]) #AQ individual parameters; sets priors
       lp = hlogprior(theta, papnsas[i]) #AQ is this function from emcee or Brian's own
       ll = vflike(theta, vfs[i], xs[i], ys[i], ivars[i], proj) #AQ: I am guessing this is Brian's function
       #if not np.isfinite(lp):
       #    #print('prior', i, end = ' ')
       #    good += 1
       #    continue
       #if not np.isfinite(ll):
       #    #print('like', i, end = ' ')
       #    good += 1
       #    continue
       lprior += lp
       llike += ll

   #if good != 0:
   #    print(good)
   #    print('\n')
   return lprior + llike