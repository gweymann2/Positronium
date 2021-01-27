import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

def gaussian_func(x, A, mu, sigma, a, b):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))+a*x+b

data= pd.read_csv("CalibrationCo.txt",delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])

##
channel_1=data.loc[(data.iloc[:,0]>=1),:]
channel_2=data.loc[(data.iloc[:,2]>=1),:]
channel_3=data.loc[(data.iloc[:,4]>=1),:]

plt.figure(figsize=(16,10))
plt.hist(channel_1['E1'],bins=100)
plt.hist(channel_2['E2'],bins=100)
plt.hist(channel_3['E3'],bins=100)
plt.show()

E1 = np.array(channel_1['E1'])
E2 = np.array(channel_2['E2'])
E3 = np.array(channel_3['E3'])

E = [E1, E2, E3]

for i in range(3):

    hist,bins = np.histogram(E[i], 5000)
    binc = 0.5*(bins[:-1]+bins[1:])
    err_hist=np.sqrt(hist)

    m = np.where(hist == np.max(hist))
    sel = (binc>(binc[m]*0.93))*(binc<(binc[m]*1.07))
    mean = np.mean(binc[sel])

    par, par_var=curve_fit(gaussian_func, binc[sel], hist[sel], p0=[np.max(hist), mean, mean*0.01, 1, 1],sigma=err_hist[sel],absolute_sigma=True)
    sel_off = (binc>(par[1]-3*par[2]))*(binc<(par[1]+3*par[2]))
#    mean = np.mean(binc[sel_off])

    par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off], hist[sel_off], p0=par, sigma=err_hist[sel_off],absolute_sigma=True)

    
    plt.figure(i, figsize=(16,10))
    plt.scatter(binc[sel_off],hist[sel_off])
    plt.plot(binc[sel_off], gaussian_func(binc[sel_off], *par_off))
    plt.savefig('FitPeak_channel{}'.format(i))


    
#    sel_off2 = (binc>(par[1]+3*par_off[2]))
#    m = np.where(hist == np.max(hist[sel_off2]))
#    sel_new = (binc>(binc[m]*0.93))*(binc<(binc[m]*1.07))
#    mean_new = np.mean(binc[sel_new])
#
#
#    par_new, par_var_new=curve_fit(gaussian_func, binc[sel_new], hist[sel_new], p0=[np.max(hist[sel_off2]), mean_new, mean_new*0.01, 1, 1],sigma=err_hist[sel_new],absolute_sigma=True)
#    sel_off_new = (binc>(par_new[1]-3*par_new[2]))*(binc<(par_new[1]+3*par_new[2]))
##    mean = np.mean(binc[sel_off_new])
##    print(mean)
#    print(binc[sel_off_new], hist[sel_off_new])
#    print(par_new)
#
#    par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off_new], hist[sel_off_new], p0=par_new, sigma=err_hist[sel_off_new],absolute_sigma=True)
#
#    plt.figure(figsize=(16,10))
#    plt.scatter(binc[sel_off_new],hist[sel_off_new])
#    plt.plot(binc[sel_off_new], gaussian_func(binc[sel_off_new], *par_off))
#    plt.savefig('FitPeak2_channel{}'.format(i))

    

#Cuts on the arrays, last peak

#sel = (binc>620000)*(binc<700000)

#mean = np.mean(binc[sel])

#par, par_var=curve_fit(gaussian_func, binc[sel], hist[sel], p0=[1, mean, mean*0.01, 1, 1],sigma=err_hist[sel],absolute_sigma=True)

#print(*par)
#print(np.sqrt(np.diag(par_var)))
#plt.figure(figsize=(16,10))
#plt.scatter(binc[sel],hist[sel])
#plt.plot(binc[sel], gaussian_func(binc[sel], *par))
#
#sel_off = (binc>(par[1]-3*par[2]))*(binc<(par[1]+3*par[2]))
#mean = np.mean(binc[sel_off])
#
#par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off], hist[sel_off], p0=par, sigma=err_hist[sel_off],absolute_sigma=True)
#
#plt.figure(figsize=(16,10))
#plt.scatter(binc[sel_off],hist[sel_off])
#plt.plot(binc[sel_off], gaussian_func(binc[sel_off], *par_off))
#plt.show()
