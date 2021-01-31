import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas

def gaussian_func(x, A, mu, sigma):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))

x = [-74.4, -67.31428571, -53.14285714, -46.05714286, -38.97142857, -31.88571429, -24.8, -17.71428571, -10.62857143]
y = [175, 401, 1187, 2700, 4315, 4181, 2516, 876, 161 ]
err = np.sqrt(y)

par, par_var = curve_fit(gaussian_func, x, y,p0=[4300,-40,20],sigma=err,absolute_sigma=True)
print(par)

plt.figure(figsize=(8,5))
#plt.hist(t12,bins=100,label='Title')
plt.scatter(x,y,label='Peak',color='orange')
plt.errorbar(x, y, err,fmt='.',color='orange',ecolor='lightgray', elinewidth=2, capsize=0)
plt.plot(np.linspace(-80,0,1000), gaussian_func(np.linspace(-80,0,1000), *par),label='Peak fit',color='orange')


plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
plt.title('Histograms for each channels ',fontsize=17)
plt.xlabel('Time [ns]',fontsize=17)
plt.ylabel('Events',fontsize=17)
plt.show()
# sel_off = (x>(par[1]-3*np.abs(par[2])))*(x<(par[1]+3*np.abs(par[2])))
# par_new, par_var_new = curve_fit(gaussian_func, x[sel_off], y[sel_off], p0=par, sigma=err[sel_off],absolute_sigma=True)
# print(par_new)
