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

def chi2red(sel,are1or2gaussians,*p):
    if are1or2gaussians==1:
        chi2=sum(((hist[sel] - gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
        dof=len(binc[sel])-5
    else:
        chi2=sum(((hist[sel] - gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
        dof=len(binc[sel])-8
    #print(chi2, dof)
    return round(chi2/dof,2)


ch = ['12','23','13']
Title =['Channel 1: BERNA', 'Channel 2: DEMO', 'Channel 3: BLACK']
energy = ['E1','E2','E3']
coeff = [495.3662907422378, 499.20938900287547, 457.96722725099886]
data=[]
en = [258984.8434375515, 258388.2395968352, 238664.36333105486]
devEn = [115.89422582052933, 140.38598828118543, 112.05306994484995]

for el in ch:
    data.append(pd.read_csv('../Coinc'+el+'noDel.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4']))

coinc12 = pd.read_csv('../Coinc12del1.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])
coinc13 = pd.read_csv('../Coinc13del1.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])
coinc23 = pd.read_csv('../Coinc23del2.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])

dataD = [coinc12,coinc13,coinc23]

# plt.figure(figsize=(8,5))
# plt.scatter(data[0]['t1'],data[0]['E1'],label='1st peak',color='lightblue')
# plt.show()

for n,j in enumerate(energy):
    hist,bins = np.histogram(data[n][j], 2000)
    m = np.where(hist == np.max(hist))

    data[n]=data[n].loc[(data[n][j]>(bins[m][0]-0.1*bins[m][0]))&(data[n][j]<(bins[m][0]+0.1*bins[m][0])),:]

t12 = data[0]['t1'] - data[0]['t2']
t23 = data[1]['t2'] - data[1]['t3']
t13 = data[2]['t1'] - data[2]['t3']



dataT = [t12, t23, t13]

for n,j in enumerate(energy):
    histT,binsT = np.histogram(dataT[n], 1000)

    bincT = 0.5*(binsT[:-1]+binsT[1:])

    # plt.figure(figsize=(8,5))
    # plt.hist(dataT[n],bins=100,label='Title')
    # plt.show()
    # plt.hist(data[n][j],bins=100,label='Title')

    ma = np.where(histT == np.max(histT))
    # print(bincT[ma])
    defin = (bincT>(bincT[ma][0]-np.abs(1.1*bincT[ma][0])))*(bincT<(bincT[ma][0]+np.abs(1.1*bincT[ma][0])))
    #defin = (bincT>(bincT[ma][0]-0.3*np.abs(bincT[ma][0])))*(bincT<(bincT[ma][0]+0.3*np.abs(bincT[ma][0])))


    histTcut = np.delete(histT[defin], np.where(histT[defin] == 0))
    bincTcut = np.delete(bincT[defin], np.where(histT[defin] == 0))
    # print(bincTcut)
    # print(histTcut)
    # plt.figure(figsize=(8,5))
    # plt.scatter(bincTcut,histTcut,label='1st peak',color='lightblue')
    # plt.show()
    err_histT=np.sqrt(histTcut)

    par, par_var = curve_fit(gaussian_func, bincTcut, histTcut,p0=[np.max(histTcut),bincT[ma][0],bincT[ma][0]],sigma=err_histT,absolute_sigma=True)
    # print(par)
    sel_off = (bincT>(par[1]-3*np.abs(par[2])))*(bincT<(par[1]+3*np.abs(par[2])))
    #print(bincT[sel_off])
    par_new, par_var_new = curve_fit(gaussian_func, bincTcut, histTcut, p0=par, sigma=err_histT,absolute_sigma=True)

    # chi2=chi2red(sel_off,2,*par)
    # a=chi2**0.5
    # print(a)
    # #
    # par_off, par_var_off = curve_fit(gaussian_func, bincTcut, histTcut, p0=par_new, sigma=a*err_histT,absolute_sigma=True)

    mu_err = np.sqrt(np.diag(par_var_new)[1])
    sigma_err = np.sqrt(np.diag(par_var_new)[2])

    print(n+1, par_new[0],par_new[1], mu_err, par_new[2], sigma_err)

    plt.figure(10*n+1, figsize=(8,5))
    g1=plt.scatter(bincTcut,histTcut,label='1st peak',color='lightblue')
    plt.errorbar(bincTcut, histTcut, err_histT,fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)

    g2=plt.plot(np.linspace(bincT[ma][0]-np.abs(bincT[ma][0]),bincT[ma][0]+np.abs(bincT[ma][0]),1000), gaussian_func(np.linspace(bincT[ma][0]-np.abs(bincT[ma][0]),bincT[ma][0]+np.abs(bincT[ma][0]),1000), *par),label='peak fit',color='lightblue')

    plt.title('Coincidence {}'.format(j),fontsize=17)
    plt.xlabel('Time [s]',fontsize=17)
    plt.ylabel('Count of photons',fontsize=17)
    plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
    plt.savefig('Coincidence_'+ch[n]+'.png')

    #data[n]=data[n].loc[(data[n][j]>(bins[m][0]-0.1*bins[m][0]))&(data[n][j]<(bins[m][0]+0.1*bins[m][0])),:]
#     data[n][j] = data[n][j]/coeff[n]
#print(data[1])

for n,j in enumerate(energy):
    histD,binsD = np.histogram(dataD[n][j], 1000)
    mD = np.where(histD == np.max(histD))

    dataD[n]=dataD[n].loc[(dataD[n][j]>(binsD[mD][0]-0.1*binsD[mD][0]))&(dataD[n][j]<(binsD[mD][0]+0.1*binsD[mD][0])),:]

T12 = t12 - t12D
T13 = t13 - t13D
T23 = t23 - t23D

plt.figure(figsize=(8,5))
#plt.hist(t12,bins=100,label='Title')
plt.hist(t23,bins=100,label='Title')
plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
plt.title('Histograms for each channels ',fontsize=17)
plt.xlabel('Time [ns]',fontsize=17)
plt.ylabel('Events',fontsize=17)
plt.show()

#print(data[0])
#data=data.loc[(data['Channel']==1),:]
