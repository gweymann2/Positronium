import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas

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

for n,j in enumerate(energy):
    hist,bins = np.histogram(data[n][j], 2000)
    m = np.where(hist == np.max(hist))
    print(data[n])
    data[n]=data[n].loc[(data[n][j]>(bins[m][0]-0.1*bins[m][0]))&(data[n][j]<(bins[m][0]+0.1*bins[m][0])),:]

    plt.figure(figsize=(8,5))
    print(data[n])
    #plt.scatter(bins,hist,label='1st peak',color='lightblue')


#     data[n][j] = data[n][j]/coeff[n]
#print(data[1])

for n,j in enumerate(energy):
    histD,binsD = np.histogram(dataD[n][j], 1000)
    mD = np.where(histD == np.max(histD))

    dataD[n]=dataD[n].loc[(dataD[n][j]>(binsD[mD][0]-0.1*binsD[mD][0]))&(dataD[n][j]<(binsD[mD][0]+0.1*binsD[mD][0])),:]

t12 = data[0]['t1'] - data[0]['t2']
t23 = data[1]['t2'] - data[1]['t3']
t13 = data[2]['t1'] - data[2]['t3']

t12D = dataD[0]['t1'] - dataD[0]['t2']
t23D = dataD[1]['t2'] - dataD[1]['t3']
t13D = dataD[2]['t1'] - dataD[2]['t3']

T12 = t12 - t12D
T13 = t13 - t13D
T23 = t23 - t23D

print(t12)

plt.figure(figsize=(8,5))
plt.hist(t12,bins=1500,label='Title')
#plt.hist(t23,bins=100,label='Title')
plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
plt.title('Histograms for each channels ',fontsize=17)
plt.xlabel('t',fontsize=17)
plt.ylabel('Count of photons',fontsize=17)
plt.show()
