import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas

def gaussian_func(x, A, mu, sigma, a, b):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))+a*x+b

Elements=['Co', 'Na', 'Cs']
Channel=['ch1', 'ch2', 'ch3']
Title =['Channel 1: BERNA', 'Channel 2: DEMO', 'Channel 3: BLACK']

for k,el in enumerate(Elements):
    data= pd.read_csv('../muoni/Calibration'+el+'.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])
    #data= pd.read_csv('../muoni/CalibrationCo.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])

    ##
    channel_1=data.loc[(data.iloc[:,0]>=1),:]
    channel_2=data.loc[(data.iloc[:,2]>=1),:]
    channel_3=data.loc[(data.iloc[:,4]>=1),:]

    plt.figure(10*k,figsize=(16,10))
    plt.hist(channel_1['E1'],bins=1000,label=Title[0])
    plt.hist(channel_3['E3'],bins=1000,label=Title[1])
    plt.hist(channel_2['E2'],bins=1000,label=Title[2])
    plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
    plt.title('Histograms for each channels '+el,fontsize=17)
    plt.xlabel('ADC',fontsize=17)
    plt.ylabel('Count of photons',fontsize=17)
    for tickLabel in plt.gca().get_xticklabels()+plt.gca().get_yticklabels():
        tickLabel.set_fontsize(15)
    #plt.show()
    plt.savefig('FullSpectra_'+el)

    E1 = np.array(channel_1['E1'])
    E2 = np.array(channel_2['E2'])
    E3 = np.array(channel_3['E3'])

    E = [E1, E2, E3]
    parameters1 = [['Channel', 'Amplitude', 'Average', 'Sigma', 'a', 'b']]
    #parameters2 = [['Channel', 'Amplitude', 'Average', 'Sigma', 'a', 'b']]
    parameters2 = []

    for i,ch in enumerate(Channel):

        hist,bins = np.histogram(E[i], 1000)
        binc = 0.5*(bins[:-1]+bins[1:])
        err_hist=np.sqrt(hist)
        conr = (binc>2*10**5)
        m = np.where(hist == np.max(hist[conr]))
        sel = (binc>(binc[m]*0.93))*(binc<(binc[m]*1.07))
        mean = np.mean(binc[sel])

        par, par_var=curve_fit(gaussian_func, binc[sel], hist[sel], p0=[np.max(hist), mean, mean*0.01, 1, 1],sigma=err_hist[sel],absolute_sigma=True)
        sel_off = (binc>(par[1]-3*par[2]))*(binc<(par[1]+3*par[2]))
    #    mean = np.mean(binc[sel_off])

        par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off], hist[sel_off], p0=par, sigma=err_hist[sel_off],absolute_sigma=True)

        parameters1.append([i+1, par_off[0],par_off[1],par_off[2],par_off[3],par_off[4]])

        plt.figure(10*k+i+1, figsize=(16,10))
        g1=plt.scatter(binc[sel_off],hist[sel_off],label='1st peak',color='lightblue')
        plt.errorbar(binc[sel_off], hist[sel_off], err_hist[sel_off],fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)
        g2=plt.plot(binc[sel_off], gaussian_func(binc[sel_off], *par_off),label='1st peak fit',color='lightblue')

        if (el=='Co'):

            sel_off2 = (binc>(par[1]+3*par_off[2]))
            m = np.where(hist == np.max(hist[sel_off2]))
            sel_new = (binc>(binc[m]*0.93))*(binc<(binc[m]*1.07))
            mean_new = np.mean(binc[sel_new])


            par_new, par_var_new=curve_fit(gaussian_func, binc[sel_new], hist[sel_new], p0=[np.max(hist[sel_off2]), mean_new, mean_new*0.01, 1, 1],sigma=err_hist[sel_new],absolute_sigma=True)
            sel_off_new = (binc>(par_new[1]-3*par_new[2]))*(binc<(par_new[1]+3*par_new[2]))
        #    mean = np.mean(binc[sel_off_new])
        #    print(mean)
        #    print(binc[sel_off_new], hist[sel_off_new])
        #    print(par_new)

            par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off_new], hist[sel_off_new], p0=par_new, sigma=err_hist[sel_off_new],absolute_sigma=True)

            parameters2.append([i+1, par_off[0], par_off[1],par_off[2],par_off[3],par_off[4]])

            g3=plt.scatter(binc[sel_off_new],hist[sel_off_new],label='2nd peak',color='orange')
            plt.errorbar(binc[sel_off_new], hist[sel_off_new], err_hist[sel_off_new],fmt='.',color='orange',ecolor='lightgray', elinewidth=2, capsize=0)
            g4=plt.plot(binc[sel_off_new], gaussian_func(binc[sel_off_new], *par_off),label='2nd peak fit',color='orange')

        plt.title(Title[i]+' '+el,fontsize=17)
        plt.xlabel('ADC',fontsize=17)
        plt.ylabel('Count of photons',fontsize=17)
        for tickLabel in plt.gca().get_xticklabels()+plt.gca().get_yticklabels():
            tickLabel.set_fontsize(15)
        plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
        plt.savefig('FitPeak2_channel_'+el+ch.format(i))
    #    plt.show()
    paramTot = parameters1 + parameters2

    np.savetxt("Parameters"+el+".csv", paramTot, delimiter=" ", fmt='%s')

    #parameters1_table = AsciiTable(parameters1)
    #parameters2_table = AsciiTable(parameters2)
    #print(parameters1_table.table)
    #print(parameters2_table.table)
