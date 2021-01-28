import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas

def chi2red(sel,*p):
    chi2=sum(((hist[sel] - gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
    dof=len(binc[sel])-5
    #print(chi2, dof)
    return round(chi2/dof,2)

def gaussian_func(x, A, mu, sigma, a, b):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))+a*x+b

def _2_gaussian_func(x, A1, mu1, sigma1, A2, mu2, sigma2, a, b):
    return A1 * np.exp( - (x - mu1)**2 / (2 * sigma1**2))+a*x+b+A2 * np.exp( - (x - mu2)**2 / (2 * sigma2**2))

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
    parameters1 = [['Channel', 'Amplitude', 'Average', 'AverageError', 'Sigma', 'SigmaError', 'a', 'b']]
    #parameters2 = [['Channel', 'Amplitude', 'Average', 'Sigma', 'a', 'b']]
    parameters2 = []

    for i,ch in enumerate(Channel):

        hist,bins = np.histogram(E[i], 1000)
        binc = 0.5*(bins[:-1]+bins[1:])
        err_hist=np.sqrt(hist)
        conr = (binc>2*10**5)
        m = np.where(hist == np.max(hist[conr]))
        if (el == 'Co'):
            sel_2max = (binc>(binc[m]*1.07))
            maxim = np.where(hist == np.max(hist[sel_2max]))
            defin = (binc>(binc[m]*0.93))*(binc<binc[maxim]*1.07)
            mean_2max = np.mean(binc[sel_2max])
            #print(binc[m], binc[maxim])
            par_2max, par_va_2max = curve_fit(_2_gaussian_func, binc[defin], hist[defin], p0=[np.max(hist), binc[m][0], binc[m][0]*0.01, np.max(hist), binc[maxim][0], binc[maxim][0]*0.01, 1, 1 ],sigma=err_hist[defin],absolute_sigma=True)
            #print(*par_2max)
            #sel_off_2max = (binc>(par_2max[1]-3*par_2max[2]))*(binc<(par_2max[4]+3*par_2max[5]))
            #par_2max_off, par_va_2max_off = curve_fit(_2_gaussian_func, binc[defin], hist[defin])

            plt.figure(10*k+i+1, figsize=(16,10))
            g1=plt.scatter(binc[defin],hist[defin],label='1st peak',color='lightblue')
            plt.errorbar(binc[defin], hist[defin], err_hist[defin],fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)
            g2=plt.plot(binc[defin], _2_gaussian_func(binc[defin], *par_2max),label='peak fit'+' ($\chi^2$/dof = %s'+')',color='lightblue')

        if (el == 'Na' or 'Cs'):
            sel = (binc>(binc[m]*0.93))*(binc<(binc[m]*1.07))
            mean = np.mean(binc[sel])

            par, par_var=curve_fit(gaussian_func, binc[sel], hist[sel], p0=[np.max(hist), mean, mean*0.01, 1, 1],sigma=err_hist[sel],absolute_sigma=True)
            sel_off = (binc>(par[1]-3*par[2]))*(binc<(par[1]+3*par[2]))
        #    mean = np.mean(binc[sel_off])

            par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off], hist[sel_off], p0=par, sigma=err_hist[sel_off],absolute_sigma=True)
            mu_err = np.sqrt(np.diag(par_var_off)[1])
            sigma_err = np.sqrt(np.diag(par_var_off)[2])
            parameters1.append([i+1, par_off[0],par_off[1], mu_err, par_off[2], sigma_err, par_off[3],par_off[4]])

            plt.figure(10*k+i+1, figsize=(16,10))
            g1=plt.scatter(binc[sel_off],hist[sel_off],label='1st peak',color='lightblue')
            plt.errorbar(binc[sel_off], hist[sel_off], err_hist[sel_off],fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)
            g2=plt.plot(binc[sel_off], gaussian_func(binc[sel_off], *par_off),label='1st peak fit'+' ($\chi^2$/dof = %s'%chi2red(sel_off,*par_off)+')',color='lightblue')

            if (el=='Na'): # or 'Na'

                sel_off2 = (binc>(par[1]+3*par_off[2]))
                ma = np.where(hist == np.max(hist[sel_off2]))
                c=binc[ma]
                d=np.max(c)
                # print(binc[sel_off2])
                # print(c)
                # print(type(c))
                sel_new = (binc>(d*0.90))*(binc<(d*1.10))
                mean_new = np.mean(binc[sel_new])


                par_new, par_var_new=curve_fit(gaussian_func, binc[sel_new], hist[sel_new], p0=[np.max(hist[sel_off2]), mean_new, mean_new*0.01, 1, 1],sigma=err_hist[sel_new],absolute_sigma=True)
                sel_off_new = (binc>(par_new[1]-3*par_new[2]))*(binc<(par_new[1]+3*par_new[2]))

                par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off_new], hist[sel_off_new], p0=par_new, sigma=err_hist[sel_off_new],absolute_sigma=True)
                mu_err = np.sqrt(np.diag(par_var_off)[1])
                sigma_err = np.sqrt(np.diag(par_var_off)[2])
                parameters2.append([i+1, par_off[0], par_off[1], mu_err, par_off[2], sigma_err, par_off[3],par_off[4]])

                g3=plt.scatter(binc[sel_off_new],hist[sel_off_new],label='2nd peak',color='orange')
                plt.errorbar(binc[sel_off_new], hist[sel_off_new], err_hist[sel_off_new],fmt='.',color='orange',ecolor='lightgray', elinewidth=2, capsize=0)
                g4=plt.plot(binc[sel_off_new], gaussian_func(binc[sel_off_new], *par_off),label='2nd peak fit'+' ($\chi^2$/dof = %s'%chi2red(sel_off,*par_off)+')',color='orange')

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
