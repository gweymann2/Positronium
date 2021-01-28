import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas as pd

def Cal(E,a):
    i=a*E
    return(i)

def Reso(E,k):
    i=2.35*np.sqrt(k/E)
    return(i)

def Calibration(ADC,err_ADC,E_exp,Channel_number):  ##Channel number is a word; ex: "channel 2"
    """Calibration fit C=f(E)"""

    par_Cal, par_Cal_var=curve_fit(Cal,E_exp,ADC,sigma=err_ADC,absolute_sigma=True)

    a=round(par_Cal[0],3)                   ### a coef to go like C=a*E
    err_a=round(np.sqrt(np.diag(par_Cal_var)[0]),3)


    #plt.figure(Channel_number)
    plt.errorbar(E_exp,ADC,yerr=err_ADC,fmt='b.',capsize=4)
    plt.plot(E_exp,Cal(E_exp,par_Cal[0]),c='k',label="a=%s"%a+"$\pm$%s"%err_a)
    plt.title('Calibration detector{} ADC=f(E)'.format(Channel_number),fontsize=20)
    plt.xlabel('E (keV)')
    plt.ylabel('ADC')
    # plt.annotate(r'$\a$',xy=(1000,np.mean(ADC)), xycoords='data',xytext=(+10, +30), textcoords='offset points', fontsize=16)
    # plt.annotate(f'={round(a,3)}+-{round(err_a,3)}',xy=(1200,np.mean(ADC)), xycoords='data',xytext=(+10, +30), textcoords='offset points', fontsize=16)
    plt.legend()
    plt.savefig('Calibration_channel_{}'.format(Channel_number))
    plt.show()

    return ([a,err_a])

def Resolution(ADC,err_ADC,ADC_Width,err_ADC_width,a,err_a,Channel_number):    ##a is the coeff ADC=a*E so in kev-1 and Channel number is a word; ex: "channel 2"

    """Resolution fit"""
    ###First define E and E_width with their respective uncertainties
    E=ADC/a
    err_E=np.sqrt((err_ADC/a)**2+(ADC*err_a/a**2)**2)

    E_Width=ADC_width/a
    err_E_Width=np.sqrt((err_ADC_Width/a)**2+(ADC_Width*err_a/a**2)**2)

    ###Define the resolution as R
    R=E_Width/E
    err_R=np.sqrt((err_E_Width/E)**2+(err_E*E_Width/E**2)**2)

    par_Res, par_Res_var=curve_fit(Res,E,R,sigma=err_R,absolute_sigma=True)

    #plt.figure(Channel_number)
    plt.errorbar(E,R,xerr=err_E,yerr=err_R,fmt='b.',capsize=4)
    plt.plot(E,Res(np.linspace(min(E),max(E),1000),par_Res[0]),c='k',label="ajustement function")
    plt.title('Resolution R=f(E)',fontsize=20)
    plt.xlabel('E (keV)')
    plt.ylabel('R')
    plt.legend()
    plt.show()
    plt.savefig('Resolution_channel_{}'.format(Channel_number))

    k=par_Res[0]**2
    err_k=np.sqrt(np.diag(par_Res_var)[0])*par_Res[0]*np.sqrt(2)

    return(Channel_number," k=", k,"+-",err_k,"keV")




Elements=['Co', 'Na', 'Cs']
Channels=['Ch1', 'Ch2', 'Ch3']

##creating empty dataframes to regroup the Mean and width values of the pics by channels instead of sources
Ch1=pd.DataFrame(columns=['Channel','Amplitude','Average', 'AverageError','Sigma', 'SigmaError','a','b'])
Ch2=pd.DataFrame(columns=['Channel','Amplitude','Average', 'AverageError','Sigma', 'SigmaError','a','b'])
Ch3=pd.DataFrame(columns=['Channel','Amplitude','Average','AverageError','Sigma', 'SigmaError','a','b'])

for el in Elements:
    data= pd.read_csv('Parameters'+el+'.csv',delimiter=" ")

    data1=data.loc[(data['Channel']==1),:]  ##"2lines" dataframe conrresponding to the channel1 (because 2 pics )
    data2=data.loc[(data['Channel']==2),:] ##"2lines" dataframe conrresponding to the channel2
    data3=data.loc[(data['Channel']==3),:] ##"2lines" dataframe conrresponding to the channel3

    Ch1=pd.concat(([Ch1,data1]))
    Ch2=pd.concat(([Ch2,data2]))
    Ch3=pd.concat(([Ch3,data3]))

Ch1=(np.array(Ch1)).astype(float)   ##Vectorizing the dataframes + making sure the elements are floats
Ch2=(np.array(Ch2)).astype(float)
Ch3=(np.array(Ch3)).astype(float)


E_exp=np.array([1173.,1332.,511.,1275.,662.]) ##expected energy of the pics from the spectrums Co Na Cs in keV

print(Calibration(Ch1[:,2],Ch1[:,3],E_exp,1))
print(Calibration(Ch2[:,2],Ch2[:,3],E_exp,2))
print(Calibration(Ch3[:,2],Ch3[:,3],E_exp,3))

# Resolution(Ch1[:,2],Ch1[:,3],Ch1[:,4],Ch1[:,5],a,err_a,Channel_number)
