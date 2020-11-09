#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 23:49:34 2020

Reads in some InSAR and model data and makes some plots..

@author: mlees
"""

import sys
import distutils.util as ut

if len(sys.argv) <= 1:
    print('compare_with_deformation_data error; terminal. Incorrect number of input arguments. Correct usage: python compare_with_deformation_data.py model_directory')
    sys.exit(1)

directory=sys.argv[1]

extraloc = [var for var in sys.argv[1:] if var.split('=')[0]=='extraloc']
if not extraloc:
    extraloc=None
else:
    extraloc=bool(ut.strtobool(extraloc[0].split('=')[1]))
    print('\textraloc specified; also plotting the periodic bullseye.')


foot_to_cm = 30.48

import os
cwd = os.getcwd()
os.chdir(directory)


import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
import seaborn as sns
from InSAR_postSBAS import *

if ospath.exists('/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    mac=1
    linux=0
    print("\tWe're in Mac, looking for downloaded InSAR data accordingly.")
    file='/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
elif ospath.exists('/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    linux=1
    mac=0
    print("\tWe're in Linux, looking for downloaded InSAR data accordingly.")
    file='/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
else:
    print("\tUnable to find downloaded InSAR data. Check for bigdata/InSAR folder. Aborting.")
    sys.exit()


InSAR_data = import_InSAR_csv(file)


Data = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
WellHlat = 36.32750
WellHlon = -119.58056

datesinsarH,InSAR_H = extract_series_from_latlon(WellHlat,WellHlon,InSAR_data)



if extraloc:
    borsa_lat = 36.34
    borsa_lon = -119.47
    datesinsarBULL,InSAR_BULL = extract_series_from_latlon(borsa_lat,borsa_lon,InSAR_data)


Swanson_1998_quote_dates = date2num([date(1980,1,1),date(1993,1,1)])
Swanson_1998_quote_data = 30.48*np.array([0,-2])

Poland_75_dates = date2num([date(1954,1,1),date(1957,1,1),date(1958,7,1),date(1962,1,1),date(1966,4,1),date(1970,1,1)])
Poland_75_data = -30.48*  np.array([0,0.83,1.34,3.25,4.92,6.31])

Highway_198_dates = date2num([date(1972,6,1),date(2004,1,1)])
Highway_198_data = 100*np.array([0,-0.84])
Highway_198_data_uncertainty = np.array([0,Poland_75_data[4]-Poland_75_data[2]])

Jacobus_data = pd.read_excel('/Users/mlees/Documents/RESEARCH/ground_deformation/GPS/From_Matt_Jacobus/Jacobus_Data_MASTER.xlsx')
Jacobus_data['YRMO']= pd.to_datetime(Jacobus_data['YRMO'],format='%Y%b')
S224P2data = Jacobus_data[Jacobus_data['STATION']=='S224P2']

# Import Envisat data
Envisat = import_InSAR_csv('/Users/mlees/Dropbox/Mer, RK, ML, RS NASA/InSAR/Datasets/Envisat/Envisat.csv')
Envisat_dates,Envisat_data = extract_series_from_latlon(36.32750,-119.58056,Envisat)


save = True
#%%

sns.set_style('whitegrid')
sns.set_context('talk')

fig,ax1 = plt.subplots(figsize=(18,12))

modelled_data_rezeroed = np.array(100 * rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980'))

ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5,label='Modelled (spin up?)')

ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'b--',label='Modelled (believable)')


ax1.plot_date(Poland_75_dates,Poland_75_data + modelled_data_rezeroed[Data['dates']=='1954-01-01'],'k.--',label='Poland 1975 levelling surveys')


ax1.plot_date(Highway_198_dates,Highway_198_data + modelled_data_rezeroed[Data['dates']=='1972-06-01'],'k^')
ax1.errorbar(Highway_198_dates[1],Highway_198_data[1] +modelled_data_rezeroed[Data['dates']=='1972-06-01'],xerr=None,yerr=0.5*Highway_198_data_uncertainty[1],fmt='k^',label='Highway 198 data',capsize=5)


ax1.arrow(Swanson_1998_quote_dates[1],30,-(Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0]),0,shape='full',head_width=15,head_length=300,facecolor='black',color='black',width=3)
ax1.arrow(Swanson_1998_quote_dates[0],30,Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0],0,shape='full',head_width=15,head_length=300,facecolor='black',color='black',width=3,label='Swanson (1998) period of subsidence')
txt = ax1.text(Swanson_1998_quote_dates[0] + (Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0])/2,40,'Swanson et al. (1998) period of subsidence',horizontalalignment='center',size='x-small')


ax1.plot_date(datesinsarH,0.1*InSAR_H + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]],'k-',label='Sentinel InSAR subsidence')



ax1.plot_date(Envisat_dates,Envisat_data + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==Envisat_dates[0])[0][0]],'r-',label='Envisat InSAR subsidence')


ax1.plot_date(S224P2data['YRMO'],foot_to_cm* ( S224P2data['ELEf'] - S224P2data['ELEf'][S224P2data['YRMO']==dt(2016,2,1)].values) +  modelled_data_rezeroed[Data['dates']==dt(2016,2,1)],label='S224P2 data')


plt.ylabel('Subsidence (cm)')

#plt.xlim([date.toordinal(date(2010,1,1)),date.toordinal(date(2020,1,1))])
#plt.ylim([-80,10])
plt.title('%s' % directory.split('/')[-1])

plt.legend()

if save:
    plt.savefig('figures/compare_with_measurementsFULL.png',bbox_inches='tight')
    plt.savefig('figures/compare_with_measurementsFULL.pdf',bbox_inches='tight')
    plt.savefig('figures/compare_with_measurementsFULL.svg',bbox_inches='tight')


yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]

txt.remove()
plt.xlim([date(2015,1,1),date(2020,1,1)])
yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]
plt.ylim([modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]-0.25*yrange,modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]]+0.25*yrange])


if extraloc:
    ax1.plot_date(datesinsarBULL,0.1*InSAR_BULL + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarBULL[0])[0][0]],'g-',label='Sentinel InSAR subsidence at Borsa Bullseye')

plt.legend()

plt.title('%s' % directory.split('/')[-1])

if save:
    plt.savefig('figures/compare_with_TRE_Altamira.png',bbox_inches='tight')
    plt.savefig('figures/compare_with_TRE_Altamira.svg',bbox_inches='tight')

os.chdir(cwd)