#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 23:49:34 2020

Reads in some InSAR and model data and makes some plots..

@author: mlees
"""

import sys

if len(sys.argv) != 2:
    print('compare_with_deformation_data error; terminal. Incorrect number of input arguments. Correct usage: python compare_with_deformation_data.py model_directory')
    sys.exit(1)

directory=sys.argv[1]

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

Swanson_1998_quote_dates = date2num([date(1980,1,1),date(1993,1,1)])
Swanson_1998_quote_data = 30.48*np.array([0,-2])

Poland_75_dates = date2num([date(1963,1,1),date(1966,4,1),date(1970,1,1)])
Poland_75_data = np.array([0,-3.3*0.65*30.48,-3.3*0.65*30.48 -3.6*0.3333*30.48])

Highway_198_dates = date2num([date(1965,1,1),date(2004,1,1)])
Highway_198_data = 100*np.array([0,-2.5])
Highway_198_data_uncertainty = np.array([0,Poland_75_data[2]])

save = True
#%%

sns.set_style('whitegrid')
sns.set_context('talk')

fig,ax1 = plt.subplots(figsize=(18,12))

modelled_data_rezeroed = np.array(100 * rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980'))

ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5,label='Modelled (spin up?)')

ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'b--',label='Modelled (believable)')


ax1.plot_date(Poland_75_dates,Poland_75_data + modelled_data_rezeroed[Data['dates']=='1963-01-01'],'k.--',label='Poland 1975 levelling surveys')


ax1.plot_date(Highway_198_dates,Highway_198_data + modelled_data_rezeroed[Data['dates']=='1965-01-01'],'k^')
ax1.errorbar(Highway_198_dates[1],Highway_198_data[1] +modelled_data_rezeroed[Data['dates']=='1965-01-01'],xerr=None,yerr=0.5*Highway_198_data_uncertainty[1],fmt='k^',label='Highway 198 data',capsize=5)


ax1.arrow(Swanson_1998_quote_dates[1],30,-(Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0]),0,shape='full',head_width=15,head_length=300,facecolor='black',color='black',width=3)
ax1.arrow(Swanson_1998_quote_dates[0],30,Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0],0,shape='full',head_width=15,head_length=300,facecolor='black',color='black',width=3,label='Swanson (1998) period of subsidence')
txt = ax1.text(Swanson_1998_quote_dates[0] + (Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0])/2,40,'Swanson et al. (1998) period of subsidence',horizontalalignment='center',size='x-small')


ax1.plot_date(datesinsarH,0.1*InSAR_H + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]],'k-',label='Sentinel InSAR subsidence')



plt.ylabel('Subsidence (cm)')

#plt.xlim([date.toordinal(date(2010,1,1)),date.toordinal(date(2020,1,1))])
#plt.ylim([-80,10])

plt.legend()

if save:
    plt.savefig('figures/compare_with_measurementsFULL.png',bbox_inches='tight')

yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]

txt.remove()
plt.xlim([date(2015,1,1),date(2020,1,1)])
yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]
plt.ylim([modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]-0.25*yrange,modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]]+0.25*yrange])


if save:
    plt.savefig('figures/compare_with_TRE_Altamira.png',bbox_inches='tight')
    plt.savefig('figures/compare_with_TRE_Altamira.svg',bbox_inches='tight')

os.chdir(cwd)