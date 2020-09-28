#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 23:49:34 2020

Reads in some InSAR and model data and makes some plots..

@author: mlees
"""

import sys

if len(sys.argv) != 2:
    print('compare_with_deformation_data error; terminal. Incorrect number of input arguments. Correct usage: python execute_model.py parameter_file.par')
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

InSAR_data = import_InSAR_csv('/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv')


Data = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
WellHlat = 36.32750
WellHlon = -119.58056

datesinsarH,InSAR_H = extract_series_from_latlon(WellHlat,WellHlon,InSAR_data)

Highway_198_dates = [date.toordinal(date(1965,1,1)),date.toordinal(date(2004,1,1))]
Highway_198_data = 100*np.array([0,-2.5])

Swanson_1998_quote_dates = [date.toordinal(date(1980,1,1)),date.toordinal(date(1993,1,1))]
Swanson_1998_quote_data = 30.48*np.array([0,-2])

Poland_75_dates = [date.toordinal(date(1963,1,1)),date.toordinal(date(1966,4,1)),date.toordinal(date(1970,1,1))]
Poland_75_data = np.array([0,-3.3*0.65*30.48,-3.3*0.65*30.48 -3.6*0.3333*30.48])


save = True
#%%

sns.set_style('whitegrid')
sns.set_context('talk')

fig,ax1 = plt.subplots(figsize=(18,12))

modelled_data_rezeroed = np.array(100 * rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980'))

ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5,label='Modelled (spin up?)')

ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'b--',label='Modelled (believable)')


ax1.plot_date(Poland_75_dates,Poland_75_data + 180,'k.--',label='Poland 1975 levelling surveys')

ax1.plot_date(Highway_198_dates,Highway_198_data + 100,'k^',label='Highway 198 data')


ax1.plot_date(Swanson_1998_quote_dates,Swanson_1998_quote_data ,'b^',label='Swanson 1998 data (v approximate)')


ax1.plot_date(datesinsarH,0.1*InSAR_H + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]],'k-',label='Sentinel InSAR subsidence')



plt.ylabel('Subsidence (cm)')

#plt.xlim([date.toordinal(date(2010,1,1)),date.toordinal(date(2020,1,1))])
#plt.ylim([-80,10])

plt.legend()

if save:
    plt.savefig('compare_with_measurementsFULL.png',bbox_inches='tight')

plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - np.min(modelled_data_rezeroed)
plt.ylim([np.min(modelled_data_rezeroed)-0.25*yrange,modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]]+0.25*yrange])


if save:
    plt.savefig('compare_with_TRE_Altamira.png',bbox_inches='tight')
    plt.savefig('compare_with_TRE_Altamira.svg',bbox_inches='tight')

os.chdir(cwd)