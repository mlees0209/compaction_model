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

deephead = [var for var in sys.argv[1:] if var.split('=')[0]=='deephead']
if not deephead:
    deephead=None
else:
    deephead=bool(ut.strtobool(deephead[0].split('=')[1]))
    print('\thead specified; also plotting the deephead.')



foot_to_cm = 30.48

import os
cwd = os.getcwd()
os.chdir(directory)


import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
sys.path.append('/home/mlees/InSAR_postprocessing/')
import seaborn as sns
from InSAR_postSBAS import *

if ospath.exists('/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    mac=1
    linux=0
    knightblade=0
    print("\tWe're in Mac, looking for downloaded InSAR data accordingly.")
    file='/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
elif ospath.exists('/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    linux=1
    mac=0
    knightblade=0
    print("\tWe're in Linux, looking for downloaded InSAR data accordingly.")
    file='/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
elif ospath.exists('/data1/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    print("\tWe're in Knightblade, looking for downloading InSAR data accordingly.")
    mac=0
    linux=0
    knightblade=1
    file='/data1/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
else:
    print("\tUnable to find downloaded InSAR data. Check for bigdata/InSAR folder. Aborting.")
    sys.exit()


InSAR_data = import_InSAR_csv(file)

if deephead:
    print('Reading head.')
    H= pd.read_csv('input_data/input_time_series_Lower_Aquifer.csv',header=None)

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
Highway_198_data = 100*np.array([0,-1.7])
#Highway_198_data_uncertainty = np.array([0,Poland_75_data[4]-Poland_75_data[2]])

# Jacobus_data = pd.read_excel('/Users/mlees/Documents/RESEARCH/ground_deformation/GPS/From_Matt_Jacobus/Jacobus_Data_MASTER.xlsx')
# Jacobus_data['YRMO']= pd.to_datetime(Jacobus_data['YRMO'],format='%Y%b')
# S224P2data = Jacobus_data[Jacobus_data['STATION']=='S224P2']

# Import Envisat data
if linux:
    envisat_file ='/home/mlees/bigdata/InSAR/Processed_datasets/Tom_ENVISAT/Envisat.csv'
if mac:
    envisat_file ='/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/Tom_ENVISAT/Envisat.csv'
if knightblade:
    envisat_file='/data1/mlees/bigdata/InSAR/Processed_datasets/Tom_ENVISAT/Envisat.csv'


# Envisat = import_InSAR_csv(envisat_file)
# Envisat_dates,Envisat_data = extract_series_from_latlon(36.32750,-119.58056,Envisat)
# print('Finding envisat incidence angle.')
# envisat_incidence = Envisat['incidence_angle'].values[np.argmax([np.sum(Envisat.values[i,-24:] == Envisat_data) for i in range(len(Envisat))])]
# print('\tIncidence angle is %.2f deg.' % envisat_incidence)
# print()

# KSB0486_lat=36.2811
# KSB0486_lon = -119.6367
# Envisat_dates_KSB0486,Envisat_data_KSB0486 = extract_series_from_latlon(KSB0486_lat,KSB0486_lon,Envisat)
# print('Finding envisat incidence angle at KSB0486.')
# envisat_incidence_KSB0486 = Envisat['incidence_angle'].values[np.argmax([np.sum(Envisat.values[i,-24:] == Envisat_data_KSB0486) for i in range(len(Envisat))])]
# print('\tIncidence angle is %.2f deg.' % envisat_incidence_KSB0486)
# print()


# Import ALOS data
if linux:
    ALOS_file ='/home/mlees/bigdata/InSAR/Processed_datasets/ALOS/SmithKnight19_ALOS_SBAS.csv'
if mac:
    ALOS_file ='/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/ALOS/SmithKnight19_ALOS_SBAS.csv'
if knightblade:
    ALOS_file='/data1/mlees/bigdata/InSAR/Processed_datasets/Tom_ENVISAT/Envisat.csv'

ALOS = import_InSAR_csv(ALOS_file)
ALOS_dates,ALOS_data = extract_series_from_latlon(36.32750,-119.58056,ALOS)


save = True
#%%

if not deephead:
    sns.set_style('whitegrid')
sns.set_context('talk')

fig,ax1 = plt.subplots(figsize=(18,12))

modelled_data_rezeroed = np.array(100 * rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980'))

#modelled_data_adjusted_rezeroed = np.array(100 * rezero_series(Data_variablethickness['Total'],np.array([date2num(date) for date in Data_variablethickness['dates']]),'Jun-1980')

Sentinel_rezeroed = np.array(0.1 * rezero_series(InSAR_H,datesinsarH,'Jan-2016'))


ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'--',color='grey',linewidth=0.5,label='Modelled (spin up?)')

ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'b--',label='Modelled (believable)')

ax1.plot_date(Poland_75_dates,Poland_75_data + modelled_data_rezeroed[Data['dates']=='1966-04-01'] - Poland_75_data[Poland_75_dates==date2num(date(1966,4,1))],'k.--',label='Poland 1975 levelling surveys')


ax1.plot_date(Highway_198_dates,Highway_198_data + modelled_data_rezeroed[Data['dates']=='1972-06-01'],'k^',label='Highway 198 levelling surveys')
#ax1.errorbar(Highway_198_dates[1],Highway_198_data[1] +modelled_data_rezeroed[Data['dates']=='1972-06-01'],xerr=None,yerr=0.5*Highway_198_data_uncertainty[1],fmt='k^',label='Highway 198 data',capsize=5)


ax1.arrow(Swanson_1998_quote_dates[1],30,-(Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0]),0,shape='full',head_width=15,head_length=300,facecolor='black',color='black',width=3)
ax1.arrow(Swanson_1998_quote_dates[0],30,Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0],0,shape='full',head_width=15,head_length=300,facecolor='black',color='black',width=3,label='Swanson (1998) period of subsidence')
txt = ax1.text(Swanson_1998_quote_dates[0] + (Swanson_1998_quote_dates[1]-Swanson_1998_quote_dates[0])/2,40,'Swanson et al. (1998) period of subsidence',horizontalalignment='center',size='x-small')


ax1.plot_date(datesinsarH,Sentinel_rezeroed + modelled_data_rezeroed[np.argmin(np.abs(date2num(Data['dates'])- date2num(date(2016,1,1))))],'k-',label='Sentinel InSAR subsidence')


rezero_idx_env =1
ax1.plot_date(ALOS_dates,-(ALOS_data - ALOS_data[rezero_idx_env]) + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==ALOS_dates[rezero_idx_env])[0][0]],'p-',label='ALOS InSAR subsidence (vert)')
# ax1.plot_date(Envisat_dates_KSB0486,0.1/np.cos(np.deg2rad(envisat_incidence_KSB0486)) * (Envisat_data_KSB0486 - Envisat_data_KSB0486[rezero_idx_env]) + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==Envisat_dates_KSB0486[rezero_idx_env])[0][0]],'r-',label='Envisat InSAR subsidence (KSB-486 vert)')


# ax1.plot_date(S224P2data['YRMO'],foot_to_cm* ( S224P2data['ELEf'] - S224P2data['ELEf'][S224P2data['YRMO']==dt(2016,2,1)].values) +  modelled_data_rezeroed[Data['dates']==dt(2016,2,1)],label='S224P2 data')


plt.ylabel('Subsidence (cm)')

#plt.xlim([date.toordinal(date(2010,1,1)),date.toordinal(date(2020,1,1))])
#plt.ylim([-80,10])
plt.title('%s' % directory.split('/')[-1])

leg = plt.legend()

if deephead:
    ax2 = ax1.twinx()
    leg.remove()
    ax2.plot_date(H[0],H[1],'k-',linewidth=0.5,label='Deep aquifer head')
    plt.ylabel('Head (masl)')
    
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    
    #ax1.get_legend().remove()
    ax2.legend(lines + lines2, labels + labels2,fancybox=True)


if save:
    plt.savefig('figures/compare_with_measurementsFULL.png',bbox_inches='tight')
    plt.savefig('figures/compare_with_measurementsFULL.pdf',bbox_inches='tight')
    plt.savefig('figures/compare_with_measurementsFULL.svg',bbox_inches='tight')


yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]

if deephead:
    ax2.remove()

txt.remove()
plt.xlim([date(2015,1,1),date(2020,1,1)])
yrange = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]
plt.ylim([modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[-1])[0][0]]-0.25*yrange,modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarH[0])[0][0]]+0.25*yrange])


if extraloc:
    ax1.plot_date(datesinsarBULL,0.1*InSAR_BULL + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarBULL[0])[0][0]],'g-',label='Sentinel InSAR subsidence at Borsa Bullseye')

leg = plt.legend()

plt.title('%s' % directory.split('/')[-1])

if deephead:
    leg.remove()
    ax2 = ax1.twinx()
    
    ax2.plot_date(H[0],H[1],'k-',linewidth=0.5,label='Deep aquifer head')
    plt.ylabel('Head (masl)')
    plt.ylim([-20,40])
    
    ax2.legend(lines + lines2, labels + labels2,fancybox=True)


if save:
    if deephead:
        plt.savefig('figures/compare_with_TRE_Altamira_deephead.png',bbox_inches='tight')
        plt.savefig('figures/compare_with_TRE_Altamira_deephead.svg',bbox_inches='tight')
    else:
        plt.savefig('figures/compare_with_TRE_Altamira.png',bbox_inches='tight')
        plt.savefig('figures/compare_with_TRE_Altamira.svg',bbox_inches='tight')

#%%

# yrange_envisat = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==Envisat_dates[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==Envisat_dates[-1])[0][0]]

# if deephead:
#     ax2.remove()

# plt.xlim([date(2005,1,1),date(2011,1,1)])
# plt.ylim([modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==Envisat_dates[-1])[0][0]]-0.25*yrange_envisat,modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==Envisat_dates[0])[0][0]]+0.25*yrange_envisat])


# if extraloc:
#     ax1.plot_date(datesinsarBULL,0.1*InSAR_BULL + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarBULL[0])[0][0]],'g-',label='Sentinel InSAR subsidence at Borsa Bullseye')

# leg = plt.legend()

# plt.title('%s' % directory.split('/')[-1])

# if deephead:
#     leg.remove()
#     ax2 = ax1.twinx()
    
#     ax2.plot_date(H[0],H[1],'k-',linewidth=0.5,label='Deep aquifer head')
#     plt.ylabel('Head (masl)')
#     plt.ylim([-20,40])
    
#     ax2.legend(lines + lines2, labels + labels2,fancybox=True)


# if save:
#     if deephead:
#         plt.savefig('figures/compare_with_ENVISAT_deephead.png',bbox_inches='tight')
#         plt.savefig('figures/compare_with_ENVISAT_deephead.svg',bbox_inches='tight')
#     else:
#         plt.savefig('figures/compare_with_ENVISAT.png',bbox_inches='tight')
#         plt.savefig('figures/compare_with_ENVISAT.svg',bbox_inches='tight')

#%%

yrange_ALOS = modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==ALOS_dates[0])[0][0]] - modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==ALOS_dates[-1])[0][0]]

if deephead:
    ax2.remove()

plt.xlim([date(2005,1,1),date(2011,1,1)])
plt.ylim([modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==ALOS_dates[-1])[0][0]]-0.25*yrange_ALOS,modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==ALOS_dates[0])[0][0]]+0.25*yrange_ALOS])


if extraloc:
    ax1.plot_date(datesinsarBULL,0.1*InSAR_BULL + modelled_data_rezeroed[np.where([date2num(date) for date in Data['dates']]==datesinsarBULL[0])[0][0]],'g-',label='Sentinel InSAR subsidence at Borsa Bullseye')

leg = plt.legend()

plt.title('%s' % directory.split('/')[-1])

if deephead:
    leg.remove()
    ax2 = ax1.twinx()
    
    ax2.plot_date(H[0],H[1],'k-',linewidth=0.5,label='Deep aquifer head')
    plt.ylabel('Head (masl)')
    plt.ylim([-20,40])
    
    ax2.legend(lines + lines2, labels + labels2,fancybox=True)


if save:
    if deephead:
        plt.savefig('figures/compare_with_ALOS_deephead.png',bbox_inches='tight')
        plt.savefig('figures/compare_with_ALOS_deephead.svg',bbox_inches='tight')
    else:
        plt.savefig('figures/compare_with_ALOS.png',bbox_inches='tight')
        plt.savefig('figures/compare_with_ALOS.svg',bbox_inches='tight')




#%% Extra plot for "adjusted"


# print('Making version with scaled lower aquifer prior to 2015')
# Data_variablethickness = Data.copy()
# Data_variablethickness['Lower Aquifer'][date2num(Data_variablethickness['dates']) <= date2num(dt(2015,1,1))] = Data['Lower Aquifer'][date2num(Data_variablethickness['dates']) <= date2num(dt(2015,1,1))] * 0.5
# Data_variablethickness['Lower Aquifer'][date2num(Data_variablethickness['dates']) > date2num(dt(2015,1,1))] = Data_variablethickness['Lower Aquifer'][date2num(Data_variablethickness['dates']) == date2num(dt(2015,1,1))].values + (Data['Lower Aquifer'][date2num(Data_variablethickness['dates']) > date2num(dt(2015,1,1))] - Data['Lower Aquifer'][date2num(Data_variablethickness['dates']) == date2num(dt(2015,1,1))].values)
# Data_variablethickness['Total'] = Data_variablethickness['Upper Aquifer'] + Data_variablethickness['Corcoran Clay'] + Data_variablethickness['Lower Aquifer']


# modelled_data_rezeroed_adusted = np.array(100 * rezero_series(Data_variablethickness['Total'],np.array([date2num(date) for date in Data_variablethickness['dates']]),'Jun-1980'))
# Sentinel_rezeroed = np.array(0.1 * rezero_series(InSAR_H,datesinsarH,'Jan-2016'))


# fig,ax1 = plt.subplots(figsize=(18,12))


# # ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'b--',label='Modelled (believable)')

# ax1.plot_date([date2num(date) for date in Data_variablethickness['dates']][0:365*20],100*rezero_series(Data_variablethickness['Total'],np.array([date2num(date) for date in Data_variablethickness['dates']]),'Jun-1980')[0:365*20],'--',color='grey',linewidth=0.5,label='Modelled (spin up?)')

# ax1.plot_date([date2num(date) for date in Data_variablethickness['dates']][365*20:],100*rezero_series(Data_variablethickness['Total'],np.array([date2num(date) for date in Data_variablethickness['dates']]),'Jun-1980')[365*20:],'r--',label='Modelled adjusted (believable)')

# ax1.plot_date(Highway_198_dates,Highway_198_data + modelled_data_rezeroed_adusted[Data_variablethickness['dates']=='1972-06-01'],'k^',label='Highway 198 data')

# ax1.plot_date(Poland_75_dates,Poland_75_data + modelled_data_rezeroed_adusted[Data['dates']=='1966-04-01'] - Poland_75_data[Poland_75_dates==date2num(date(1966,4,1))],'k.--',label='Poland 1975 levelling surveys')

# rezero_idx_env =2
# ax1.plot_date(Envisat_dates,0.1/np.cos(np.deg2rad(envisat_incidence)) * (Envisat_data - Envisat_data[rezero_idx_env]) + modelled_data_rezeroed_adusted[np.where([date2num(date) for date in Data['dates']]==Envisat_dates[rezero_idx_env])[0][0]],'p-',label='Envisat InSAR subsidence (vert)')
# ax1.plot_date(Envisat_dates_KSB0486,0.1/np.cos(np.deg2rad(envisat_incidence_KSB0486)) * (Envisat_data_KSB0486 - Envisat_data_KSB0486[rezero_idx_env]) + modelled_data_rezeroed_adusted[np.where([date2num(date) for date in Data['dates']]==Envisat_dates_KSB0486[rezero_idx_env])[0][0]],'o-',label='Envisat InSAR subsidence (KSB-486 vert)')

# ax1.plot_date(datesinsarH,Sentinel_rezeroed + modelled_data_rezeroed_adusted[np.argmin(np.abs(date2num(Data['dates'])- date2num(date(2016,1,1))))],'k-',label='Sentinel InSAR subsidence')

# plt.legend()

# plt.title('%s (with expanding lower aquifer in 2015)' % directory.split('/')[-1])

# if save:
#     plt.savefig('figures/compare_with_measurementsFULL_adjusted.png',bbox_inches='tight')

os.chdir(cwd)
