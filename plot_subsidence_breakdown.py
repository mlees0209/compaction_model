#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:29:37 2020

Script creates two plots: one breaks down the 2015-2020 deformation by aquifer, and the other breaks it down by individual clay layers. Optional parameter is enddate, which can be as late as you jolly well like!

@author: mlees
"""

#%%

import sys

if len(sys.argv) <= 2:
    print('plot_subsidence_breakdown.py error; terminal. Incorrect number of input arguments. Correct usage: python breakdown_subsidence.py model_directory endyear=XXXX')
    sys.exit(1)

directory=sys.argv[1]

args=sys.argv[1:]

endyear = [var for var in args if var.split('=')[0]=='endyear']
if not endyear:
    print('No endyear specified; ending on Sep 1st 2019.')
    endyear=2019
else:
    endyear=int(endyear[0].split('=')[1])
    print('Endyear specified; zeroing on Sep 1st %i.' % endyear)
startyear=2015

save=True

import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
import seaborn as sns
import os
from InSAR_postSBAS import *

cwd = os.getcwd()
os.chdir(directory)

print('Reading in data.')
Data = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
Data_lower = pd.read_csv('Lower Aquifer_Total_Deformation_Out.csv',parse_dates=[0])
Data_upper = pd.read_csv('Upper Aquifer_Total_Deformation_Out.csv',parse_dates=[0])


#%%
print('Plotting aquifer breakdown.')
sns.set_context('poster')
TOT_rezeroed = np.array(100 * rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Sep-%s' % startyear))

UPPER_TOT_rezeroed = np.array(100 * rezero_series(Data_upper['total'],np.array([date2num(date) for date in Data_upper['dates']]),'Sep-%s' % startyear))
LOWER_TOT_rezeroed = np.array(100 * rezero_series(Data_lower['total'],np.array([date2num(date) for date in Data_lower['dates']]),'Sep-%s' % startyear))

plt.figure(figsize=(18,12))
plt.plot_date(Data['dates'],TOT_rezeroed,'-',label='Total')
plt.plot_date(Data['dates'],UPPER_TOT_rezeroed,'-',label='Upper aquifer only')
plt.plot_date(Data['dates'],LOWER_TOT_rezeroed,'-',label='Lower aquifer only')
plt.xlim([dt(startyear,1,1),dt(endyear,12,31)])
yrange = np.abs(TOT_rezeroed[np.where(Data['dates']==dt(startyear,1,1))] - TOT_rezeroed[np.where(Data['dates']==dt(endyear,12,31))])
print('\tyrange is ',yrange)
plt.ylim([TOT_rezeroed[np.where(Data['dates']==dt(endyear,12,31))]-0.25*yrange,0+0.25*yrange])
plt.ylabel('Subsidence (cm)')
plt.legend()
if save:
    plt.savefig('figures/subsidence_breakdown_aquifer.png',bbox_inches='tight')
    plt.savefig('figures/subsidence_breakdown_aquifer.pdf',bbox_inches='tight')
    plt.savefig('figures/subsidence_breakdown_aquifer.svg',bbox_inches='tight')

#%%
print('Plotting interbed breakdown')
layers_lower=[layer for layer in list(Data_lower.columns) if layer.startswith('total_')]
print('\tinterbeds in lower are ', layers_lower)
layers_upper=[layer for layer in list(Data_upper.columns) if layer.startswith('total_')]
print('\tinterbeds in upper are ', layers_upper)


#head_in_lower = pd.read_csv('input_head_data/input_time_series_Lower_Aquifer.csv')


fig,ax1 = plt.subplots(figsize=(18,12))
plt.plot_date(Data['dates'],TOT_rezeroed,'k-',label='Tot')

for layer in layers_lower:
    plt.plot_date(Data_lower['dates'],np.array(100 * rezero_series(Data_lower[layer],np.array([date2num(date) for date in Data_lower['dates']]),'Jul-2015')),'-.',label=layer)

for layer in layers_upper:
    plt.plot_date(Data_upper['dates'],np.array(100 * rezero_series(Data_upper[layer],np.array([date2num(date) for date in Data_upper['dates']]),'Jul-2015')),'--',label=layer)

plt.xlim([dt(startyear,1,1),dt(endyear,12,31)])

plt.ylim([TOT_rezeroed[np.where(Data['dates']==dt(endyear,12,31))]-0.25*yrange,0+0.25*yrange])
plt.ylabel('Subsidence (cm)')

# Sort the nightmare of a multi-column legend
h, l = ax1.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*2
handles = [h[0]] + ph[:1] + h[1:len(layers_lower)+1] + ph[1:] + h[1+len(layers_lower):]
labels = [l[0]] + ["Lower"] + l[1:len(layers_lower)+1] + ["Upper"] + l[1+len(layers_lower):]
leg = plt.legend(handles, labels, ncol=2)

# for vpack in leg._legend_handle_box.get_children():
#     for hpack in vpack.get_children()[:1]:
#         hpack.get_children()[0].set_width(0)


plt.title('broken down by interbed thickness')

# ax2 = ax1.twinx()
# ax2.plot_date(head_in_lower.values[:,0],head_in_lower.values[:,1])
if save:
    plt.savefig('figures/subsidence_breakdown_layers.png',bbox_inches='tight')
    plt.savefig('figures/subsidence_breakdown_layers.pdf',bbox_inches='tight')
    plt.savefig('figures/subsidence_breakdown_layers.svg',bbox_inches='tight')
