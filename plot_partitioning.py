#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 17:22:40 2020

@author: mlees
"""

import sys

if len(sys.argv) != 2:
    print('plot_partitioning error; terminal. Incorrect number of input arguments. Correct usage: python plot_partitioning.py runname')
    sys.exit(1)

directory=sys.argv[1]

import os
cwd = os.getcwd()
os.chdir(directory)

save=True

import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
sys.path.append('/home/mlees/InSAR_postprocessing/')
import seaborn as sns
from InSAR_postSBAS import *
from model_functions import *


print('Reading layer names...')
param_filename = 'paramfile.par'
f = open("%s" % (param_filename), "r")
paramfilelines = f.readlines()
f.close()
paramfilelines = [line.strip() for line in paramfilelines]
paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
paramfilelines[:] = [x for x in paramfilelines if x]
no_layers = read_parameter('no_layers',int,1,paramfilelines)
layer_names=read_parameter('layer_names',str,no_layers,paramfilelines)

Data = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
years = np.unique([a.strftime('%Y') for a in Data['dates']])

output={}
output['year']=years

results={}
for layer in layer_names:
    results[layer] = [float(Data[layer][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data[layer][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]


# upper = [float(Data['Upper Aquifer'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Upper Aquifer'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
# lower = [float(Data['Lower Aquifer'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Lower Aquifer'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
# cclay = [float(Data['Corcoran Clay'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Corcoran Clay'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
results['tot'] = [float(Data['Total'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Total'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]

pcs={}
for layer in layer_names:
    pcs[layer] = [100*(results[layer][i] / results['tot'][i]) for i in range(len(results[layer]))]

#%%
fig,ax1 = plt.subplots(figsize=(18,12))


#ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5)

ax1.plot_date([date2num(date) for date in Data['dates']],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980'),'k-',label='Modelled Subsidence')
plt.ylabel('Deformation (cm)')

ax2 = plt.twinx()

no_layers = len(layer_names)

i_tmp=0
for layer in layer_names:
    ax2.bar(365 + date2num(years[:-2])+(60 + (i_tmp -1) * 365/( no_layers*1.07)),pcs[layer],width=365/(no_layers*1.07),label=layer)
    i_tmp+=1

# ax2.bar(365 + date2num(years[:-2])+60,pc_upper,width=365/3.2,label='upper')
# ax2.bar(365 + date2num(years[:-2])+(60 + 365/3.2),pc_cclay,width=365/3.2,label='cclay')

plt.ylabel('% contribution')
plt.title('%s' % directory.split('/')[-1])

import matplotlib.dates as mdates
years_label = mdates.YearLocator()   # every year
#years_fmt = mdates.DateFormatter('%Y')
ax1.xaxis.set_minor_locator(years_label)


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax2.legend(lines + lines2, labels + labels2, fancybox=True)

if save:
    plt.savefig('figures/aquifer_partitioning_plot_PERCENT.png',bbox_inches='tight')
    plt.savefig('figures/aquifer_partitioning_plot_PERCENT.pdf',bbox_inches='tight')

#%%
# fig,ax1 = plt.subplots(figsize=(18,12))


# ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5)

# ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'k-',label='Modelled Subsidence')
# plt.ylabel('Deformation (cm)')

# ax2 = plt.twinx()
# ax2.bar(365 + date2num(years[:-2])+(60 - 365/3.2),100*np.array(lower),width=365/3.2,label='lower')
# ax2.bar(365 + date2num(years[:-2])+60,100*np.array(upper),width=365/3.2,label='upper')
# ax2.bar(365 + date2num(years[:-2])+(60 + 365/3.2),100*np.array(cclay),width=365/3.2,label='cclay')

# plt.ylabel('layer thickness change (cm)')
# plt.title('%s' % directory.split('/')[-1])

# import matplotlib.dates as mdates
# years_label = mdates.YearLocator()   # every year
# #years_fmt = mdates.DateFormatter('%Y')
# ax1.xaxis.set_minor_locator(years_label)


# lines, labels = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax2.legend(lines + lines2, labels + labels2, fancybox=True)

# if save:
#     plt.savefig('figures/aquifer_partitioning_plot_ABSOLUTE.png',bbox_inches='tight')
#     plt.savefig('figures/aquifer_partitioning_plot_ABSOLUTE.pdf',bbox_inches='tight')

os.chdir(cwd)


