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
import seaborn as sns
from InSAR_postSBAS import *


Data = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
years = np.unique([a.strftime('%Y') for a in Data['dates']])

output={}
output['year']=years


upper = [float(Data['Upper Aquifer'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Upper Aquifer'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
lower = [float(Data['Lower Aquifer'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Lower Aquifer'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
cclay = [float(Data['Corcoran Clay'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Corcoran Clay'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
tot = np.array(upper)+np.array(lower)+np.array(cclay)
pc_upper = [100*(upper[i] / tot[i]) for i in range(len(upper))]
pc_lower = [100*(lower[i] / tot[i]) for i in range(len(upper))]
pc_cclay = [100*(cclay[i] / tot[i]) for i in range(len(upper))]

#%%
fig,ax1 = plt.subplots(figsize=(18,12))


ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5)

ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'k-',label='Modelled Subsidence')
plt.ylabel('Deformation (cm)')

ax2 = plt.twinx()
ax2.bar(date2num(years[:-2]),pc_lower,width=365/3,label='lower')
ax2.bar(date2num(years[:-2])+365/3,pc_upper,width=365/3,label='upper')
ax2.bar(date2num(years[:-2])+2*365/3,pc_cclay,width=365/3,label='cclay')

plt.ylabel('% contribution')



lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax2.legend(lines + lines2, labels + labels2, fancybox=True)

if save:
    plt.savefig('figures/aquifer_partitioning_plot.png',bbox_inches='tight')
    plt.savefig('figures/aquifer_partitioning_plot.png',bbox_inches='tight')

os.chdir(cwd)

