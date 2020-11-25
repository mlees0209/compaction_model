#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:54:40 2020

Plot the contribution of different clay layer thicknesses to subsidence across the model.

@author: mlees
"""

import sys

if len(sys.argv) != 2:
    print('plot_partitioning_interbeds.py error; terminal. Incorrect number of input arguments. Correct usage: python plot_partitioning_interbeds.py model_directory')
    sys.exit(1)

directory=sys.argv[1]
#directory = '/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/Model_Runs/Ohm_Cutoffs_Family/Run34'
args=sys.argv[1:]

save=True

import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
sys.path.append('/home/mlees/InSAR_postprocessing/')

import seaborn as sns
import os
from InSAR_postSBAS import *
sys.path.append('/home/mlees/Land_Subsidence/Local_Scale/compaction_model')
sys.path.append('/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL')
from model_functions import *

cwd = os.getcwd()
os.chdir(directory)

print('Reading in data.')
Data = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
Data_lower = pd.read_csv('Lower Aquifer_Total_Deformation_Out.csv',parse_dates=[0])
Data_upper = pd.read_csv('Upper Aquifer_Total_Deformation_Out.csv',parse_dates=[0])

dat_all = {}

for col in list(Data_lower.columns):
    if col.startswith('total_'):
        dat_all['Lower Aquifer_%.2fclay' % float(col.split('_')[1].split(' ')[0])] = Data_lower[col].values

for col in list(Data_upper.columns):
    if col.startswith('total_'):
        dat_all['Upper Aquifer_%.2fclay' % float(col.split('_')[1].split(' ')[0])] = Data_upper[col].values


#%%
print('Finding interbed thicknesses')

param_filename = 'paramfile.par'
f = open("%s" % param_filename, "r")
paramfilelines = f.readlines()
f.close()
paramfilelines = [line.strip() for line in paramfilelines]
paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
paramfilelines[:] = [x for x in paramfilelines if x]

interbeds_distributions1=read_parameter('interbeds_distributions',dict,2,paramfilelines)
interbeds_distributions1=np.array(interbeds_distributions1)
if np.shape(interbeds_distributions1)[0]==1:
    interbeds_distributions1=interbeds_distributions1[0]
    minidics = [dict([(float(re.split(',|:',interbeds_distributions1[2*i + 1])[2*j]),int( re.split(',|:',interbeds_distributions1[2*i + 1])[2*j+1])) for j in range(int(len( re.split(',|:',interbeds_distributions1[2*i + 1]))/2))]) for i in range(sum(interbeds_switch.values()))]
    interbeds_distributions = dict([(interbeds_distributions1[2*i],minidics[i]) for i in range(sum(interbeds_switch.values()))])
    print('\tinterbeds_distributions=%s' % interbeds_distributions)
else:
    interbeds_distributions = {}
    for abc in interbeds_distributions1:
        interbeds_distributions[abc[0]] = dict([(float(re.split(':|,',abc[1])[2*i]),float(re.split(':|,',abc[1])[2*i+1])) for i in range(len(re.split(',',abc[1])))])
    print('\tinterbeds_distributions=%s' % interbeds_distributions)

layername_dict={}

for key in list(interbeds_distributions.keys()):
    for key2 in list(interbeds_distributions[key].keys()):
        layername_dict['%s_%.2fclay' % (key,key2)] = key2
#%%
lessthan5=[]
fiveto10=[]
greater10=[]

for key, value in layername_dict.items():
    if value < 5:
        lessthan5.append(key)
    if (value >= 5) & (value < 10):
        fiveto10.append(key)
    if value >=10:
        greater10.append(key)
        
#%% Now get the summed deformation from less than 5 etc
lessthan5_defm=np.zeros_like(Data['Total'])
lessthan5_thickness=0

for layer in lessthan5:
    lessthan5_defm += dat_all[layer]
    lessthan5_thickness+=layername_dict[layer] * interbeds_distributions[layer.split('_')[0]][float(layer.split('clay')[0].split('_')[1])]
    
fiveto10_defm=np.zeros_like(Data['Total'])
fiveto10_thickness=0

for layer in fiveto10:
    fiveto10_defm += dat_all[layer]
    fiveto10_thickness+=layername_dict[layer] * interbeds_distributions[layer.split('_')[0]][float(layer.split('clay')[0].split('_')[1])]

greater10_defm=np.zeros_like(Data['Total'])
greater10_thickness=0
for layer in greater10:
    greater10_defm += dat_all[layer]
    greater10_thickness+=layername_dict[layer] * interbeds_distributions[layer.split('_')[0]][float(layer.split('clay')[0].split('_')[1])]

print('Plotting deformation series for different clay thicknesses.')
plt.figure()
plt.plot_date(Data['dates'],Data['Total'],label='total def')
plt.plot_date(Data['dates'],lessthan5_defm,label='Interbeds thinner than 5 m; total thickness = %.2f' % lessthan5_thickness)
plt.plot_date(Data['dates'],fiveto10_defm,label='5-10 m interbeds; total thickness = %.2f' % fiveto10_thickness)
plt.plot_date(Data['dates'],greater10_defm,label='Interbeds thicker than 10 m; total thickness = %.2f' % greater10_thickness)

plt.legend()
if save:
    plt.savefig('figures/defm_by_interbeds.png',bbox_inches='tight')
    plt.savefig('figures/defm_by_interbeds.pdf',bbox_inches='tight')
    plt.savefig('figures/defm_by_interbeds.svg',bbox_inches='tight')


years = np.unique([a.strftime('%Y') for a in Data['dates']])

less5_ann = [float(lessthan5_defm[Data['dates']=='%s-09-30' % (int(year)+1)]) - float(lessthan5_defm[Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
fiveto10_ann = [float(fiveto10_defm[Data['dates']=='%s-09-30' % (int(year)+1)]) - float(fiveto10_defm[Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
greater10_ann = [float(greater10_defm[Data['dates']=='%s-09-30' % (int(year)+1)]) - float(greater10_defm[Data['dates']=='%s-10-01' % year]) for year in years[:-2]]

tot = [float(Data['Total'][Data['dates']=='%s-09-30' % (int(year)+1)]) - float(Data['Total'][Data['dates']=='%s-10-01' % year]) for year in years[:-2]]
#tot = np.array(less5_ann)+np.array(fiveto10_ann)+np.array(greater10_ann)

pc_less5 = [100*(less5_ann[i] / tot[i]) for i in range(len(less5_ann))]
pc_5to10 = [100*(fiveto10_ann[i] / tot[i]) for i in range(len(fiveto10_ann))]
pc_greater10 = [100*(greater10_ann[i] / tot[i]) for i in range(len(greater10_ann))]

#%%

print('Plotting % contribution of different interbeds.')
fig,ax1 = plt.subplots(figsize=(18,12))



ax2 = plt.twinx()
ax2.bar(365 + date2num(years[:-2])+(60 - 365/3.2),pc_less5,width=365/3.2,label='Interbeds thinner than 5 m; total thickness = %.2f' % lessthan5_thickness,color='lightblue')
ax2.bar(365 + date2num(years[:-2])+60,pc_5to10,width=365/3.2,label='5-10 m interbeds; total thickness = %.2f' % fiveto10_thickness,color='blue')
ax2.bar(365 + date2num(years[:-2])+(60 + 365/3.2),pc_greater10,width=365/3.2,label='Interbeds thicker than 10 m; total thickness = %.2f' % greater10_thickness,color='darkblue')
# ax2.plot_date(date2num(years[:-2]),pc_less5,'--',label='pc_less5')
# ax2.plot_date(date2num(years[:-2])+365/3,pc_5to10,'--',label='pc_5to10')
# ax2.plot_date(date2num(years[:-2])+2*365/3,pc_greater10,'--',label='pc_greater10')
plt.ylabel('% contribution')
plt.title('%s' % directory.split('/')[-1])

ax1.plot_date([date2num(date) for date in Data['dates']][0:365*20],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[0:365*20],'-',color='grey',linewidth=0.5)
ax1.plot_date([date2num(date) for date in Data['dates']][365*20:],100*rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Jun-1980')[365*20:],'k-',label='Modelled Subsidence')
plt.ylabel('Deformation (cm)')
ax1.set_zorder(2)
ax1.set_facecolor("none")

import matplotlib.dates as mdates
years = mdates.YearLocator()   # every year
ax1.xaxis.set_minor_locator(years)


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax2.legend(lines + lines2, labels + labels2, fancybox=True)

if save:
    plt.savefig('figures/partitioning_interbeds.png',bbox_inches='tight')
    plt.savefig('figures/partitioning_interbeds.pdf',bbox_inches='tight')
    plt.savefig('figures/partitioning_interbeds.svg',bbox_inches='tight')


os.chdir(cwd)

