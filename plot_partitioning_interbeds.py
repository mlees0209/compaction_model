#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:54:40 2020

Plot the contribution of different clay layer thicknesses to subsidence across the model.

@author: mlees
"""

import sys

# if len(sys.argv) <= 2:
#     print('plot_subsidence_breakdown.py error; terminal. Incorrect number of input arguments. Correct usage: python breakdown_subsidence.py model_directory endyear=XXXX')
#     sys.exit(1)

#directory=sys.argv[1]
directory = '/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/Model_Runs/Ohm_Cutoffs_Family/Run34'
#args=sys.argv[1:]

# endyear = [var for var in args if var.split('=')[0]=='endyear']
# if not endyear:
#     print('No endyear specified; ending on Sep 1st 2019.')
#     endyear=2019
# else:
#     endyear=int(endyear[0].split('=')[1])
#     print('Endyear specified; zeroing on Sep 1st %i.' % endyear)
startyear=2015

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
for layer in lessthan5:
    lessthan5_defm += dat_all[layer]

fiveto10_defm=np.zeros_like(Data['Total'])
for layer in fiveto10:
    fiveto10_defm += dat_all[layer]

greater10_defm=np.zeros_like(Data['Total'])
for layer in greater10:
    greater10_defm += dat_all[layer]


plt.figure()
plt.plot_date(Data['dates'],Data['Total'],label='total def')
plt.plot_date(Data['dates'],lessthan5_defm,label='<5 m clay def')
plt.plot_date(Data['dates'],fiveto10_defm,label='5-10 m clay def')
plt.plot_date(Data['dates'],greater10_defm,label='>10 m clay def')

plt.legend()
plt.show()