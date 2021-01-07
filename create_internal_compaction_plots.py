#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 13:43:15 2020

Create internal compaction plots for an entire simulation.

@author: mlees
"""

import sys

# if len(sys.argv) != 2:
#     print('compare_with_deformation_data error; terminal. Incorrect number of input arguments. Correct usage: python compare_with_deformation_data.py model_directory')
#     sys.exit(1)

# directory=sys.argv[1]
directory='/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/Model_Runs/Head_Timeseries_Family/Run32.5_Oct26head'

sys.path.append('/home/mlees/Land_Subsidence/Local_Scale/compaction_model')
sys.path.append('/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL')
from model_functions import *

print('Finding info on layers to plot.')

param_filename = 'paramfile.par'
f = open("%s/%s" % (directory,param_filename), "r")
paramfilelines = f.readlines()
f.close()
paramfilelines = [line.strip() for line in paramfilelines]
paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
paramfilelines[:] = [x for x in paramfilelines if x]

no_layers = read_parameter('no_layers',int,1,paramfilelines)
layer_types=read_parameter('layer_types',str,no_layers,paramfilelines)
layer_thicknesses=read_parameter('layer_thicknesses',float,no_layers,paramfilelines)

# Read in all the interbeds
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

# Add any aquitards
for layer in list(layer_types.keys()):
    if layer_types[layer]=='Aquitard':
        layername_dict[layer]=layer_thicknesses[layer]

print('Found %i clay layers to plot.' % len(layername_dict.keys()))


#%%
print("Reading dates (this always takes longer than you'd ever expect...).")
dates_dict = {}
for layer in list(layer_types.keys()):
    print('\t%s' % layer)
    layer = layer.replace(' ','_') # for some reason I saved outputs with an underscore.
    dates_str_tmp = np.core.defchararray.rstrip(np.genfromtxt('%s/head_outputs/%s_groundwater_solution_dates.csv' % (directory,layer),dtype=str,delimiter=','))
    import datetime
    print('\tConverting date formats.')
    if layer_types[layer.replace('_',' ')]=='Aquifer':
        dates_dict[layer.replace('_',' ')] = [datetime.datetime.strptime(string, '%c').strftime('%X %d-%b-%Y') for string in dates_str_tmp]
#    dates_dict[layer] = [datetime.datetime.strptime(string, '%c').strftime('%d-%b-%Y') for string in dates_str_tmp]


#%%

clay_layer = list(layername_dict.keys())[0]

plt.figure(figsize=(18,12))
t = dates_dict[clay_layer.split('_')[0]]


x_lims = [date2num(datetime.datetime.strptime(min(t),'%X %d-%b-%Y')),date2num(datetime.datetime.strptime(max(t),'%X %d-%b-%Y'))]

#%%
# Import the data

from netCDF4 import Dataset
Dat_tmp = Dataset("%s/s_outputs/%s_db.nc" % (directory,clay_layer.replace(' ','_')), "r", format="CF-1.7")
db_tmp = np.array(Dat_tmp.variables['z'][:])
Z_tmp = np.array(Dat_tmp.variables['y'][:])
#%%
import matplotlib.dates as mdates

y_lims=[0,max(Z_tmp)]


sns.set_style('white')
sns.set_context('talk')
plt.figure(figsize=(18,12))
plt.imshow(np.array(db_tmp),aspect='auto',cmap='RdBu',vmin=-np.max(np.abs(np.array(db_tmp)[5:,:])),vmax=np.max(np.abs(db_tmp)[5:,:]),extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]]) # note the min/max are set starting at the 5th timestep because the early timesteps can have large changes due to the initial condition and the boundary condition being discontinuous at these times
plt.gca().xaxis_date()
date_format = mdates.DateFormatter('%Y')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gcf().autofmt_xdate()
plt.colorbar(label='db (m)')
plt.ylabel('Z (m)')
#plt.savefig('%s/figures/%s/%sclay_compaction_internal.png' % (outdestination, layer,thickness),bbox_inches='tight')
#plt.close()
#plt.figure(figsize=(18,12))
#plt.imshow(np.array(db[layer]['total_%.2f clays' % thickness]).T,aspect='auto',cmap='RdBu',norm=colors.TwoSlopeNorm(vmin=-np.max(np.abs(np.array(db[layer]['total_%.2f clays' % thickness])[5:,:])), vcenter=-0.1*np.max(np.abs(np.array(db[layer]['total_%.2f clays' % thickness])[5:,:])), vmax=0)
#,extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]]) # note the min/max are set starting at the 5th timestep because the early timesteps can have large changes due to the initial condition and the boundary condition being discontinuous at these times
#plt.gca().xaxis_date()
#date_format = mdates.DateFormatter('%Y')
#plt.gca().xaxis.set_major_formatter(date_format)
#plt.gcf().autofmt_xdate()
#plt.colorbar(label='db (m)')
#plt.ylabel('Z (m)')
#plt.savefig('%s/figures/%s/%sclay_compaction_internalHIGCONTRAST.png' % (outdestination, layer,thickness),bbox_inches='tight')
#plt.close()




#%%

aquifer=re.split('_(\d+)',layername)[0]

Dat = Dataset("head_outputs/%s_head_data.nc" % layername, "r", format="CF-1.7")
head = Dat.variables['z'][:]
z = Dat.variables['y'][:]

InelasticFlagDat = Dataset("%sinelastic_flag_GWFLOW.nb" % layername, "r", format="CF-1.7")
inelasticflag = InelasticFlagDat.variables['z'][:]

print('Reading dates.')
dates_str_tmp = np.core.defchararray.rstrip(np.genfromtxt('head_outputs/%s_groundwater_solution_dates.csv' % aquifer,dtype=str,delimiter=','))
import datetime
print('Converting date formats.')
dates_str = [datetime.datetime.strptime(string, '%c').strftime('%d-%b-%Y') for string in dates_str_tmp]


