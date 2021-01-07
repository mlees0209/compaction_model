#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:59:47 2020

After the fact, create a compaction video for a given clay interbed.

@author: mlees
"""
#!/usr/bin/env python3


import sys

if len(sys.argv) <= 2:
    print('create_compaction_video_layer.py error; terminal. Incorrect number of input arguments. Correct usage: python create_compaction_video_layer.py layername dt=XXXX. Layername should be fmt LAYER_THICKclay; dt=DAYS')
    sys.exit(1)


import os


import pandas as pd
sys.path.append('/home/mlees/Land_Subsidence/Local_Scale/compaction_model')
sys.path.append('/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL')
from model_functions import *
import seaborn as sns
import re


from netCDF4 import Dataset
import numpy as np

layername=sys.argv[1]
args=sys.argv[1:]

dt = [var for var in args if var.split('=')[0]=='dt']
if not dt:
    print('No dt specified; dt=30 days.')
    dt=30
else:
    dt=int(dt[0].split('=')[1])
    print('dt specified; dt=%i days.' % dt)

startyear = [var for var in args if var.split('=')[0]=='startyear']
if not startyear:
    print('No startyear specified; starting on day 0.')
    startyear=None
else:
    startyear=int(startyear[0].split('=')[1])
    print('startyear specified as %s.' % startyear)

endyear = [var for var in args if var.split('=')[0]=='endyear']
if not endyear:
    print('No endyear specified; starting on day 0.')
    endyear=None
else:
    endyear=int(endyear[0].split('=')[1])
    print('endyear specified as %s.' % endyear)

datelabels = [var for var in args if var.split('=')[0]=='datelabels']
if not datelabels:
    datelabels='year'
else:
    datelabels=str(datelabels[0].split('=')[1])
    print('datelabels specified as %s.' % datelabels)


save=True

aquifer=re.split('_(\d+)',layername)[0]

print('Reading db,z,time.')
Dat = Dataset("s_outputs/%s_db.nc" % layername, "r", format="CF-1.7")
db = Dat.variables['z'][:]
z = Dat.variables['y'][:]
dz = np.diff(z)[0]
time_sim = Dat.variables['time'][:]
duration = (max(time_sim)-min(time_sim)) / 365
delta_t = np.diff(time_sim)[0]

print('Reading total aquifer thickness...')
param_filename = 'paramfile.par'
f = open("%s" % param_filename, "r")
paramfilelines = f.readlines()
f.close()
paramfilelines = [line.strip() for line in paramfilelines]
paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
paramfilelines[:] = [x for x in paramfilelines if x]
no_layers = read_parameter('no_layers',int,1,paramfilelines)
layer_names=read_parameter('layer_names',str,no_layers,paramfilelines)
if no_layers==1:
    layer_names = np.array([layer_names])
layer_types=read_parameter('layer_types',str,no_layers,paramfilelines)
layer_thicknesses=read_parameter('layer_thicknesses',float,no_layers,paramfilelines)
total_model_thickness = np.sum(list(layer_thicknesses.values()))

print('Forming scaling factor.')
scaling_factor = 100 * total_model_thickness/dz * (365/delta_t) # obtains units of cm per year per unit thickness aquifer
db_plot = db * scaling_factor


print('Reading inelastic flag...')
if os.path.isfile("s_outputs/%sinelastic_flag_COMPACTION.nb" % layername):
    print('\tFound as .nb file. Reading in...')
    InelasticFlagDat = Dataset("s_outputs/%sinelastic_flag_COMPACTION.nb" % layername, "r", format="CF-1.7")     
    Inelastic_Flag = InelasticFlagDat.variables['z'][:].astype(bool)
#elif os.path.isfile("head_outputs/%sinelastic_flag_GWFLOW.csv" % layername):
#    print('\tFound as .csv file. Reading in...')
#    inelasticflag = np.genfromtxt('head_outputs/%sinelastic_flag_GWFLOW.csv' % layername,delimiter=',')
else:
    print('\tUnable to find inelastic flag file as .nb and .csv not yet coded. Something has gone wrong; aborting.')
    sys.exit(1)

#%%

print('\t\tCreating video for %s' % layername)
create_compaction_video('figures',layername,db_plot,z,time_sim,Inelastic_Flag,delt=dt,startyear=startyear,endyear=endyear,datelabels=datelabels)