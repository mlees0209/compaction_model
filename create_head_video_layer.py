#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 12:47:42 2020

After the fact: creates a head video for a given layer.

@author: mlees
"""


import sys

if len(sys.argv) <= 1:
    print('create_head_video_layer.py error; terminal. Incorrect number of input arguments. Correct usage: python create_head_video_layer.py layername dt=XXXX. Layername should be fmt LAYER_THICKclay; dt=DAYS')
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

print('Reading head.')
if os.path.isfile("head_outputs/%s_head_data.nc" % layername):
    print('Head found as .nc file. Reading.')
    Dat = Dataset("head_outputs/%s_head_data.nc" % layername, "r", format="CF-1.7")
    head = Dat.variables['z'][:]
    z = Dat.variables['y'][:]
elif os.path.isfile("head_outputs/%s_head_data.csv" % layername):
    print('Head found as .csv file. Reading.')
    head=np.genfromtxt('head_outputs/%s_head_data.csv' % layername,delimiter=',')     
    zmax = float(layername.split('clay')[0].split('_')[-1])
    z = np.linspace(start=0,stop=zmax,num=np.shape(head)[0])
else:
    print('\tUnable to find head file as .nc or .csv. Something has gone wrong; aborting.')
    sys.exit(1)   

print('Reading inelastic flag...')
if os.path.isfile("head_outputs/%sinelastic_flag_GWFLOW.nb" % layername):
    print('\tFound as .nb file. Reading in...')
    InelasticFlagDat = Dataset("head_outputs/%sinelastic_flag_GWFLOW.nb" % layername, "r", format="CF-1.7")     
    if 'z' in list(InelasticFlagDat.variables.keys()):
        inelasticflag = InelasticFlagDat.variables['z'][:]
    elif 'W.nb' in list(InelasticFlagDat.variables.keys()):
        inelasticflag = InelasticFlagDat.variables['W.nb'][:]
    else:
        print('ERROR: cannot find variable name for Inelastic Flag. They are:. Quitting...')
        print(InelasticFlagDat.variables.keys())
        sys.exit(1)

elif os.path.isfile("head_outputs/%sinelastic_flag_GWFLOW.csv" % layername):
    print('\tFound as .csv file. Reading in...')
    inelasticflag = np.genfromtxt('head_outputs/%sinelastic_flag_GWFLOW.csv' % layername,delimiter=',')
else:
    print('\tUnable to find inelastic flag file as .nb or .csv. Something has gone wrong; aborting.')
    sys.exit(1)

print('Reading dates.')
dates_str_tmp = np.core.defchararray.rstrip(np.genfromtxt('head_outputs/%s_groundwater_solution_dates.csv' % layername,dtype=str,delimiter=','))
import datetime
print('Converting date formats.')
if dates_str_tmp[0][-1]=='M':
    dates_str = [datetime.datetime.strptime(string, '%a %d %b %Y %X %p').strftime('%d-%b-%Y') for string in dates_str_tmp]
else:
    dates_str = [datetime.datetime.strptime(string, '%c').strftime('%d-%b-%Y') for string in dates_str_tmp]


#%%

print('\t\tCreating video for %s' % layername)
create_head_video_elasticinelastic(head,z,inelasticflag==1,dates_str,'figures',layername,delt=dt,startyear=startyear,endyear=endyear,datelabels=datelabels)
