#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 12:47:42 2020

After the fact: creates a head video for a given layer.

@author: mlees
"""


import sys

if len(sys.argv) != 2:
    print('create_head_video_layer.py error; terminal. Incorrect number of input arguments. Correct usage: python create_head_video_layer.py layername. Layername should be fmt LAYER_THICKclay')
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



#%%

print('\t\tCreating video for %s' % layername)
create_head_video_elasticinelastic(head,z,inelasticflag==1,dates_str,'figures',layername,delt=30)
