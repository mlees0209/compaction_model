#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:01:31 2020

partition between two dates

@author: mlees
"""

import sys

if len(sys.argv) != 4:
    print('partition_aquifer_datepair error; terminal. Incorrect number of input arguments. Correct usage: python plot_partitioning.py runname startdate enddate. date formats YYYY-MM-DD.')
    sys.exit(1)

directory=sys.argv[1]
startdate=sys.argv[2]
enddate = sys.argv[3]

import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
import seaborn as sns
import numpy as np

Data = pd.read_csv('%s/Total_Deformation_Out.csv' % directory,parse_dates=[0])
years = np.unique([a.strftime('%Y') for a in Data['dates']])

output={}
output['year']=years


upper = float(Data['Upper Aquifer'][Data['dates']==enddate]) - float(Data['Upper Aquifer'][Data['dates']==startdate])
lower = float(Data['Lower Aquifer'][Data['dates']==enddate]) - float(Data['Lower Aquifer'][Data['dates']==startdate])
cclay = float(Data['Corcoran Clay'][Data['dates']==enddate]) - float(Data['Corcoran Clay'][Data['dates']==startdate])
tot = float(Data['Total'][Data['dates']==enddate]) - float(Data['Total'][Data['dates']==startdate])

print('Deformation between %s and %s is %.2f m. Partitioned between:' % (startdate, enddate, tot))
print('\tUpper Aquifer: %.2f m (%.2f %%)' % (upper, 100*(upper/tot)))
print('\tLower Aquifer: %.2f m (%.2f %%)' % (lower, 100*(lower/tot)))
print('\tCorcoran Clay: %.2f m (%.2f %%)' % (cclay, 100*(cclay/tot)))