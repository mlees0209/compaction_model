#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 16:41:38 2020

Simple script to compare N runs. OPTIONAL INPUT ARGUMENT 'startyear=XXXX' allows the run to be rezerod on a year of your choice; default is Jun 1980 for some reason.

@author: mlees
"""
import sys

if len(sys.argv) <= 2:
    print('compare_with_deformation_data error; terminal. Incorrect number of input arguments. Correct usage: python compare_models.py model1 model2 model3 ... modelN.')
    sys.exit(1)


import os


import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
import seaborn as sns
from InSAR_postSBAS import *

no_runs = len(sys.argv)-1

run_names = sys.argv[1:]
startyear = [var for var in run_names if var.split('=')[0]=='startyear']
if not startyear:
    print('No startyear specified; zeroing on Sep 1st 1980.')
    startyear=1980
else:
    startyear=int(startyear[0].split('=')[1])
    print('Startyear specified; zeroing on Sep 1st %i.' % startyear)
run_names = [var for var in run_names if var.split('=')[0]!='startyear']
print("Plotting for runs ", run_names)

cwd = os.getcwd()
Data={}

sns.set_context('poster')
sns.set_style('darkgrid')

fig,ax1 = plt.subplots(figsize=(18,12))


for run in run_names:
    print('\tPlotting %s.' % run)
    os.chdir(run)
    Data[run] = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
    ax1.plot_date([date2num(date) for date in Data[run]['dates']][0:365*20],100*rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[0:365*20],'-',color='grey',linewidth=0.5,label='_nolegend_')

    ax1.plot_date([date2num(date) for date in Data[run]['dates']][365*20:],100*rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[365*20:],'-',label=run)
    os.chdir(cwd)

plt.legend()
plt.ylabel('Deformation (cm)')
s='COMPARISON'
for run in run_names:
    s = s+'_'+run
plt.savefig('%s_STARTYEAR%s.png' % (s,startyear), bbox_inches='tight')

os.chdir(cwd)
