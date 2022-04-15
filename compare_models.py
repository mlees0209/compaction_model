#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 16:41:38 2020

Simple script to compare N runs. OPTIONAL INPUT ARGUMENT 'startyear=XXXX' allows the run to be rezerod on a year of your choice; default is Jun 1980 for some reason. OPTIONAL INPUT ARGUMENT 'realdata=TRE_Altamira' adds the Tre Altamira data to the plot.

@author: mlees
"""
import sys

if len(sys.argv) <= 1:
    print('compare_models error; terminal. Incorrect number of input arguments. Correct usage: python compare_models.py model1 model2 model3 ... modelN.')
    sys.exit(1)


import os


import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
import seaborn as sns
from InSAR_postSBAS import *

#no_runs = len(sys.argv)-1

run_names = sys.argv[1:]
startyear = [var for var in run_names if var.split('=')[0]=='startyear']
if not startyear:
    print('No startyear specified; zeroing on Sep 1st 1980.')
    startyear=1980
else:
    startyear=int(startyear[0].split('=')[1])
    print('Startyear specified; zeroing on Sep 1st %i.' % startyear)

realdata= [var for var in run_names if var.split('=')[0]=='realdata']
if realdata:
    if realdata[0].split('=')[1]=='TRE_Altamira':
        print('realdata specified as TRE Altamira; importing.')
        if ospath.exists('/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
            mac=1
            linux=0
            print("\tWe're in Mac, looking for downloaded InSAR data accordingly.")
            file='/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
            InSAR_data = import_InSAR_csv(file)
            WellHlat = 36.32750
            WellHlon = -119.58056
            datesinsarH,InSAR_H = extract_series_from_latlon(WellHlat,WellHlon,InSAR_data)
        elif ospath.exists('/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
            linux=1
            mac=0
            print("\tWe're in Linux, looking for downloaded InSAR data accordingly.")
            file='/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
            InSAR_data = import_InSAR_csv(file)
            WellHlat = 36.32750
            WellHlon = -119.58056
            datesinsarH,InSAR_H = extract_series_from_latlon(WellHlat,WellHlon,InSAR_data)

        else:
            print("\tUnable to find downloaded InSAR data. Check for bigdata/InSAR folder. Ignoring.")
            
    else:
        print('realdata specified but unknown option. Try TRE_Altamira')

set_context = [var for var in run_names if var.split('=')[0]=='set_context']
if not set_context:
    context='poster'
else:
    context=set_context[0].split('=')[1]
    print('\tset_context specified as %s. Setting sns.set_context accordingly.' % context)

run_names = [var for var in run_names if var.split('=')[0]!='startyear']
run_names = [var for var in run_names if var.split('=')[0]!='realdata']
run_names = [var for var in run_names if var.split('=')[0]!='set_context']

print("Plotting for runs ", run_names)

cwd = os.getcwd()
Data={}

sns.set_context(context)
sns.set_style('darkgrid')

fig,ax1 = plt.subplots(figsize=(16,8))


for run in run_names:
    print('\tPlotting %s.' % run)
    os.chdir(run)
    Data[run] = pd.read_csv('Total_Deformation_Out.csv',parse_dates=[0])
#    ax1.plot_date([date2num(date) for date in Data[run]['dates'] if date < date2num(datetime(1965,1,1))],100*rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[0:365*20],'-',color='grey',linewidth=0.5,label='_nolegend_')

    ax1.plot_date([date2num(date) for date in Data[run]['dates']],100*rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear),'-',label=run)
    os.chdir(cwd)

if realdata:
    if realdata[0].split('=')[1]=='TRE_Altamira':
        print('\tAdding InSAR.')
        ax1.plot_date(datesinsarH,0.1*InSAR_H + 100* rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[np.where([date2num(date) for date in Data[run]['dates']]==datesinsarH[0])[0][0]],'k-',label='Sentinel InSAR subsidence',linewidth=0.5)
        upperbound_tmp = 0.1*InSAR_H + 100* rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[np.where([date2num(date) for date in Data[run]['dates']]==datesinsarH[0])[0][0]] + 1.6        
        lowerbound_tmp = 0.1*InSAR_H + 100* rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[np.where([date2num(date) for date in Data[run]['dates']]==datesinsarH[0])[0][0]] - 1.6
        ax1.fill_between(datesinsarH, lowerbound_tmp, upperbound_tmp,facecolor='grey',alpha=0.5,edgecolor='None',label='InSAR uncertainty',zorder=200)
        ax1.set_ylabel('between y1 and 0')



plt.legend()
plt.ylabel('Deformation (cm)')
s='COMPARISON'
for run in run_names:
    s = s+'_'+run
plt.savefig('%s_STARTYEAR%s.png' % (s,startyear), bbox_inches='tight')
plt.savefig('%s_STARTYEAR%s.pdf' % (s,startyear), bbox_inches='tight')
plt.savefig('%s_STARTYEAR%s.svg' % (s,startyear), bbox_inches='tight')

# Do 2015-2020 version
#yrange = 10
#plt.ylim([100*rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[np.argwhere(date2num(Data[run]['dates']) == date2num(dt(2019,6,1)))[0][0]]-yrange,100*rezero_series(Data[run]['Total'],np.array([date2num(date) for date in Data[run]['dates']]),'Sep-%s' % startyear)[np.argwhere(date2num(Data[run]['dates']) == date2num(dt(2015,1,1)))[0][0]]+yrange])
    
#plt.xlim([date(2015,1,1),date(2020,1,1)])
#plt.savefig('%s_STARTYEAR%s_20152020.pdf' % (s,startyear), bbox_inches='tight')


os.chdir(cwd)
