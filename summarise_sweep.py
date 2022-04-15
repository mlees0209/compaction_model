#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:11:42 2021

Creates a big csv file which can be used to assess the parameter space sweep

@author: mlees
"""

import sys
import os

import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
sys.path.append('/home/mlees/InSAR_postprocessing/')
sys.path.append('/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/')
sys.path.append('/home/mlees/Land_Subsidence/Local_Scale/MODEL/')

import seaborn as sns
from InSAR_postSBAS import *
from model_functions import *

import datetime as datetime

sweep_directory='/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/Model_Runs/ParamSpaceSweep_Apr1'

foot_to_cm = 30.48

#%% Import Sentinel InSAR

if os.path.exists('/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    mac=1
    linux=0
    knightblade=0
    print("\tWe're in Mac, looking for downloaded InSAR data accordingly.")
    file='/Users/mlees/Documents/RESEARCH/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
elif ospath.exists('/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    linux=1
    mac=0
    knightblade=0
    print("\tWe're in Linux, looking for downloaded InSAR data accordingly.")
    file='/home/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
elif ospath.exists('/data1/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical'):
    print("\tWe're in Knightblade, looking for downloading InSAR data accordingly.")
    mac=0
    linux=0
    knightblade=1
    file='/data1/mlees/bigdata/InSAR/Processed_datasets/TRE_Altamira_Vertical/CALIFORNIA_DWR_studyarea.csv'
else:
    print("\tUnable to find downloaded InSAR data. Check for bigdata/InSAR folder. Aborting.")
    sys.exit()


InSAR_data = import_InSAR_csv(file)
WellHlat = 36.32750
WellHlon = -119.58056

datesinsarH,InSAR_H = extract_series_from_latlon(WellHlat,WellHlon,InSAR_data)

datesinsarH_primaryperiod = [datey for datey in datesinsarH if (datey>=date2num(datetime.date(2015,8,1))) & (datey <=date2num(datetime.date(2019,9,1)))]

Sentinel_rezeroed = np.array(0.1 * rezero_series(InSAR_H,datesinsarH,'Aug-2015'))

Sentinel_rezeroed_primary = Sentinel_rezeroed[np.isin(datesinsarH,datesinsarH_primaryperiod)]

Poland_75_dates = date2num([dt(1954,1,1),dt(1957,1,1),dt(1958,7,1),dt(1962,1,1),dt(1966,4,1),dt(1970,1,1)])
Poland_75_data = -30.48*  np.array([0,0.83,1.34,3.25,4.92,6.31])
Poland_6670_rate = 365 * (Poland_75_data[-1] - Poland_75_data[-2]) / (Poland_75_dates[-1] - Poland_75_dates[-2])

Highway_198_dates = date2num([date(1972,6,1),date(2004,1,1)])

#%%
#directory='/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/Model_Runs/Parameter_Space_Sweep_Jan2121'


model_runs= next(os.walk(sweep_directory))[1]
model_runs = [run for run in model_runs if run.startswith('C')]

results = pd.DataFrame(columns=['RunID','Runname','Clay','Sskv','Sske','Kv','d_i','precons','Poland 66-70 rate (cm/yr)','Highway 198 modelled subsidence (cm)','ALOS subsidence rate (cm/yr)','Period15-17 average rate (cm/yr)','Period07-09 average rate (cm/yr)'])

#run = model_runs[0]

i=0

for run in model_runs:
    print('Reading layer names...')
    param_filename = 'paramfile.par'
    f = open("%s/%s/%s" % (sweep_directory,run,param_filename), "r")
    paramfilelines = f.readlines()
    f.close()
    paramfilelines = [line.strip() for line in paramfilelines]
    paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
    paramfilelines[:] = [x for x in paramfilelines if x]
    no_layers = read_parameter('no_layers',int,1,paramfilelines)
    layer_names=read_parameter('layer_names',str,no_layers,paramfilelines)
    layer_types=read_parameter('layer_types',str,no_layers,paramfilelines)
    aquitards = [name for name,value in layer_types.items() if value=='Aquitard']
    interbeds_switch=read_parameter('interbeds_switch',bool,list(layer_types.values()).count('Aquifer'),paramfilelines)
    interbedded_layers= [name for name,value in interbeds_switch.items() if value==True]    
    layers_requiring_solving = interbedded_layers + aquitards
    
    groundwater_flow_solver_type=read_parameter('groundwater_flow_solver_type',str,len(layers_requiring_solving),paramfilelines)
    clay_Sse = read_parameter('clay_Sse',float,sum(groundwater_flow_solver_type[layer]=='elastic-inelastic' or compaction_solver_compressibility_type[layer]=='elastic-inelastic' for layer in layer_names),paramfilelines)
    clay_Ssv = read_parameter('clay_Ssv',float,sum(groundwater_flow_solver_type[layer]=='elastic-inelastic' or compaction_solver_compressibility_type[layer]=='elastic-inelastic' for layer in layer_names),paramfilelines)
    clay_value = run.split('y')[1].split('S')[0]  
    di_value = run.split('di')[1].split('pre')[0]
    if 'preconshead' in run:
        precons_value = run.split('preconshead')[1]
    else:
        precons_value=0
    
    vertical_conductivity = read_parameter('vertical_conductivity',float,len(layers_requiring_solving),paramfilelines)

    if os.path.isfile('%s/%s/Total_Deformation_Out.csv' % (sweep_directory,run)):
        Deformation_Out_tmp = pd.read_csv('%s/%s/Total_Deformation_Out.csv' % (sweep_directory,run),parse_dates=[0])
    
        Deformation_Out_primaryperiodtmp = Deformation_Out_tmp['Total'][np.isin(date2num(Deformation_Out_tmp['dates']),datesinsarH_primaryperiod)]
        Deformation_Out_primaryperiodtmp = rezero_series(Deformation_Out_primaryperiodtmp.values,datesinsarH_primaryperiod,'Aug-2015')
        
        # primary_RMS = np.sqrt(np.mean( (100*Deformation_Out_primaryperiodtmp-Sentinel_rezeroed_primary)**2))
        # print('Sentinel RMS = %.2f cm' % primary_RMS)
        PolandPeriod_dates = np.arange(Poland_75_dates[-2],Poland_75_dates[-1],step=1)
        PolandPeriod_modelled_average_rate = 365 * np.polyfit(PolandPeriod_dates,100*Deformation_Out_tmp['Total'][np.isin(date2num(Deformation_Out_tmp['dates']),PolandPeriod_dates)],1)[0] # get the linear modelled rate over the envisat period.
        print('1966-1970 average rate = %.2f cm/yr' % PolandPeriod_modelled_average_rate)


        Period1517_dates = np.arange(date2num(datetime.date(2015,3,1)),date2num(datetime.date(2017,3,1)),step=1)
        Period1517_modelled_average_rate = 365 * np.polyfit(Period1517_dates,100*Deformation_Out_tmp['Total'][np.isin(date2num(Deformation_Out_tmp['dates']),Period1517_dates)],1)[0] # get the linear modelled rate over the envisat period.
        print('2015-17 average rate = %.2f cm/yr' % Period1517_modelled_average_rate)

        
        deformation_modelled_Highway198_tmp =100 * (Deformation_Out_tmp['Total'][date2num(Deformation_Out_tmp['dates'])==Highway_198_dates[1]].values - Deformation_Out_tmp['Total'][date2num(Deformation_Out_tmp['dates'])==Highway_198_dates[0]].values)
        
        print('highway 198 sub = %.2f cm' % deformation_modelled_Highway198_tmp)

        # # Plot them to check
        # plt.figure()
        # plt.plot_date(datesinsarH_primaryperiod,Sentinel_rezeroed_primary,'r-')
        # plt.plot_date(datesinsarH_primaryperiod,100*Deformation_Out_primaryperiodtmp,'b-')
        Deformation_Out_ALOStmp = Deformation_Out_tmp['Total'][(date2num(Deformation_Out_tmp['dates']) >= date2num(datetime.date(2007,6,1))) & (date2num(Deformation_Out_tmp['dates']) <= date2num(datetime.date(2010,7,2)))]
        Deformation_OUT_ALOS_dates = date2num(Deformation_Out_tmp['dates'][(date2num(Deformation_Out_tmp['dates']) >= date2num(datetime.date(2007,6,1))) & (date2num(Deformation_Out_tmp['dates']) <= date2num(datetime.date(2010,7,2)))])
        
        ALOS_modelled_average_rate = 365 * np.polyfit(Deformation_OUT_ALOS_dates,100*Deformation_Out_ALOStmp,1)[0] # get the linear modelled rate over the envisat period.
        print('ALOS average rate = %.2f cm/yr' % ALOS_modelled_average_rate)

        Period0709_dates = np.arange(date2num(datetime.date(2007,9,1)),date2num(datetime.date(2009,9,1)),step=1)

        Period0709_modelled_average_rate = 365 * np.polyfit(Period0709_dates,100*Deformation_Out_tmp['Total'][np.isin(date2num(Deformation_Out_tmp['dates']),Period0709_dates)],1)[0] # get the linear modelled rate over the envisat period.
        print('2007-09 average rate = %.2f cm/yr' % Period0709_modelled_average_rate)


        results.loc[len(results.index)] = [i,run,clay_value,clay_Ssv['Upper Aquifer'],clay_Sse['Upper Aquifer'],vertical_conductivity['Upper Aquifer'],di_value,precons_value,PolandPeriod_modelled_average_rate,deformation_modelled_Highway198_tmp[0],ALOS_modelled_average_rate,Period1517_modelled_average_rate,Period0709_modelled_average_rate]
        print('i=%i' % i)
        print('')
        print('')

        i+=1

    else:
        results.loc[len(results.index)] = [i,run,clay_value,clay_Ssv['Upper Aquifer'],clay_Sse['Upper Aquifer'],vertical_conductivity['Upper Aquifer'],di_value,precons_value,'NA','NA','NA','NA','NA']

        i+=1


results.to_csv('%s/param_sweep_summary.csv' % sweep_directory,index=False)
