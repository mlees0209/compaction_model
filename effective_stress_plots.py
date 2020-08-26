#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:22:05 2020
Checks that things are going well with effective stress calculations (useful when including overburden)
@author: mlees
"""

# Usage: python effective_stress_plots DIR 
# DIR is the output directory of the run where I want to make these plots.




#%% Import some data... 

import sys
import os
import numpy as np
import re
import pandas as pd
from matplotlib.dates import date2num
from datetime import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import csv

save=True
rho=1000
g=9.81

if len(sys.argv) != 2:
    print('effective_stress_plots error; terminal. Incorrect number of input arguments. Correct usage: python effective_stress_plots DIR')
    sys.exit(1)

directory=sys.argv[1]
#directory='/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/debug/check_overburden_impact_on_gwflow/debug_OverburdenONGwflow_elasticinelastic/'
if directory[-1] != '/':
    directory+='/'

print('Finding stress output files.')
effective_stress_filenames_tmp = [directory+file for file in os.listdir(directory) if 'effective_stress' in file]
#head_filenames_tmp = [directory+'/head_outputs/'+file for file  in os.listdir(directory+'/head_outputs') if 'head_data' in file]
overburden_stress_filenames_tmp = [directory+file for file in os.listdir(directory) if 'overburden_stress' in file]
inelastic_flag_filenames_tmp = [directory+file for file in os.listdir(directory) if 'inelastic_flag' in file]

if len(effective_stress_filenames_tmp) == len(overburden_stress_filenames_tmp):
    head_filenames_tmp = [directory+'head_outputs/'+file for file in os.listdir(directory+'head_outputs') if (any([ word.replace('_effective_stress.csv','').replace(directory,'') in file for word in effective_stress_filenames_tmp])) and (not 'groundwater_solution_dates' in file)]
else:
    print("Error, number of effective_stress and overburden_stress filenames don't match up.")
    sys.exit(1)
if len(head_filenames_tmp)==len(effective_stress_filenames_tmp):
    layernames = [name.replace('_head_data.csv','').replace(directory,'').replace('head_outputs/','') for name in head_filenames_tmp]
    print('Files found. Stress output files found for the following layers. These will be plotted:')
    print(layernames)
else:
    print("Error, number of effective_stress and head filenames don't match up.")
    sys.exit(1)

head_filenames={layer:fileloc for (layer,fileloc) in zip(layernames,head_filenames_tmp)}

effective_stress_filenames={}
for layer in layernames:
    effective_stress_filenames[layer] = np.array(effective_stress_filenames_tmp)[np.where([layer in bs for bs in effective_stress_filenames_tmp])[0]][0]
#overburden_stress_filenames={layer:fileloc for (layer,fileloc) in zip(layernames,overburden_stress_filenames)}

overburden_stress_filenames={}
for layer in layernames:
    overburden_stress_filenames[layer] = np.array(overburden_stress_filenames_tmp)[np.where([ layer in bs for bs in overburden_stress_filenames_tmp])[0]][0]


inelastic_flag_filenames={}
for layer in layernames:
    inelastic_flag_filenames[layer] = np.array(inelastic_flag_filenames_tmp)[np.where([ layer in bs for bs in inelastic_flag_filenames_tmp])[0]][0]

#%% For each layer, make a 4-panel plot showing the effective stress thing

print('')
for layer in layernames:
    print('Importing data for '+layer)
    head = np.genfromtxt(head_filenames[layer],delimiter=',') 
    head = head - head[0,0]
    overburden = np.genfromtxt(overburden_stress_filenames[layer],delimiter=',') 
    overburden = overburden - overburden[0,0]
    effective_stress = np.genfromtxt(effective_stress_filenames[layer],delimiter=',') 
    effective_stress = effective_stress - np.min(effective_stress[:,0]) 
    inelastic_flag = np.genfromtxt(inelastic_flag_filenames[layer],delimiter=',')
    inelastic_idxs = np.argwhere(inelastic_flag==1)
    
    vmin=np.min([np.min(overburden),np.min(effective_stress),rho*g*np.min(head)])
    vmax=np.max([np.min(overburden),np.max(effective_stress),rho*g*np.max(head)])
    
    print('\tReading the dates')
    datesfilename = re.sub('[0-9]+', '', layer)
    datesfilename = datesfilename.replace('_.clay','').replace('.clay','')
    datesfilename=datesfilename+'_groundwater_solution_dates.csv'
    dates_tmp=[]
    with open(directory+'head_outputs/'+datesfilename) as csvfile:
        datesreader = csv.reader(csvfile, delimiter=',')
        for row in datesreader:
            dates_tmp.append(row)
    dates_tmp = date2num([pd.to_datetime(date) for date in dates_tmp[0]])

    print('\tData imported')
    print('\tPlotting 4-panel plot')
    
    x_lims = list(map(dt.fromordinal,[int(min(dates_tmp)),int(max(dates_tmp))]))
    x_lims = date2num(x_lims)
    x_values = np.linspace(x_lims[0],x_lims[1],num=np.shape(head)[1])
    y_lims = [0,np.shape(head)[0]]
        
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1,sharex=True,figsize=(12,15))
    f.suptitle(layer) # Make title for whole thing
    
    a= ax1.imshow(rho * g * head,aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]],vmin=vmin,vmax=vmax)
    plt.gca().xaxis_date()
    date_format = mdates.DateFormatter('%Y')
    plt.gca().xaxis.fmajor_formatter(date_format)
    plt.gcf().autofmt_xdate()
    ax1.set_title('Pressure Head')
#    plt.colorbar(a,ax=ax1,label='rho g h (Nm^-2)')
    
    #y_lims=[np.min(lower_aquifer_overburdenon_elasticinelastic_overburden),np.max(lower_aquifer_overburdenon_elasticinelastic_overburden)]
    b= ax2.imshow(overburden,aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]],vmin=vmin,vmax=vmax)
    ax2.set_title('Overburden stress')
#    plt.colorbar(b,ax=ax2,label='overburden stress (Nm^-2)')
    
    #y_lims=[np.min(lower_aquifer_overburdenon_elasticinelastic_effectivestress),np.max(lower_aquifer_overburdenon_elasticinelastic_effectivestress)]

    
    
    c= ax3.imshow(effective_stress,aspect='auto',extent = [x_lims[0], x_lims[1],0,np.shape(effective_stress)[0]],vmin=vmin,vmax=vmax)
    #ax3.imshow(effective_stress_inelastic,aspect='auto',extent = [x_lims[0], x_lims[1],0,np.shape(effective_stress)[0]],vmin=vmin,vmax=vmax)
    
    [X,Y] = np.meshgrid(x_values,np.linspace(y_lims[0],y_lims[1],num=np.shape(head)[0]))
    cs = ax3.contour(X,Y,inelastic_flag,levels=[0.95],colors='black',linewidths=1,linestyles='dashed')
    
    h1,_ = cs.legend_elements()
    ax3.legend([h1[0]], ['Inelastic regions'],fontsize='x-small',loc='upper left')
    ax3.set_title('Effective stress (from model)')
#    plt.colorbar(c,ax=ax3,label='effective stress (Nm^-2)')
    
    effective_stress_calc = overburden - rho * g * head
    #y_lims=[np.min(effective_stress_calc),np.max(effective_stress_calc)]
    d= ax4.imshow(effective_stress_calc,aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]],vmin=vmin,vmax=vmax)
#    print('drawing around inelastic bits')
#    for coords in inelastic_idxs[0:200]:
#        highlight_cell(x_values[coords[0]],coords[1],ax=ax4,color='black',linewidth=1)
    ax4.set_title('Effective stress (from data)')
#    plt.colorbar(d,ax=ax4,)
    
    f.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.02, hspace=0.1)

# add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8

    cb_ax = f.add_axes([0.83, 0.1, 0.02, 0.8])
    cbar = f.colorbar(d, cax=cb_ax,label=r'effective stress (Nm$^{-2}$)')

    
    if save:
        plt.savefig(directory+'figures/'+'Effective_Stress_Debug_Fourpanel_%s.pdf' % layer,bbox_inches='tight')