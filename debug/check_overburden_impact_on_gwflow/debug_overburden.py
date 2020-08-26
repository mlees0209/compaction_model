#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 17:44:06 2020

Does plotting etc to check the impact of allowing overburden stress to vary.

@author: mlees
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt
from matplotlib.dates import date2num
from matplotlib.dates import num2date

save=False

lower_aquifer_overburdenoff = np.genfromtxt('debug_OverburdenOFFGwflow_elastic/head_outputs/Lower_Aquifer_10.0clay_head_data.csv',delimiter=',')
lower_aquifer_overburdenon = np.genfromtxt('debug_OverburdenONGwflow_elastic/head_outputs/Lower_Aquifer_10.0clay_head_data.csv',delimiter=',')
lower_aq_elasticflag = np.genfromtxt('debug_OverburdenONGwflow_elastic/Lower_Aquifer_inelastic_flag_DEBUG.csv',delimiter=',')
lower_aq_elasticflag_no = np.genfromtxt('debug_OverburdenOFFGwflow_elastic/Lower_Aquifer_inelastic_flag_DEBUG.csv',delimiter=',')

#%%

dates_loweraq = np.arange(7.12356000e+05,7.37502000e+05,0.25)
input_upper = np.genfromtxt('debug_OverburdenOFFGwflow_elastic/input_head_data/input_time_series_Upper_Aquifer.csv',delimiter=',')

input_lower = np.genfromtxt('debug_OverburdenOFFGwflow_elastic/input_head_data/input_time_series_Lower_Aquifer.csv',delimiter=',')

#%% Prediction with a simple equation; for whatever reason, this didn't seem to work.......perhaps it's the lateral flow which does in fact differ between the two simuiations.



rho = 1000
g=9.81
dsigma = rho*g*input_upper[:,1]
dsigma = dsigma - dsigma[0]
porosity = 0.3
alpha_w = 1.5e-6 / (rho*g*porosity)
alpha = 1.5e-5/ (rho*g) - (alpha_w * porosity)
predicted_dh = 1/(rho * g) * (alpha/alpha_w) * (1/porosity - 1) * dsigma 
#predicted_dh = predicted_dh - predicted_dh[0]

#remainder = lower_aquifer_overburdenon[node,:]- lower_aquifer_overburdenoff[node,:]
#remainder_byinput = remainder[np.where(np.isin(dates_loweraq,input_upper[:,0]))] / (input_upper[1:,1]- input_upper[0,1])

#%%

node=21


plt.figure(figsize=(18,12))
plt.plot_date(np.array(dates_loweraq)[np.where(lower_aq_elasticflag[node,1:])],np.array(lower_aquifer_overburdenon[node,1:])[np.where(lower_aq_elasticflag[node,1:])],'r.',label='overburden on (inelastic)',markersize=0.5)
plt.plot_date(np.array(dates_loweraq)[np.where(1 - lower_aq_elasticflag[node,1:])],np.array(lower_aquifer_overburdenon[node,1:])[np.where(1 - lower_aq_elasticflag[node,1:])],'r^',label='overburden on (elastic)',markersize=1)
plt.plot_date(np.array(dates_loweraq)[np.where(lower_aq_elasticflag_no[node,1:])],np.array(lower_aquifer_overburdenoff[node,1:])[np.where(lower_aq_elasticflag_no[node,1:])],'b.',label='overburden off (inelastic)',markersize=0.5)
plt.plot_date(np.array(dates_loweraq)[np.where(1 -lower_aq_elasticflag_no[node,1:])],np.array(lower_aquifer_overburdenoff[node,1:])[np.where(1 - lower_aq_elasticflag_no[node,1:])],'b^',label='overburden off (elastic)',markersize=1)

plt.ylabel('Head in clay layer (m)')

#plt.plot_date(dates_loweraq[np.where(np.isin(dates_loweraq,input_upper[:,0]))],lower_aquifer_overburdenoff[node,:][np.where(np.isin(dates_loweraq,input_upper[:,0]))] + predicted_dh[1:],label='Prediction of overburden on')


lines, labels = plt.gca().get_legend_handles_labels()

ax2 = plt.gca().twinx()
ax2.plot_date(input_upper[:,0],input_upper[:,1],'k--',label='upper aq head')
plt.ylabel('Head in unconfined aquifer (m)')
#ax2.plot_date(input_lower[:,0],input_lower[:,1],'k-.',label='lower aquifer head')
lines2, labels2 = ax2.get_legend_handles_labels()

ax2.legend(lines + lines2, labels + labels2, loc=[0.71,0.21],fancybox=True)

if save:
    plt.savefig('BasicIllustration.pdf')

plt.xlim(date2num([dt.strptime('Jan-1970','%b-%Y').date(),dt.strptime('Jul-1980','%b-%Y').date()]))

if save:
    plt.savefig('BasicIllustration_ZOOM.pdf')

#%% Does lateral flow differ in the two cases, or is it a  simple addition?

plt.figure(figsize=(18,12))

no_timesteps=5

colors=['b','g','r','c','m','y','k']
j=0

for idx in np.arange(1,len(lower_aq_elasticflag_no[0,:]),len(lower_aq_elasticflag_no[0,:])/no_timesteps):
    print(int(idx))
    plt.plot(lower_aquifer_overburdenoff[:,int(idx)],'-',c=colors[j],label='no overburden t = %i' % int(idx))
    plt.plot(lower_aquifer_overburdenon[:,int(idx)],'--',c=colors[j],label='overburden t = %i' % int(idx))
    j+=1

plt.legend(ncol=no_timesteps)
plt.xlabel('Node')
plt.ylabel('Head (m)')

if save:
    plt.savefig("BasicIllustration_Timeevolution.pdf")
    
    
    
#%% Now pull in the elastic/inelastic version
    
lower_aquifer_overburdenon_elasticinelastic = np.genfromtxt('debug_OverburdenONGwflow_elasticinelastic/head_outputs/Lower_Aquifer_10.0clay_head_data.csv',delimiter=',')
lower_aq_elasticflag_elasticinelastic = np.genfromtxt('debug_OverburdenONGwflow_elasticinelastic/Lower_Aquifer_inelastic_flag_DEBUG.csv',delimiter=',')

lower_aquifer_overburdenoff_elasticinelastic = np.genfromtxt('debug_OverburdenOFFGwflow_elasticinelastic/head_outputs/Lower_Aquifer_10.0clay_head_data.csv',delimiter=',')
lower_aq_elasticflag_elasticinelastic_noverburden = np.genfromtxt('debug_OverburdenOFFGwflow_elasticinelastic/Lower_Aquifer_inelastic_flag_DEBUG.csv',delimiter=',')



#%%

node=5

plt.figure(figsize=(18,12))
#plt.plot_date(dates_loweraq,lower_aquifer_overburdenon[node,1:],label='overburden on, all elastic',markersize=1)
plt.plot_date(np.array(dates_loweraq)[np.where(lower_aq_elasticflag_elasticinelastic[node,1:])],lower_aquifer_overburdenon_elasticinelastic[node,1:][np.where(lower_aq_elasticflag_elasticinelastic[node,1:])],'^',c=colors[1],label='overburden on, elastic inelastic on; inelastic flag',markersize=1)
plt.plot_date(np.array(dates_loweraq)[np.where(1-lower_aq_elasticflag_elasticinelastic[node,1:])],lower_aquifer_overburdenon_elasticinelastic[node,1:][np.where(1-lower_aq_elasticflag_elasticinelastic[node,1:])],'.',c=colors[1],label='overburden on, elastic inelastic on; elastic flag',markersize=0.5)

plt.plot_date(np.array(dates_loweraq)[np.where(lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],lower_aquifer_overburdenoff_elasticinelastic[node,1:][np.where(lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],'^',c=colors[2],label='overburden off, elastic inelastic on; inelastic flag',markersize=1)
plt.plot_date(np.array(dates_loweraq)[np.where(1-lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],lower_aquifer_overburdenoff_elasticinelastic[node,1:][np.where(1-lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],'.',c=colors[2],label='overburden off, elastic inelastic on; elastic flag',markersize=0.5)

plt.ylabel("Head in clay node %i (m)" % node)

lines, labels = plt.gca().get_legend_handles_labels()

ax2 = plt.gca().twinx()
ax2.plot_date(input_upper[:,0],input_upper[:,1],'k--',label='upper aq head')
plt.ylabel('Head in unconfined aquifer (m)')
#ax2.plot_date(input_lower[:,0],input_lower[:,1],'k-.',label='lower aquifer head')
lines2, labels2 = ax2.get_legend_handles_labels()


ax2.legend(lines + lines2, labels + labels2, markerscale=8,ncol=3,fancybox=True)

if save:
    plt.savefig('BasicIllustration_ElasticInelastic_n5.pdf')

node=21

plt.figure(figsize=(18,12))
plt.plot_date(np.array(dates_loweraq)[np.where(lower_aq_elasticflag_elasticinelastic[node,1:])],lower_aquifer_overburdenon_elasticinelastic[node,1:][np.where(lower_aq_elasticflag_elasticinelastic[node,1:])],'^',c=colors[1],label='overburden on, elastic inelastic on; inelastic flag',markersize=1)
plt.plot_date(np.array(dates_loweraq)[np.where(1-lower_aq_elasticflag_elasticinelastic[node,1:])],lower_aquifer_overburdenon_elasticinelastic[node,1:][np.where(1-lower_aq_elasticflag_elasticinelastic[node,1:])],'.',c=colors[1],label='overburden on, elastic inelastic on; elastic flag',markersize=0.5)

plt.plot_date(np.array(dates_loweraq)[np.where(lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],lower_aquifer_overburdenoff_elasticinelastic[node,1:][np.where(lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],'^',c=colors[2],label='overburden off, elastic inelastic on; inelastic flag',markersize=1)
plt.plot_date(np.array(dates_loweraq)[np.where(1-lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],lower_aquifer_overburdenoff_elasticinelastic[node,1:][np.where(1-lower_aq_elasticflag_elasticinelastic_noverburden[node,1:])],'.',c=colors[2],label='overburden off, elastic inelastic on; elastic flag',markersize=0.5)

plt.ylabel("Head in clay node %i (m)" % node)

lines, labels = plt.gca().get_legend_handles_labels()

ax2 = plt.gca().twinx()
ax2.plot_date(input_upper[:,0],input_upper[:,1],'k--',label='upper aq head')
plt.ylabel('Head in unconfined aquifer (m)')
#ax2.plot_date(input_lower[:,0],input_lower[:,1],'k-.',label='lower aquifer head')
lines2, labels2 = ax2.get_legend_handles_labels()


ax2.legend(lines + lines2, labels + labels2, markerscale=8,ncol=3,fancybox=True)

if save:
    plt.savefig('BasicIllustration_ElasticInelastic_n21.pdf')
    
    
#%% Check that effective stress is going ok
import matplotlib.dates as mdates

lower_aquifer_overburdenon_elasticinelastic_overburden = np.genfromtxt('debug_OverburdenONGwflow_elasticinelastic_Vthick/Lower_Aquifer_10.0clay_overburden_stress.csv',delimiter=',')

lower_aquifer_overburdenon_elasticinelastic_effectivestress = np.genfromtxt('debug_OverburdenONGwflow_elasticinelastic_Vthick/Lower_Aquifer_10.0clay_effective_stress.csv',delimiter=',')


#%% Check things are going well with effective stress. 

x_lims = list(map(dt.fromordinal,[int(min(dates_loweraq)),int(max(dates_loweraq))]))
x_lims = date2num(x_lims)
y_lims = [0,10]


#y_lims=rho* g *np.array([np.min(lower_aquifer_overburdenon_elasticinelastic - 40.9233),np.max(lower_aquifer_overburdenon_elasticinelastic- 40.9233)])


f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1,sharex=True)
a= ax1.imshow(rho * g * (lower_aquifer_overburdenon_elasticinelastic - 40.9233),aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
plt.gca().xaxis_date()
date_format = mdates.DateFormatter('%Y')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gcf().autofmt_xdate()
ax1.set_title('Pressure Head')
plt.colorbar(a,ax=ax1,label='rho g h (Nm^-2)')

#y_lims=[np.min(lower_aquifer_overburdenon_elasticinelastic_overburden),np.max(lower_aquifer_overburdenon_elasticinelastic_overburden)]
b= ax2.imshow(lower_aquifer_overburdenon_elasticinelastic_overburden,aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
ax2.set_title('Overburden stress')
plt.colorbar(b,ax=ax2,label='overburden stress (Nm^-2)')

#y_lims=[np.min(lower_aquifer_overburdenon_elasticinelastic_effectivestress),np.max(lower_aquifer_overburdenon_elasticinelastic_effectivestress)]
c= ax3.imshow(lower_aquifer_overburdenon_elasticinelastic_effectivestress,aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
ax3.set_title('Effective stress (from model)')
plt.colorbar(c,ax=ax3,label='effective stress (Nm^-2)')

effective_stress_calc = lower_aquifer_overburdenon_elasticinelastic_effectivestress - rho * g * (lower_aquifer_overburdenon_elasticinelastic - 40.9233)
#y_lims=[np.min(effective_stress_calc),np.max(effective_stress_calc)]
d= ax4.imshow(effective_stress_calc,aspect='auto',extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]])
ax4.set_title('Effective stress (from data)')
plt.colorbar(d,ax=ax4,label='effective stress (Nm^-2)')
