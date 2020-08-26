#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 18:47:23 2020
Compares mine and Ryan's models
@author: mlees
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import date2num

save=True

modflow_sub = np.genfromtxt('MODFLOW/modflow_subsidence_data_output_daily.csv')
modflow_dates = np.genfromtxt('MODFLOW/modflow_subsidence_dates_output_daily.csv')



ML_sub = pd.read_csv('ML_Model_Benchmark/Total_Deformation_Out.csv',parse_dates=[0])

#%%

plt.figure(figsize=(18,12))
plt.plot_date(modflow_dates+date2num(ML_sub['dates'][0]),-modflow_sub,'--',label='MODFLOW')
plt.plot_date(date2num(ML_sub['dates']),ML_sub['Total'],'-',label='ML model')
plt.legend()
plt.xlabel('Date')
plt.ylabel("Subsidence (m)")

if save:
    plt.savefig('Benchmark_Output_Comparisons.pdf',bbox_inches='tight')