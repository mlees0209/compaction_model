#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 20:35:14 2020

Plot legacy component

@author: mlees
"""

import sys
import pandas as pd
sys.path.append('/home/mlees/InSAR_processing/postprocessing_scripts/')
sys.path.append('/Users/mlees/Documents/RESEARCH/InSAR_processing/postprocessing_scripts/')
import seaborn as sns
from InSAR_postSBAS import *

directory=sys.argv[1]

directory_legacy = directory + '_LEGACY'

print('Reading in data.')
Data = pd.read_csv('%s/Total_Deformation_Out.csv' % directory,parse_dates=[0])
Data_legacy = pd.read_csv('%s/Total_Deformation_Out.csv' % directory_legacy,parse_dates=[0])


print('Rezeroing on 2015.')
Data_rezero = np.array(100 * rezero_series(Data['Total'],np.array([date2num(date) for date in Data['dates']]),'Oct-2015'))
Data_legacyrezero = np.array(100 * rezero_series(Data_legacy['Total'],np.array([date2num(date) for date in Data_legacy['dates']]),'Oct-2015'))

legacypc = 100. *  Data_legacyrezero[np.argwhere(date2num(Data_legacy['dates'].values)==date2num(date(2020,1,1)))][0][0] / Data_rezero[np.argwhere(date2num(Data['dates'].values)==date2num(date(2020,1,1)))][0][0]

print('Legacy subsidence accounts for %.2f %% of total subsidence.' % legacypc )
                                                        
print('Plotting')
sns.set_context('poster')
plt.figure(figsize=(18,12))
plt.plot_date([date2num(date) for date in Data['dates']],Data_rezero,'-',label='Total')
plt.plot_date([date2num(date) for date in Data_legacy['dates']],Data_legacyrezero,'-',label='Legacy')

plt.xlim([date(2015,1,1),date(2020,1,1)])

#yrange = Data_rezero[np.where([date2num(date) for date in Data['dates']]==date2num(date(2015,10,1)))[0][0]] - Data_rezero[np.where([date2num(date) for date in Data['dates']]==date2num(date(2021,10,1)))[0][0]]
plt.ylim([np.min(Data_legacyrezero),-0.25 * np.min(Data_legacyrezero)])

plt.ylabel('Subsidence (cm)')

plt.legend()

plt.title('Legacy subsidence accounts for %.2f %%; latent accounts for %.2f %%' % (legacypc,100-legacypc))

plt.savefig('%s/figures/legacy_subsidence_plot.png' % directory,bbox_inches='tight')