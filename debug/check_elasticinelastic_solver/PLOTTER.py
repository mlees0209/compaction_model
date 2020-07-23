#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 20:01:51 2020

@author: mlees
"""

import numpy as np
import pandas as pd
import csv
from matplotlib.dates import date2num
from datetime import datetime as dt
import os
import matplotlib.pyplot as plt
import seaborn as sns

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()



reader = csv.reader(open('Corcoran_Clay_head_data.csv', "r"), delimiter=",")
x = list(reader)
hmat = np.array(x).astype("float")

reader = csv.reader(open('Corcoran_Clay_inelastic_flag_DEBUG.csv', "r"), delimiter=",")
x = list(reader)
inelastic_flag = np.array(x).astype("float")
inelastic_flag = inelastic_flag==1

with open('Corcoran_Clay_groundwater_solution_dates.csv', newline='') as f:
    reader = csv.reader(f)
    dates_str = list(reader)[0]

t = date2num([dt.strptime(date, '%d-%b-%Y').date() for date in dates_str])
if not os.path.isdir('video_frames'):
    os.mkdir('video_frames')

# Make the video frames; use ffmpeg -f image2 -i %*.png vid.mp4 to make the vid itself
plt.ioff()
sns.set_context('talk')
for i in np.arange(0,len(t),50):
    printProgressBar(i,len(t))
    plt.figure(figsize=(12,12))
    plt.plot(hmat[:,0],np.linspace(0,40,20),'k--',label='t=0 position')
    plt.plot(np.min(hmat[:,:i+1],axis=1),np.linspace(0,40,20),'b--',label='Min head')
    plt.plot(hmat[:,i][inelastic_flag[:,i]],np.linspace(0,40,20)[inelastic_flag[:,i]],'r.')
    plt.plot(hmat[:,i][~inelastic_flag[:,i]],np.linspace(0,40,20)[~inelastic_flag[:,i]],'g.')

    plt.title('t=%.0f' % (np.round(i/365) + 1970))
    plt.xlabel('Head in the clay (m)')
    plt.ylabel('Z (m)')
    plt.xlim([60,-30])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig('video_frames/frame%05d.png' % i,bbox_inches='tight')
    plt.close()
