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

Sske=1.5*10**(-5)
Sskv=4.6*10**(-4)
dz = 2
dh0 = hmat[:,1] - hmat[:,0]
ds0 = dz*( np.dot(inelastic_flag[:,0] * Sskv, dh0) + np.dot(~inelastic_flag[:,0] * Sske, dh0))

def process(token):
    return token['text']


ds = [dz*( np.dot(inelastic_flag[:,i] * Sskv, hmat[:,i+1] - hmat[:,i]) + np.dot(~inelastic_flag[:,i] * Sske, hmat[:,i+1] - hmat[:,i])) for i in range(len(t)-1)]
s = [np.sum(ds[:i]) for i in range(len(t)-1)]

plt.plot_date(t[1:],s)