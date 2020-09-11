#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 12:23:45 2020

@author: mlees
"""

from netCDF4 import Dataset
Dat = Dataset("Upper_Aquifer_5.0clay_head_data.nc", "r", format="CF-1.7")
head = Dat.variables['z'][:]