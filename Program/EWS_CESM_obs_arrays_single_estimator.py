#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:05:55 2024

Script creates figures of the ratio of EWS along the observational transects (SAMBA, RAPID, OSNAP)

Figure 2 of paper, and Figures S1 - S5

@author: 6008399
"""

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ruptures as rpt
from scipy.interpolate import CubicSpline
from scipy.interpolate import CubicHermiteSpline
import statsmodels.api as sm
import pandas as pd
from pandas.plotting import autocorrelation_plot
from pandas import DataFrame
from sklearn.linear_model import LinearRegression
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import matplotlib.colors as mcolors
import cartopy.mpl.ticker as cticker

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/Cesm_collapse/netcdf/'
directory_figures = r'/Users/6008399/Documents/PhD/Cesm_collapse/Figures/'

#%% Ratio of variance using single estimator

var_salt_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_34S.nc','r')
var_temp_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_34S.nc','r')
       
ratio_var_salt_samba  = var_salt_samba.variables['SIGMA_SALT'][:]
ratio_var_temp_samba  = var_temp_samba.variables['SIGMA_TEMP'][:]

p_var_salt_samba  = var_salt_samba.variables['p_value_LEVENE'][:]
p_var_temp_samba  = var_temp_samba.variables['p_value_LEVENE'][:]

lat_samba = var_salt_samba.variables['lat'][:]
lon_samba = var_salt_samba.variables['lon'][:]
depth_samba = var_salt_samba.variables['depth'][:]

p_var_salt_sig_samba = ma.masked_all((len(depth_samba), len(lon_samba)))
p_var_temp_sig_samba = ma.masked_all((len(depth_samba), len(lon_samba)))

#Create mask where the p-value is smaller than 0.05
for i in range(len(depth_samba)):
    for j in range(len(lon_samba)):
        if p_var_salt_samba[i,j] < 0.05:# and ratio_var_salt_samba[i,j] > 1:
            p_var_salt_sig_samba[i,j] = 1
        else:
            continue
        
        if p_var_temp_samba[i,j] < 0.05:# and ratio_var_temp_samba[i,j] > 1:
            p_var_temp_sig_samba[i,j] = 1
        else:
            continue

var_salt_2_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_34S_branch_1650_1500.nc','r')
var_temp_2_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_34S_branch_1650_1500.nc','r')
      
ratio_var_salt_2_samba  = var_salt_2_samba.variables['SIGMA_SALT'][:]
ratio_var_temp_2_samba  = var_temp_2_samba.variables['SIGMA_TEMP'][:]

p_var_salt_2_samba  = var_salt_2_samba.variables['p_value_LEVENE'][:]
p_var_temp_2_samba  = var_temp_2_samba.variables['p_value_LEVENE'][:]

lat_samba = var_salt_2_samba.variables['lat'][:]
lon_samba = var_salt_2_samba.variables['lon'][:]
depth_samba = var_salt_2_samba.variables['depth'][:]

p_var_salt_sig_2_samba = ma.masked_all((len(depth_samba), len(lon_samba)))
p_var_temp_sig_2_samba = ma.masked_all((len(depth_samba), len(lon_samba)))

#Create mask where the p-value is smaller than 0.05
for i in range(len(depth_samba)):
    for j in range(len(lon_samba)):
        if p_var_salt_2_samba[i,j] < 0.05:# and ratio_var_salt_2_samba[i,j] > 1:
            p_var_salt_sig_2_samba[i,j] = 1
        else:
            continue
        
        if p_var_temp_2_samba[i,j] < 0.05:# and ratio_var_temp_2_samba[i,j] > 1:
            p_var_temp_sig_2_samba[i,j] = 1
        else:
            continue

var_salt_3_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_34S_branch_1650_600.nc','r')
var_temp_3_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_34S_branch_1650_600.nc','r')
      
ratio_var_salt_3_samba  = var_salt_3_samba.variables['SIGMA_SALT'][:]
ratio_var_temp_3_samba  = var_temp_3_samba.variables['SIGMA_TEMP'][:]

p_var_salt_3_samba  = var_salt_3_samba.variables['p_value_LEVENE'][:]
p_var_temp_3_samba  = var_temp_3_samba.variables['p_value_LEVENE'][:]

lat_samba = var_salt_3_samba.variables['lat'][:]
lon_samba = var_salt_3_samba.variables['lon'][:]
depth_samba = var_salt_3_samba.variables['depth'][:]

p_var_salt_sig_3_samba = ma.masked_all((len(depth_samba), len(lon_samba)))
p_var_temp_sig_3_samba = ma.masked_all((len(depth_samba), len(lon_samba)))

#Create mask where the p-value is smaller than 0.05
for i in range(len(depth_samba)):
    for j in range(len(lon_samba)):
        if p_var_salt_3_samba[i,j] < 0.05:# and ratio_var_salt_3_samba[i,j] > 1:
            p_var_salt_sig_3_samba[i,j] = 1
        else:
            continue
        
        if p_var_temp_3_samba[i,j] < 0.05:# and ratio_var_temp_3_samba[i,j] > 1:
            p_var_temp_sig_3_samba[i,j] = 1
        else:
            continue

#%% Ratio of AC1 using single estimator
ac_salt_NO_samba     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_SAMBA_noACrestriction.nc','r')
ac_salt_samba     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_34S_branch_1500_600.nc','r')

ac_salt_2_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_34S_branch_1650_1500.nc','r')
ac_salt_3_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_34S_branch_1650_600.nc','r')
ac_salt_test_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_34S_branch_1500_600test.nc','r')

ac_salt_jan_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_34S_monthly_January.nc','r')
ac_salt_jun_samba = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_34S_monthly_June.nc','r')

ac_temp_samba     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_34S_branch_1500_600.nc','r')
ac_temp_2_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_34S_branch_1650_1500.nc','r')
ac_temp_3_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_34S_branch_1650_600.nc','r')
       
lon_month_samba   = ac_salt_jan_samba.variables['lon'][:]

ratio_ac_salt_samba       = ac_salt_samba.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_NO_samba       = ac_salt_NO_samba.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_2_samba     = ac_salt_2_samba.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_3_samba     = ac_salt_3_samba.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_test_samba     = ac_salt_test_samba.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_jan_samba   = ac_salt_jan_samba.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_jun_samba   = ac_salt_jun_samba.variables['RATIO_AC1_SALT'][:]

ratio_ac_temp_samba       = ac_temp_samba.variables['RATIO_AC1_TEMP'][:]
ratio_ac_temp_2_samba     = ac_temp_2_samba.variables['RATIO_AC1_TEMP'][:]
ratio_ac_temp_3_samba     = ac_temp_3_samba.variables['RATIO_AC1_TEMP'][:]

sig_ac_salt_pos_samba     = ac_salt_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_NO_samba     = ac_salt_NO_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_2_samba   = ac_salt_2_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_3_samba   = ac_salt_3_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_test_samba   = ac_salt_test_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_jan_samba = ac_salt_jan_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_jun_samba = ac_salt_jun_samba.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_neg_samba     = ac_salt_samba.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_NO_samba     = ac_salt_NO_samba.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_2_samba   = ac_salt_2_samba.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_3_samba   = ac_salt_3_samba.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_jan_samba = ac_salt_jan_samba.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_jun_samba = ac_salt_jun_samba.variables['SIG_RATIO_AC1_SALT_neg'][:]

sig_ac_temp_pos_samba     = ac_temp_samba.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_pos_2_samba   = ac_temp_2_samba.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_pos_3_samba   = ac_temp_3_samba.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_neg_samba     = ac_temp_samba.variables['SIG_RATIO_AC1_TEMP_neg'][:]
sig_ac_temp_neg_2_samba   = ac_temp_2_samba.variables['SIG_RATIO_AC1_TEMP_neg'][:]
sig_ac_temp_neg_3_samba   = ac_temp_3_samba.variables['SIG_RATIO_AC1_TEMP_neg'][:]


#%% Ratio of lambda using single estimator
lambda_salt_samba     = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_34S.nc','r')
lambda_salt_2_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_34S_branch_1650_1500.nc','r')
lambda_salt_3_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_34S_branch_1650_600.nc','r')

lambda_temp_samba     = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_34S.nc','r')
lambda_temp_2_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_34S_branch_1650_1500.nc','r')
lambda_temp_3_samba   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_34S_branch_1650_600.nc','r')

ratio_lambda_salt_samba       = lambda_salt_samba.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_salt_2_samba     = lambda_salt_2_samba.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_salt_3_samba     = lambda_salt_3_samba.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_temp_samba       = lambda_temp_samba.variables['RATIO_LAMBDA_TEMP'][:]
ratio_lambda_temp_2_samba     = lambda_temp_2_samba.variables['RATIO_LAMBDA_TEMP'][:]
ratio_lambda_temp_3_samba     = lambda_temp_3_samba.variables['RATIO_LAMBDA_TEMP'][:]

sig_lambda_salt_pos_samba  = lambda_salt_samba.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_samba  = lambda_salt_samba.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_salt_pos_2_samba  = lambda_salt_2_samba.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_2_samba  = lambda_salt_2_samba.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_salt_pos_3_samba  = lambda_salt_3_samba.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_3_samba  = lambda_salt_3_samba.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_temp_pos_samba  = lambda_temp_samba.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_samba  = lambda_temp_samba.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

sig_lambda_temp_pos_2_samba  = lambda_temp_2_samba.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_2_samba  = lambda_temp_2_samba.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

sig_lambda_temp_pos_3_samba  = lambda_temp_3_samba.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_3_samba  = lambda_temp_3_samba.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]


#%% Difference use of branches

#SAMBA salinity (Figure 2 of paper)

cmap = colors.ListedColormap(['white'])

divnorm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=6)

fig, axs = plt.subplots(3, 3, figsize=(14, 12))
axs[0,0].contourf(lon_samba, depth_samba, ratio_var_salt_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[0,0].contourf(lon_samba, depth_samba, p_var_salt_sig_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,0].set_ylabel('Depth [m]', fontsize = 14)
axs[0,0].set_title('a) $R^S_{\mathrm{VAR}}$ (E$_2$/E$_1$)', fontsize = 16)
axs[0,0].set_ylim(depth_samba[-1], 0)
axs[0,0].set_xlim(-60,20)
axs[0,0].set_xticks([-60, -40, -20, 0, 20])
axs[0,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0,0].tick_params(axis='both', which='major', labelsize=12)
axs[0,0].set_facecolor('lightgrey')
axs[0,0].grid()

axs[0,1].contourf(lon_samba, depth_samba, ratio_var_salt_samba, cmap=cmap)
axs[0,1].contourf(lon_samba, depth_samba, ratio_ac_salt_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[0,1].contourf(lon_samba, depth_samba, sig_ac_salt_pos_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].contourf(lon_samba, depth_samba, sig_ac_salt_neg_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].set_title('b) $R^S_{\mathrm{AC}}$ (E$_2$/E$_1$)', fontsize = 16)
axs[0,1].set_ylim(depth_samba[-1], 0)
axs[0,1].set_xlim(-60,20)
axs[0,1].set_xticks([-60, -40, -20, 0, 20])
axs[0,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0,1].tick_params(axis='both', which='major', labelsize=12)
axs[0,1].set_facecolor('lightgrey')
axs[0,1].grid()

CS = axs[0,2].contourf(lon_samba, depth_samba, ratio_lambda_salt_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[0,2].contourf(lon_samba, depth_samba, sig_lambda_salt_pos_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].contourf(lon_samba, depth_samba, sig_lambda_salt_neg_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].set_title('c) $R^S_{\mathrm{RES}}$ (E$_2$/E$_1$)', fontsize = 16)
axs[0,2].set_ylim(depth_samba[-1], 0)
axs[0,2].set_xlim(-60,20)
axs[0,2].set_xticks([-60, -40, -20, 0, 20])
axs[0,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0,2].tick_params(axis='both', which='major', labelsize=12)
axs[0,2].set_facecolor('lightgrey')
axs[0,2].grid()
fig.colorbar(CS)

#fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[1,0].contourf(lon_samba, depth_samba, ratio_var_salt_3_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[1,0].contourf(lon_samba, depth_samba, p_var_salt_sig_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,0].set_ylabel('Depth [m]', fontsize = 14)
axs[1,0].set_title('d) $R^S_{\mathrm{VAR}}$ (E$_3$/E$_1$)', fontsize = 16)
axs[1,0].set_ylim(depth_samba[-1], 0)
axs[1,0].set_xlim(-60,20)
axs[1,0].set_xticks([-60, -40, -20, 0, 20])
axs[1,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1,0].tick_params(axis='both', which='major', labelsize=12)
axs[1,0].set_facecolor('lightgrey')
axs[1,0].grid()

axs[1,1].contourf(lon_samba, depth_samba, ratio_var_salt_3_samba, cmap=cmap)
axs[1,1].contourf(lon_samba, depth_samba, ratio_ac_salt_3_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[1,1].contourf(lon_samba, depth_samba, sig_ac_salt_pos_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].contourf(lon_samba, depth_samba, sig_ac_salt_neg_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].set_title('e) $R^S_{\mathrm{AC}}$ (E$_3$/E$_1$)', fontsize = 16)
axs[1,1].set_ylim(depth_samba[-1], 0)
axs[1,1].set_xlim(-60,20)
axs[1,1].set_xticks([-60, -40, -20, 0, 20])
axs[1,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1,1].tick_params(axis='both', which='major', labelsize=12)
axs[1,1].set_facecolor('lightgrey')
axs[1,1].grid()

CS = axs[1,2].contourf(lon_samba, depth_samba, ratio_lambda_salt_3_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[1,2].contourf(lon_samba, depth_samba, sig_lambda_salt_pos_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].contourf(lon_samba, depth_samba, sig_lambda_salt_neg_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].set_title('f) $R^S_{\mathrm{RES}}$ (E$_3$/E$_1$)', fontsize = 16)
axs[1,2].set_ylim(depth_samba[-1], 0)
axs[1,2].set_xlim(-60,20)
axs[1,2].set_xticks([-60, -40, -20, 0, 20])
axs[1,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1,2].tick_params(axis='both', which='major', labelsize=12)
axs[1,2].set_facecolor('lightgrey')
axs[1,2].grid()
fig.colorbar(CS)

#fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[2,0].contourf(lon_samba, depth_samba, ratio_var_salt_2_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[2,0].contourf(lon_samba, depth_samba, p_var_salt_sig_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,0].set_ylabel('Depth [m]', fontsize = 14)
axs[2,0].set_title('g) $R^S_{\mathrm{VAR}}$ (E$_3$/E$_2$)', fontsize = 16)
axs[2,0].set_ylim(depth_samba[-1], 0)
axs[2,0].set_xlim(-60,20)
axs[2,0].set_xticks([-60, -40, -20, 0, 20])
axs[2,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[2,0].tick_params(axis='both', which='major', labelsize=12)
axs[2,0].set_facecolor('lightgrey')
axs[2,0].grid()

axs[2,1].contourf(lon_samba, depth_samba, ratio_var_salt_2_samba, cmap=cmap)
axs[2,1].contourf(lon_samba, depth_samba, ratio_ac_salt_2_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[2,1].contourf(lon_samba, depth_samba, sig_ac_salt_pos_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].contourf(lon_samba, depth_samba, sig_ac_salt_neg_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].set_title('h) $R^S_{\mathrm{AC}}$ (E$_3$/E$_2$)', fontsize = 16)
axs[2,1].set_ylim(depth_samba[-1], 0)
axs[2,1].set_xlim(-60,20)
axs[2,1].set_xticks([-60, -40, -20, 0, 20])
axs[2,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[2,1].tick_params(axis='both', which='major', labelsize=12)
axs[2,1].set_facecolor('lightgrey')
axs[2,1].grid()

CS = axs[2,2].contourf(lon_samba, depth_samba, ratio_lambda_salt_2_samba, cmap='seismic', extend='both', norm = divnorm, levels=np.arange(0,5,0.25))
axs[2,2].contourf(lon_samba, depth_samba, sig_lambda_salt_pos_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].contourf(lon_samba, depth_samba, sig_lambda_salt_neg_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].set_title('i) $R^S_{\mathrm{RES}}$ (E$_3$/E$_2$)', fontsize = 16)
axs[2,2].set_ylim(depth_samba[-1], 0)
axs[2,2].set_xlim(-60,20)
axs[2,2].set_xticks([-60, -40, -20, 0, 20])
axs[2,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[2,2].tick_params(axis='both', which='major', labelsize=12)
axs[2,2].set_facecolor('lightgrey')
axs[2,2].grid()
fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'mean_ratio_VAR_AC1_LAMBDA_salt_branch12vsbranch23_SAMBA_lineardetrending_single_estimator.png', dpi=400)
plt.show()

#%%

# SAMBA temperature (Figure S1)
fig, axs = plt.subplots(3, 3, figsize=(14, 12))

axs[0,0].contourf(lon_samba, depth_samba, ratio_var_temp_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,0].contourf(lon_samba, depth_samba, p_var_temp_sig_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
axs[0,0].set_title('a) $R^T_{\mathrm{VAR}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,0].set_ylim(depth_samba[-1], 0)
axs[0,0].set_xlim(-60, 20)
axs[0,0].set_xticks([-60, -40, -20, 0, 20])
axs[0,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0,0].tick_params(axis='both', which='major', labelsize=12)
axs[0,0].set_facecolor('lightgrey')
axs[0,0].grid()

axs[0,1].contourf(lon_samba, depth_samba, ratio_var_temp_samba, cmap=cmap)
axs[0,1].contourf(lon_samba, depth_samba, ratio_ac_temp_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,1].contourf(lon_samba, depth_samba, sig_ac_temp_pos_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].contourf(lon_samba, depth_samba, sig_ac_temp_neg_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].set_title('b) $R^T_{\mathrm{AC}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,1].set_ylim(depth_samba[-1], 0)
axs[0,1].set_xlim(-60, 20)
axs[0,1].set_xticks([-60, -40, -20, 0, 20])
axs[0,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0,1].tick_params(axis='both', which='major', labelsize=12)
axs[0,1].set_facecolor('lightgrey')
axs[0,1].grid()

CS = axs[0,2].contourf(lon_samba, depth_samba, ratio_lambda_temp_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,2].contourf(lon_samba, depth_samba, sig_lambda_temp_pos_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].contourf(lon_samba, depth_samba, sig_lambda_temp_neg_samba, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].set_title('c) $R^T_{\mathrm{RES}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,2].set_ylim(depth_samba[-1], 0)
axs[0,2].set_xlim(-60, 20)
axs[0,2].set_xticks([-60, -40, -20, 0, 20])
axs[0,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0,2].tick_params(axis='both', which='major', labelsize=12)
axs[0,2].set_facecolor('lightgrey')
axs[0,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[1,0].contourf(lon_samba, depth_samba, ratio_var_temp_3_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,0].contourf(lon_samba, depth_samba, p_var_temp_sig_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_title('d) $R^T_{\mathrm{VAR}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,0].set_ylim(depth_samba[-1], 0)
axs[1,0].set_xlim(-60, 20)
axs[1,0].set_xticks([-60, -40, -20, 0, 20])
axs[1,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1,0].tick_params(axis='both', which='major', labelsize=12)
axs[1,0].set_facecolor('lightgrey')
axs[1,0].grid()

axs[1,1].contourf(lon_samba, depth_samba, ratio_var_temp_3_samba, cmap=cmap)
axs[1,1].contourf(lon_samba, depth_samba, ratio_ac_temp_3_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,1].contourf(lon_samba, depth_samba, sig_ac_temp_pos_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].contourf(lon_samba, depth_samba, sig_ac_temp_neg_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].set_title('e) $R^T_{\mathrm{AC}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,1].set_ylim(depth_samba[-1], 0)
axs[1,1].set_xlim(-60, 20)
axs[1,1].set_xticks([-60, -40, -20, 0, 20])
axs[1,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1,1].tick_params(axis='both', which='major', labelsize=12)
axs[1,1].set_facecolor('lightgrey')
axs[1,1].grid()

CS = axs[1,2].contourf(lon_samba, depth_samba, ratio_lambda_temp_3_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,2].contourf(lon_samba, depth_samba, sig_lambda_temp_pos_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].contourf(lon_samba, depth_samba, sig_lambda_temp_neg_3_samba, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].set_title('f) $R^T_{\mathrm{RES}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,2].set_ylim(depth_samba[-1], 0)
axs[1,2].set_xlim(-60, 20)
axs[1,2].set_xticks([-60, -40, -20, 0, 20])
axs[1,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1,2].tick_params(axis='both', which='major', labelsize=12)
axs[1,2].set_facecolor('lightgrey')
axs[1,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[2,0].contourf(lon_samba, depth_samba, ratio_var_temp_2_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,0].contourf(lon_samba, depth_samba, p_var_temp_sig_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,0].set_ylabel('Depth [m]', fontsize=14)
axs[2,0].set_title('g) $R^T_{\mathrm{VAR}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,0].set_ylim(depth_samba[-1], 0)
axs[2,0].set_xlim(-60, 20)
axs[2,0].set_xticks([-60, -40, -20, 0, 20])
axs[2,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[2,0].tick_params(axis='both', which='major', labelsize=12)
axs[2,0].set_facecolor('lightgrey')
axs[2,0].grid()

axs[2,1].contourf(lon_samba, depth_samba, ratio_var_temp_2_samba, cmap=cmap)
axs[2,1].contourf(lon_samba, depth_samba, ratio_ac_temp_2_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,1].contourf(lon_samba, depth_samba, sig_ac_temp_pos_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].contourf(lon_samba, depth_samba, sig_ac_temp_neg_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].set_title('h) $R^T_{\mathrm{AC}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,1].set_ylim(depth_samba[-1], 0)
axs[2,1].set_xlim(-60, 20)
axs[2,1].set_xticks([-60, -40, -20, 0, 20])
axs[2,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[2,1].tick_params(axis='both', which='major', labelsize=12)
axs[2,1].set_facecolor('lightgrey')
axs[2,1].grid()

CS = axs[2,2].contourf(lon_samba, depth_samba, ratio_lambda_temp_2_samba, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,2].contourf(lon_samba, depth_samba, sig_lambda_temp_pos_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].contourf(lon_samba, depth_samba, sig_lambda_temp_neg_2_samba, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].set_title('i) $R^T_{\mathrm{RES}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,2].set_ylim(depth_samba[-1], 0)
axs[2,2].set_xlim(-60, 20)
axs[2,2].set_xticks([-60, -40, -20, 0, 20])
axs[2,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[2,2].tick_params(axis='both', which='major', labelsize=12)
axs[2,2].set_facecolor('lightgrey')
axs[2,2].grid()
fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'mean_ratio_VAR_AC1_LAMBDA_temp_branch12vsbranch23_SAMBA_lineardetrending_single_estimator.png', dpi=400)
plt.show()


#%% RAPID

var_salt_rapid = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_26N.nc','r')
var_temp_rapid = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_26N.nc','r')
       
ratio_var_salt_rapid  = var_salt_rapid.variables['SIGMA_SALT'][:]
ratio_var_temp_rapid  = var_temp_rapid.variables['SIGMA_TEMP'][:]

p_var_salt_rapid  = var_salt_rapid.variables['p_value_LEVENE'][:]
p_var_temp_rapid  = var_temp_rapid.variables['p_value_LEVENE'][:]

lat_rapid = var_salt_rapid.variables['lat'][:]
lon_rapid = var_salt_rapid.variables['lon'][:]
depth_rapid = var_salt_rapid.variables['depth'][:]

p_var_salt_sig_rapid = ma.masked_all((len(depth_rapid), len(lon_rapid)))
p_var_temp_sig_rapid = ma.masked_all((len(depth_rapid), len(lon_rapid)))

# Create mask where the p-value is smaller than 0.05
for i in range(len(depth_rapid)):
    for j in range(len(lon_rapid)):
        if p_var_salt_rapid[i,j] < 0.05:# and ratio_var_salt_rapid[i,j] > 1:
            p_var_salt_sig_rapid[i,j] = 1
        else:
            continue
        
        if p_var_temp_rapid[i,j] < 0.05:# and ratio_var_temp_rapid[i,j] > 1:
            p_var_temp_sig_rapid[i,j] = 1
        else:
            continue

var_salt_2_rapid = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_26N_branch_1650_1500.nc','r')
var_temp_2_rapid = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_26N_branch_1650_1500.nc','r')
      
ratio_var_salt_2_rapid  = var_salt_2_rapid.variables['SIGMA_SALT'][:]
ratio_var_temp_2_rapid  = var_temp_2_rapid.variables['SIGMA_TEMP'][:]

p_var_salt_2_rapid  = var_salt_2_rapid.variables['p_value_LEVENE'][:]
p_var_temp_2_rapid  = var_temp_2_rapid.variables['p_value_LEVENE'][:]

lat_2_rapid = var_salt_2_rapid.variables['lat'][:]
lon_2_rapid = var_salt_2_rapid.variables['lon'][:]
depth_2_rapid = var_salt_2_rapid.variables['depth'][:]

p_var_salt_sig_2_rapid = ma.masked_all((len(depth_2_rapid), len(lon_2_rapid)))
p_var_temp_sig_2_rapid = ma.masked_all((len(depth_2_rapid), len(lon_2_rapid)))

# Create mask where the p-value is smaller than 0.05
for i in range(len(depth_2_rapid)):
    for j in range(len(lon_2_rapid)):
        if p_var_salt_2_rapid[i,j] < 0.05:# and ratio_var_salt_2_rapid[i,j] > 1:
            p_var_salt_sig_2_rapid[i,j] = 1
        else:
            continue
        
        if p_var_temp_2_rapid[i,j] < 0.05:# and ratio_var_temp_2_rapid[i,j] > 1:
            p_var_temp_sig_2_rapid[i,j] = 1
        else:
            continue

var_salt_3_rapid = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_26N_branch_1650_600.nc','r')
var_temp_3_rapid = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_26N_branch_1650_600.nc','r')
      
ratio_var_salt_3_rapid  = var_salt_3_rapid.variables['SIGMA_SALT'][:]
ratio_var_temp_3_rapid  = var_temp_3_rapid.variables['SIGMA_TEMP'][:]

p_var_salt_3_rapid  = var_salt_3_rapid.variables['p_value_LEVENE'][:]
p_var_temp_3_rapid  = var_temp_3_rapid.variables['p_value_LEVENE'][:]

lat_3_rapid = var_salt_3_rapid.variables['lat'][:]
lon_3_rapid = var_salt_3_rapid.variables['lon'][:]
depth_3_rapid = var_salt_3_rapid.variables['depth'][:]

p_var_salt_sig_3_rapid = ma.masked_all((len(depth_3_rapid), len(lon_3_rapid)))
p_var_temp_sig_3_rapid = ma.masked_all((len(depth_3_rapid), len(lon_3_rapid)))

# Create mask where the p-value is smaller than 0.05
for i in range(len(depth_3_rapid)):
    for j in range(len(lon_3_rapid)):
        if p_var_salt_3_rapid[i,j] < 0.05:# and ratio_var_salt_3_rapid[i,j] > 1:
            p_var_salt_sig_3_rapid[i,j] = 1
        else:
            continue
        
        if p_var_temp_3_rapid[i,j] < 0.05:# and ratio_var_temp_3_rapid[i,j] > 1:
            p_var_temp_sig_3_rapid[i,j] = 1
        else:
            continue

ac_salt_rapid     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_26N_branch_1500_600.nc','r')
ac_salt_2_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_26N_branch_1650_1500.nc','r')
ac_salt_3_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_26N_branch_1650_600.nc','r')

ac_temp_rapid     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_26N_branch_1500_600.nc','r')
ac_temp_2_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_26N_branch_1650_1500.nc','r')
ac_temp_3_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_26N_branch_1650_600.nc','r')

ratio_ac_salt_rapid       = ac_salt_rapid.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_2_rapid     = ac_salt_2_rapid.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_3_rapid     = ac_salt_3_rapid.variables['RATIO_AC1_SALT'][:]

ratio_ac_temp_rapid       = ac_temp_rapid.variables['RATIO_AC1_TEMP'][:]
ratio_ac_temp_2_rapid     = ac_temp_2_rapid.variables['RATIO_AC1_TEMP'][:]
ratio_ac_temp_3_rapid     = ac_temp_3_rapid.variables['RATIO_AC1_TEMP'][:]

sig_ac_salt_pos_rapid     = ac_salt_rapid.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_2_rapid   = ac_salt_2_rapid.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_3_rapid   = ac_salt_3_rapid.variables['SIG_RATIO_AC1_SALT_pos'][:]

sig_ac_salt_neg_rapid     = ac_salt_rapid.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_2_rapid   = ac_salt_2_rapid.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_3_rapid   = ac_salt_3_rapid.variables['SIG_RATIO_AC1_SALT_neg'][:]

sig_ac_temp_pos_rapid     = ac_temp_rapid.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_pos_2_rapid   = ac_temp_2_rapid.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_pos_3_rapid   = ac_temp_3_rapid.variables['SIG_RATIO_AC1_TEMP_pos'][:]

sig_ac_temp_neg_rapid     = ac_temp_rapid.variables['SIG_RATIO_AC1_TEMP_neg'][:]
sig_ac_temp_neg_2_rapid   = ac_temp_2_rapid.variables['SIG_RATIO_AC1_TEMP_neg'][:]
sig_ac_temp_neg_3_rapid   = ac_temp_3_rapid.variables['SIG_RATIO_AC1_TEMP_neg'][:]

lambda_salt_rapid     = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_26N.nc','r')
lambda_salt_2_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_26N_branch_1650_1500.nc','r')
lambda_salt_3_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_26N_branch_1650_600.nc','r')

lambda_temp_rapid     = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_26N.nc','r')
lambda_temp_2_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_26N_branch_1650_1500.nc','r')
lambda_temp_3_rapid   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_26N_branch_1650_600.nc','r')

ratio_lambda_salt_rapid       = lambda_salt_rapid.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_salt_2_rapid     = lambda_salt_2_rapid.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_salt_3_rapid     = lambda_salt_3_rapid.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_temp_rapid       = lambda_temp_rapid.variables['RATIO_LAMBDA_TEMP'][:]
ratio_lambda_temp_2_rapid     = lambda_temp_2_rapid.variables['RATIO_LAMBDA_TEMP'][:]
ratio_lambda_temp_3_rapid     = lambda_temp_3_rapid.variables['RATIO_LAMBDA_TEMP'][:]

sig_lambda_salt_pos_rapid  = lambda_salt_rapid.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_rapid  = lambda_salt_rapid.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_salt_pos_2_rapid  = lambda_salt_2_rapid.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_2_rapid  = lambda_salt_2_rapid.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_salt_pos_3_rapid  = lambda_salt_3_rapid.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_3_rapid  = lambda_salt_3_rapid.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_temp_pos_rapid  = lambda_temp_rapid.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_rapid  = lambda_temp_rapid.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

sig_lambda_temp_pos_2_rapid  = lambda_temp_2_rapid.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_2_rapid  = lambda_temp_2_rapid.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

sig_lambda_temp_pos_3_rapid  = lambda_temp_3_rapid.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_3_rapid  = lambda_temp_3_rapid.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

ratio_ac_temp_rapid = ma.masked_equal(ratio_ac_temp_rapid, 1.0)
ratio_var_temp_rapid = ma.masked_equal(ratio_var_temp_rapid, 1.0)
ratio_ac_salt_rapid = ma.masked_equal(ratio_ac_salt_rapid, 1.0)
ratio_var_salt_rapid = ma.masked_equal(ratio_var_salt_rapid, 1.0)


#%% Difference use of branches

# RAPID salinity
fig, axs = plt.subplots(3, 3, figsize=(14, 12))

axs[0,0].contourf(lon_rapid, depth_rapid, ratio_var_salt_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,0].contourf(lon_rapid, depth_rapid, p_var_salt_sig_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
axs[0,0].set_title('a) $R^S_{\mathrm{VAR}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,0].set_ylim(depth_rapid[-1], 0)
axs[0,0].set_xlim(-80, -20)
axs[0,0].set_xticks([-80, -60, -40, -20])
axs[0,0].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[0,0].tick_params(axis='both', which='major', labelsize=12)
axs[0,0].set_facecolor('lightgrey')
axs[0,0].grid()

axs[0,1].contourf(lon_rapid, depth_rapid, ratio_var_salt_rapid, cmap=cmap)
axs[0,1].contourf(lon_2_rapid, depth_rapid, ratio_ac_salt_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,1].contourf(lon_2_rapid, depth_rapid, sig_ac_salt_pos_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].contourf(lon_2_rapid, depth_rapid, sig_ac_salt_neg_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].set_title('b) $R^S_{\mathrm{AC}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,1].set_ylim(depth_rapid[-1], 0)
axs[0,1].set_xlim(-80, -20)
axs[0,1].set_xticks([-80, -60, -40, -20])
axs[0,1].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[0,1].tick_params(axis='both', which='major', labelsize=12)
axs[0,1].set_facecolor('lightgrey')
axs[0,1].grid()

CS = axs[0,2].contourf(lon_rapid, depth_rapid, ratio_lambda_salt_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,2].contourf(lon_rapid, depth_rapid, sig_lambda_salt_pos_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].contourf(lon_rapid, depth_rapid, sig_lambda_salt_neg_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].set_title('c) $R^S_{\mathrm{RES}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,2].set_ylim(depth_rapid[-1], 0)
axs[0,2].set_xlim(-80, -20)
axs[0,2].set_xticks([-80, -60, -40, -20])
axs[0,2].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[0,2].tick_params(axis='both', which='major', labelsize=12)
axs[0,2].set_facecolor('lightgrey')
axs[0,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[1,0].contourf(lon_2_rapid, depth_2_rapid, ratio_var_salt_3_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,0].contourf(lon_2_rapid, depth_2_rapid, p_var_salt_sig_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_title('d) $R^S_{\mathrm{VAR}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,0].set_ylim(depth_rapid[-1], 0)
axs[1,0].set_xlim(-80, -20)
axs[1,0].set_xticks([-80, -60, -40, -20])
axs[1,0].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[1,0].tick_params(axis='both', which='major', labelsize=12)
axs[1,0].set_facecolor('lightgrey')
axs[1,0].grid()

axs[1,1].contourf(lon_2_rapid, depth_rapid, ratio_var_salt_3_rapid, cmap=cmap)
axs[1,1].contourf(lon_2_rapid, depth_rapid, ratio_ac_salt_3_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,1].contourf(lon_2_rapid, depth_rapid, sig_ac_salt_pos_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].contourf(lon_2_rapid, depth_rapid, sig_ac_salt_neg_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].set_title('e) $R^S_{\mathrm{AC}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,1].set_ylim(depth_rapid[-1], 0)
axs[1,1].set_xlim(-80, -20)
axs[1,1].set_xticks([-80, -60, -40, -20])
axs[1,1].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[1,1].tick_params(axis='both', which='major', labelsize=12)
axs[1,1].set_facecolor('lightgrey')
axs[1,1].grid()

CS = axs[1,2].contourf(lon_2_rapid, depth_rapid, ratio_lambda_salt_3_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_salt_pos_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_salt_neg_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].set_title('f) $R^S_{\mathrm{RES}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,2].set_ylim(depth_rapid[-1], 0)
axs[1,2].set_xlim(-80, -20)
axs[1,2].set_xticks([-80, -60, -40, -20])
axs[1,2].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[1,2].tick_params(axis='both', which='major', labelsize=12)
axs[1,2].set_facecolor('lightgrey')
axs[1,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[2,0].contourf(lon_2_rapid, depth_rapid, ratio_var_salt_2_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,0].contourf(lon_2_rapid, depth_rapid, p_var_salt_sig_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,0].set_ylabel('Depth [m]', fontsize=14)
axs[2,0].set_title('g) $R^S_{\mathrm{VAR}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,0].set_ylim(depth_rapid[-1], 0)
axs[2,0].set_xlim(-80, -20)
axs[2,0].set_xticks([-80, -60, -40, -20])
axs[2,0].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[2,0].tick_params(axis='both', which='major', labelsize=12)
axs[2,0].set_facecolor('lightgrey')
axs[2,0].grid()

axs[2,1].contourf(lon_2_rapid, depth_rapid, ratio_var_salt_2_rapid, cmap=cmap)
axs[2,1].contourf(lon_2_rapid, depth_rapid, ratio_ac_salt_2_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,1].contourf(lon_2_rapid, depth_rapid, sig_ac_salt_pos_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].contourf(lon_2_rapid, depth_rapid, sig_ac_salt_neg_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].set_title('h) $R^S_{\mathrm{AC}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,1].set_ylim(depth_rapid[-1], 0)
axs[2,1].set_xlim(-80, -20)
axs[2,1].set_xticks([-80, -60, -40, -20])
axs[2,1].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[2,1].tick_params(axis='both', which='major', labelsize=12)
axs[2,1].set_facecolor('lightgrey')
axs[2,1].grid()

CS = axs[2,2].contourf(lon_2_rapid, depth_rapid, ratio_lambda_salt_2_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_salt_pos_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_salt_neg_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].set_title('i) $R^S_{\mathrm{RES}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,2].set_ylim(depth_rapid[-1], 0)
axs[2,2].set_xlim(-80, -20)
axs[2,2].set_xticks([-80, -60, -40, -20])
axs[2,2].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[2,2].tick_params(axis='both', which='major', labelsize=12)
axs[2,2].set_facecolor('lightgrey')
axs[2,2].grid()
fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'mean_ratio_VAR_AC1_LAMBDA_salt_branch12vsbranch23_RAPID_lineardetrending_single_estimator.png', dpi=400)
plt.show()

#%%

# RAPID temperature
fig, axs = plt.subplots(3, 3, figsize=(14, 12))

axs[0,0].contourf(lon_rapid, depth_rapid, ratio_var_temp_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,0].contourf(lon_rapid, depth_rapid, p_var_temp_sig_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
axs[0,0].set_title('a) $R^T_{\mathrm{VAR}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,0].set_ylim(depth_rapid[-1], 0)
axs[0,0].set_xlim(-80, -20)
axs[0,0].set_xticks([-80, -60, -40, -20])
axs[0,0].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[0,0].tick_params(axis='both', which='major', labelsize=12)
axs[0,0].set_facecolor('lightgrey')
axs[0,0].grid()

axs[0,1].contourf(lon_rapid, depth_rapid, ratio_var_temp_rapid, cmap=cmap)
axs[0,1].contourf(lon_2_rapid, depth_rapid, ratio_ac_temp_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,1].contourf(lon_2_rapid, depth_rapid, sig_ac_temp_pos_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].contourf(lon_2_rapid, depth_rapid, sig_ac_temp_neg_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].set_title('b) $R^T_{\mathrm{AC}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,1].set_ylim(depth_rapid[-1], 0)
axs[0,1].set_xlim(-80, -20)
axs[0,1].set_xticks([-80, -60, -40, -20])
axs[0,1].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[0,1].tick_params(axis='both', which='major', labelsize=12)
axs[0,1].set_facecolor('lightgrey')
axs[0,1].grid()

CS = axs[0,2].contourf(lon_rapid, depth_rapid, ratio_lambda_temp_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,2].contourf(lon_rapid, depth_rapid, sig_lambda_temp_pos_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].contourf(lon_rapid, depth_rapid, sig_lambda_temp_neg_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].set_title('c) $R^T_{\mathrm{RES}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,2].set_ylim(depth_rapid[-1], 0)
axs[0,2].set_xlim(-80, -20)
axs[0,2].set_xticks([-80, -60, -40, -20])
axs[0,2].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[0,2].tick_params(axis='both', which='major', labelsize=12)
axs[0,2].set_facecolor('lightgrey')
axs[0,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[1,0].contourf(lon_2_rapid, depth_rapid, ratio_var_temp_3_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,0].contourf(lon_2_rapid, depth_rapid, p_var_temp_sig_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_title('d) $R^T_{\mathrm{VAR}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,0].set_ylim(depth_rapid[-1], 0)
axs[1,0].set_xlim(-80, -20)
axs[1,0].set_xticks([-80, -60, -40, -20])
axs[1,0].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[1,0].tick_params(axis='both', which='major', labelsize=12)
axs[1,0].set_facecolor('lightgrey')
axs[1,0].grid()

axs[1,1].contourf(lon_2_rapid, depth_rapid, ratio_var_temp_3_rapid, cmap=cmap)
axs[1,1].contourf(lon_2_rapid, depth_rapid, ratio_ac_temp_3_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,1].contourf(lon_2_rapid, depth_rapid, sig_ac_temp_pos_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].contourf(lon_2_rapid, depth_rapid, sig_ac_temp_neg_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].set_title('e) $R^T_{\mathrm{AC}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,1].set_ylim(depth_rapid[-1], 0)
axs[1,1].set_xlim(-80, -20)
axs[1,1].set_xticks([-80, -60, -40, -20])
axs[1,1].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[1,1].tick_params(axis='both', which='major', labelsize=12)
axs[1,1].set_facecolor('lightgrey')
axs[1,1].grid()

CS = axs[1,2].contourf(lon_2_rapid, depth_rapid, ratio_lambda_temp_3_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_temp_pos_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_temp_neg_3_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].set_title('f) $R^T_{\mathrm{RES}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,2].set_ylim(depth_rapid[-1], 0)
axs[1,2].set_xlim(-80, -20)
axs[1,2].set_xticks([-80, -60, -40, -20])
axs[1,2].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[1,2].tick_params(axis='both', which='major', labelsize=12)
axs[1,2].set_facecolor('lightgrey')
axs[1,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[2,0].contourf(lon_2_rapid, depth_rapid, ratio_var_temp_2_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,0].contourf(lon_2_rapid, depth_rapid, p_var_temp_sig_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,0].set_ylabel('Depth [m]', fontsize=14)
axs[2,0].set_title('g) $R^T_{\mathrm{VAR}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,0].set_ylim(depth_rapid[-1], 0)
axs[2,0].set_xlim(-80, -20)
axs[2,0].set_xticks([-80, -60, -40, -20])
axs[2,0].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[2,0].tick_params(axis='both', which='major', labelsize=12)
axs[2,0].set_facecolor('lightgrey')
axs[2,0].grid()

axs[2,1].contourf(lon_2_rapid, depth_rapid, ratio_var_temp_2_rapid, cmap=cmap)
axs[2,1].contourf(lon_2_rapid, depth_rapid, ratio_ac_temp_2_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,1].contourf(lon_2_rapid, depth_rapid, sig_ac_temp_pos_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].contourf(lon_2_rapid, depth_rapid, sig_ac_temp_neg_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].set_title('h) $R^T_{\mathrm{AC}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,1].set_ylim(depth_rapid[-1], 0)
axs[2,1].set_xlim(-80, -20)
axs[2,1].set_xticks([-80, -60, -40, -20])
axs[2,1].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[2,1].tick_params(axis='both', which='major', labelsize=12)
axs[2,1].set_facecolor('lightgrey')
axs[2,1].grid()

CS = axs[2,2].contourf(lon_2_rapid, depth_rapid, ratio_lambda_temp_2_rapid, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_temp_pos_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].contourf(lon_2_rapid, depth_rapid, sig_lambda_temp_neg_2_rapid, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].set_title('i) $R^T_{\mathrm{RES}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,2].set_ylim(depth_rapid[-1], 0)
axs[2,2].set_xlim(-80, -20)
axs[2,2].set_xticks([-80, -60, -40, -20])
axs[2,2].set_xticklabels(['80$^{\circ}$W', '60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W'])
axs[2,2].tick_params(axis='both', which='major', labelsize=12)
axs[2,2].set_facecolor('lightgrey')
axs[2,2].grid()
fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'mean_ratio_VAR_AC1_LAMBDA_temp_branch12vsbranch23_RAPID_lineardetrending_single_estimator.png', dpi=400)

plt.show()

#%% OSNAP

var_salt_osnap = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_OSNAP.nc','r')
var_temp_osnap = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_OSNAP.nc','r')
       
ratio_var_salt_osnap  = var_salt_osnap.variables['SIGMA_SALT'][:]
ratio_var_temp_osnap  = var_temp_osnap.variables['SIGMA_TEMP'][:]

p_var_salt_osnap  = var_salt_osnap.variables['p_value_LEVENE'][:]
p_var_temp_osnap  = var_temp_osnap.variables['p_value_LEVENE'][:]

lat_osnap = var_salt_osnap.variables['lat'][:]
lon_osnap = var_salt_osnap.variables['lon'][:]
depth_osnap = var_salt_osnap.variables['depth'][:]

p_var_salt_sig_osnap = ma.masked_all((len(depth_osnap), len(lon_osnap)))
p_var_temp_sig_osnap = ma.masked_all((len(depth_osnap), len(lon_osnap)))

# Create mask where the p-value is smaller than 0.05
for i in range(len(depth_osnap)):
    for j in range(len(lon_osnap)):
        if p_var_salt_osnap[i,j] < 0.05:# and ratio_var_salt_osnap[i,j] > 1:
            p_var_salt_sig_osnap[i,j] = 1
        else:
            continue
        
        if p_var_temp_osnap[i,j] < 0.05:# and ratio_var_temp_osnap[i,j] > 1:
            p_var_temp_sig_osnap[i,j] = 1
        else:
            continue

var_salt_2_osnap = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_OSNAP_branch_1650_1500.nc','r')
var_temp_2_osnap = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_OSNAP_branch_1650_1500.nc','r')
      
ratio_var_salt_2_osnap  = var_salt_2_osnap.variables['SIGMA_SALT'][:]
ratio_var_temp_2_osnap  = var_temp_2_osnap.variables['SIGMA_TEMP'][:]

p_var_salt_2_osnap  = var_salt_2_osnap.variables['p_value_LEVENE'][:]
p_var_temp_2_osnap  = var_temp_2_osnap.variables['p_value_LEVENE'][:]

lat_2_osnap = var_salt_2_osnap.variables['lat'][:]
lon_2_osnap = var_salt_2_osnap.variables['lon'][:]
depth_2_osnap = var_salt_2_osnap.variables['depth'][:]

p_var_salt_sig_2_osnap = ma.masked_all((len(depth_2_osnap), len(lon_2_osnap)))
p_var_temp_sig_2_osnap = ma.masked_all((len(depth_2_osnap), len(lon_2_osnap)))

# Create mask where the p-value is smaller than 0.05
for i in range(len(depth_2_osnap)):
    for j in range(len(lon_2_osnap)):
        if p_var_salt_2_osnap[i,j] < 0.05:# and ratio_var_salt_2_osnap[i,j] > 1:
            p_var_salt_sig_2_osnap[i,j] = 1
        else:
            continue
        
        if p_var_temp_2_osnap[i,j] < 0.05:# and ratio_var_temp_2_osnap[i,j] > 1:
            p_var_temp_sig_2_osnap[i,j] = 1
        else:
            continue

var_salt_3_osnap = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_OSNAP_branch_1650_600.nc','r')
var_temp_3_osnap = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_OSNAP_branch_1650_600.nc','r')
      
ratio_var_salt_3_osnap  = var_salt_3_osnap.variables['SIGMA_SALT'][:]
ratio_var_temp_3_osnap  = var_temp_3_osnap.variables['SIGMA_TEMP'][:]

p_var_salt_3_osnap  = var_salt_3_osnap.variables['p_value_LEVENE'][:]
p_var_temp_3_osnap  = var_temp_3_osnap.variables['p_value_LEVENE'][:]

lat_3_osnap = var_salt_3_osnap.variables['lat'][:]
lon_3_osnap = var_salt_3_osnap.variables['lon'][:]
depth_3_osnap = var_salt_3_osnap.variables['depth'][:]

p_var_salt_sig_3_osnap = ma.masked_all((len(depth_3_osnap), len(lon_3_osnap)))
p_var_temp_sig_3_osnap = ma.masked_all((len(depth_3_osnap), len(lon_3_osnap)))

# Create mask where the p-value is smaller than 0.05
for i in range(len(depth_3_osnap)):
    for j in range(len(lon_3_osnap)):
        if p_var_salt_3_osnap[i,j] < 0.05:# and ratio_var_salt_3_osnap[i,j] > 1:
            p_var_salt_sig_3_osnap[i,j] = 1
        else:
            continue
        
        if p_var_temp_3_osnap[i,j] < 0.05:# and ratio_var_temp_3_osnap[i,j] > 1:
            p_var_temp_sig_3_osnap[i,j] = 1
        else:
            continue

ac_salt_osnap     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_OSNAP_branch_1500_600.nc','r')
ac_salt_2_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_OSNAP_branch_1650_1500.nc','r')
ac_salt_3_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_OSNAP_branch_1650_600.nc','r')

ac_temp_osnap     = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_OSNAP_branch_1500_600.nc','r')
ac_temp_2_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_OSNAP_branch_1650_1500.nc','r')
ac_temp_3_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_OSNAP_branch_1650_600.nc','r')

ratio_ac_salt_osnap       = ac_salt_osnap.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_2_osnap     = ac_salt_2_osnap.variables['RATIO_AC1_SALT'][:]
ratio_ac_salt_3_osnap     = ac_salt_3_osnap.variables['RATIO_AC1_SALT'][:]

ratio_ac_temp_osnap       = ac_temp_osnap.variables['RATIO_AC1_TEMP'][:]
ratio_ac_temp_2_osnap     = ac_temp_2_osnap.variables['RATIO_AC1_TEMP'][:]
ratio_ac_temp_3_osnap     = ac_temp_3_osnap.variables['RATIO_AC1_TEMP'][:]

sig_ac_salt_pos_osnap     = ac_salt_osnap.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_2_osnap   = ac_salt_2_osnap.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_pos_3_osnap   = ac_salt_3_osnap.variables['SIG_RATIO_AC1_SALT_pos'][:]

sig_ac_salt_neg_osnap     = ac_salt_osnap.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_2_osnap   = ac_salt_2_osnap.variables['SIG_RATIO_AC1_SALT_neg'][:]
sig_ac_salt_neg_3_osnap   = ac_salt_3_osnap.variables['SIG_RATIO_AC1_SALT_neg'][:]

sig_ac_temp_pos_osnap     = ac_temp_osnap.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_pos_2_osnap   = ac_temp_2_osnap.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_pos_3_osnap   = ac_temp_3_osnap.variables['SIG_RATIO_AC1_TEMP_pos'][:]

sig_ac_temp_neg_osnap     = ac_temp_osnap.variables['SIG_RATIO_AC1_TEMP_neg'][:]
sig_ac_temp_neg_2_osnap   = ac_temp_2_osnap.variables['SIG_RATIO_AC1_TEMP_neg'][:]
sig_ac_temp_neg_3_osnap   = ac_temp_3_osnap.variables['SIG_RATIO_AC1_TEMP_neg'][:]

lambda_salt_osnap     = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_OSNAP.nc','r')
lambda_salt_2_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_OSNAP_branch_1650_1500.nc','r')
lambda_salt_3_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_OSNAP_branch_1650_600.nc','r')

lambda_temp_osnap     = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_OSNAP.nc','r')
lambda_temp_2_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_OSNAP_branch_1650_1500.nc','r')
lambda_temp_3_osnap   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_OSNAP_branch_1650_600.nc','r')

ratio_lambda_salt_osnap       = lambda_salt_osnap.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_salt_2_osnap     = lambda_salt_2_osnap.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_salt_3_osnap     = lambda_salt_3_osnap.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_temp_osnap       = lambda_temp_osnap.variables['RATIO_LAMBDA_TEMP'][:]
ratio_lambda_temp_2_osnap     = lambda_temp_2_osnap.variables['RATIO_LAMBDA_TEMP'][:]
ratio_lambda_temp_3_osnap     = lambda_temp_3_osnap.variables['RATIO_LAMBDA_TEMP'][:]

sig_lambda_salt_pos_osnap  = lambda_salt_osnap.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_osnap  = lambda_salt_osnap.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_salt_pos_2_osnap  = lambda_salt_2_osnap.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_2_osnap  = lambda_salt_2_osnap.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_salt_pos_3_osnap  = lambda_salt_3_osnap.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg_3_osnap  = lambda_salt_3_osnap.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_temp_pos_osnap  = lambda_temp_osnap.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_osnap  = lambda_temp_osnap.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

sig_lambda_temp_pos_2_osnap  = lambda_temp_2_osnap.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_2_osnap  = lambda_temp_2_osnap.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

sig_lambda_temp_pos_3_osnap  = lambda_temp_3_osnap.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg_3_osnap  = lambda_temp_3_osnap.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

ratio_ac_temp_osnap = ma.masked_equal(ratio_ac_temp_osnap, 1.0)
ratio_var_temp_osnap = ma.masked_equal(ratio_var_temp_osnap, 1.0)
ratio_ac_salt_osnap = ma.masked_equal(ratio_ac_salt_osnap, 1.0)
ratio_var_salt_osnap = ma.masked_equal(ratio_var_salt_osnap, 1.0)

#%% Difference use of branches

# OSNAP salinity
fig, axs = plt.subplots(3, 3, figsize=(14, 12))

axs[0,0].contourf(lon_osnap, depth_osnap, ratio_var_salt_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,0].contourf(lon_osnap, depth_osnap, p_var_salt_sig_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
axs[0,0].set_title('a) $R^S_{\mathrm{VAR}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,0].set_ylim(depth_osnap[-1], 0)
axs[0,0].set_xlim(-60, 0)
axs[0,0].set_xticks([-60, -40, -20, 0])
axs[0,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[0,0].tick_params(axis='both', which='major', labelsize=12)
axs[0,0].set_facecolor('lightgrey')
axs[0,0].grid()

axs[0,1].contourf(lon_osnap, depth_osnap, ratio_var_salt_osnap, cmap=cmap)
axs[0,1].contourf(lon_osnap, depth_osnap, ratio_ac_salt_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,1].contourf(lon_osnap, depth_osnap, sig_ac_salt_pos_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].contourf(lon_osnap, depth_osnap, sig_ac_salt_neg_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].set_title('b) $R^S_{\mathrm{AC}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,1].set_ylim(depth_osnap[-1], 0)
axs[0,1].set_xlim(-60, 0)
axs[0,1].set_xticks([-60, -40, -20, 0])
axs[0,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[0,1].tick_params(axis='both', which='major', labelsize=12)
axs[0,1].set_facecolor('lightgrey')
axs[0,1].grid()

CS = axs[0,2].contourf(lon_osnap, depth_osnap, ratio_lambda_salt_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,2].contourf(lon_osnap, depth_osnap, sig_lambda_salt_pos_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].contourf(lon_osnap, depth_osnap, sig_lambda_salt_neg_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].set_title('c) $R^S_{\mathrm{RES}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,2].set_ylim(depth_osnap[-1], 0)
axs[0,2].set_xlim(-60, 0)
axs[0,2].set_xticks([-60, -40, -20, 0])
axs[0,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[0,2].tick_params(axis='both', which='major', labelsize=12)
axs[0,2].set_facecolor('lightgrey')
axs[0,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[1,0].contourf(lon_osnap, depth_osnap, ratio_var_salt_3_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,0].contourf(lon_osnap, depth_osnap, p_var_salt_sig_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_title('d) $R^S_{\mathrm{VAR}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,0].set_ylim(depth_osnap[-1], 0)
axs[1,0].set_xlim(-60, 0)
axs[1,0].set_xticks([-60, -40, -20, 0])
axs[1,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[1,0].tick_params(axis='both', which='major', labelsize=12)
axs[1,0].set_facecolor('lightgrey')
axs[1,0].grid()

axs[1,1].contourf(lon_osnap, depth_osnap, ratio_var_salt_3_osnap, cmap=cmap)
axs[1,1].contourf(lon_osnap, depth_osnap, ratio_ac_salt_3_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,1].contourf(lon_osnap, depth_osnap, sig_ac_salt_pos_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].contourf(lon_osnap, depth_osnap, sig_ac_salt_neg_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].set_title('e) $R^S_{\mathrm{AC}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,1].set_ylim(depth_osnap[-1], 0)
axs[1,1].set_xlim(-60, 0)
axs[1,1].set_xticks([-60, -40, -20, 0])
axs[1,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[1,1].tick_params(axis='both', which='major', labelsize=12)
axs[1,1].set_facecolor('lightgrey')
axs[1,1].grid()

CS = axs[1,2].contourf(lon_osnap, depth_osnap, ratio_lambda_salt_3_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,2].contourf(lon_osnap, depth_osnap, sig_lambda_salt_pos_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].contourf(lon_osnap, depth_osnap, sig_lambda_salt_neg_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].set_title('f) $R^S_{\mathrm{RES}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,2].set_ylim(depth_osnap[-1], 0)
axs[1,2].set_xlim(-60, 0)
axs[1,2].set_xticks([-60, -40, -20, 0])
axs[1,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[1,2].tick_params(axis='both', which='major', labelsize=12)
axs[1,2].set_facecolor('lightgrey')
axs[1,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[2,0].contourf(lon_osnap, depth_osnap, ratio_var_salt_2_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,0].contourf(lon_osnap, depth_osnap, p_var_salt_sig_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,0].set_ylabel('Depth [m]', fontsize=14)
axs[2,0].set_title('g) $R^S_{\mathrm{VAR}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,0].set_ylim(depth_osnap[-1], 0)
axs[2,0].set_xlim(-60, 0)
axs[2,0].set_xticks([-60, -40, -20, 0])
axs[2,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[2,0].tick_params(axis='both', which='major', labelsize=12)
axs[2,0].set_facecolor('lightgrey')
axs[2,0].grid()

axs[2,1].contourf(lon_osnap, depth_osnap, ratio_var_salt_2_osnap, cmap=cmap)
axs[2,1].contourf(lon_osnap, depth_osnap, ratio_ac_salt_2_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,1].contourf(lon_osnap, depth_osnap, sig_ac_salt_pos_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].contourf(lon_osnap, depth_osnap, sig_ac_salt_neg_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].set_title('h) $R^S_{\mathrm{AC}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,1].set_ylim(depth_osnap[-1], 0)
axs[2,1].set_xlim(-60, 0)
axs[2,1].set_xticks([-60, -40, -20, 0])
axs[2,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[2,1].tick_params(axis='both', which='major', labelsize=12)
axs[2,1].set_facecolor('lightgrey')
axs[2,1].grid()

CS = axs[2,2].contourf(lon_osnap, depth_osnap, ratio_lambda_salt_2_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,2].contourf(lon_osnap, depth_osnap, sig_lambda_salt_pos_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].contourf(lon_osnap, depth_osnap, sig_lambda_salt_neg_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].set_title('i) $R^S_{\mathrm{RES}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,2].set_ylim(depth_osnap[-1], 0)
axs[2,2].set_xlim(-60, 0)
axs[2,2].set_xticks([-60, -40, -20, 0])
axs[2,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[2,2].tick_params(axis='both', which='major', labelsize=12)
axs[2,2].set_facecolor('lightgrey')
axs[2,2].grid()
fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'mean_ratio_VAR_AC1_LAMBDA_salt_branch12vsbranch23_OSNAP_lineardetrending_single_estimator.png', dpi=400)

plt.show()

#%%

# OSNAP temperature
fig, axs = plt.subplots(3, 3, figsize=(14, 12))

axs[0,0].contourf(lon_osnap, depth_osnap, ratio_var_temp_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,0].contourf(lon_osnap, depth_osnap, p_var_temp_sig_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
axs[0,0].set_title('a) $R^T_{\mathrm{VAR}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,0].set_ylim(depth_osnap[-1], 0)
axs[0,0].set_xlim(-60, 0)
axs[0,0].set_xticks([-60, -40, -20, 0])
axs[0,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[0,0].tick_params(axis='both', which='major', labelsize=12)
axs[0,0].set_facecolor('lightgrey')
axs[0,0].grid()

axs[0,1].contourf(lon_osnap, depth_osnap, ratio_var_temp_osnap, cmap=cmap)
axs[0,1].contourf(lon_osnap, depth_osnap, ratio_ac_temp_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,1].contourf(lon_osnap, depth_osnap, sig_ac_temp_pos_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].contourf(lon_osnap, depth_osnap, sig_ac_temp_neg_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,1].set_title('b) $R^T_{\mathrm{AC}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,1].set_ylim(depth_osnap[-1], 0)
axs[0,1].set_xlim(-60, 0)
axs[0,1].set_xticks([-60, -40, -20, 0])
axs[0,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[0,1].tick_params(axis='both', which='major', labelsize=12)
axs[0,1].set_facecolor('lightgrey')
axs[0,1].grid()

CS = axs[0,2].contourf(lon_osnap, depth_osnap, ratio_lambda_temp_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[0,2].contourf(lon_osnap, depth_osnap, sig_lambda_temp_pos_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].contourf(lon_osnap, depth_osnap, sig_lambda_temp_neg_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[0,2].set_title('c) $R^T_{\mathrm{RES}}$ (E$_2$/E$_1$)', fontsize=16)
axs[0,2].set_ylim(depth_osnap[-1], 0)
axs[0,2].set_xlim(-60, 0)
axs[0,2].set_xticks([-60, -40, -20, 0])
axs[0,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[0,2].tick_params(axis='both', which='major', labelsize=12)
axs[0,2].set_facecolor('lightgrey')
axs[0,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[1,0].contourf(lon_osnap, depth_osnap, ratio_var_temp_3_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,0].contourf(lon_osnap, depth_osnap, p_var_temp_sig_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_title('d) $R^T_{\mathrm{VAR}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,0].set_ylim(depth_osnap[-1], 0)
axs[1,0].set_xlim(-60, 0)
axs[1,0].set_xticks([-60, -40, -20, 0])
axs[1,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[1,0].tick_params(axis='both', which='major', labelsize=12)
axs[1,0].set_facecolor('lightgrey')
axs[1,0].grid()

axs[1,1].contourf(lon_osnap, depth_osnap, ratio_var_temp_3_osnap, cmap=cmap)
axs[1,1].contourf(lon_osnap, depth_osnap, ratio_ac_temp_3_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,1].contourf(lon_osnap, depth_osnap, sig_ac_temp_pos_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].contourf(lon_osnap, depth_osnap, sig_ac_temp_neg_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,1].set_title('e) $R^T_{\mathrm{AC}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,1].set_ylim(depth_osnap[-1], 0)
axs[1,1].set_xlim(-60, 0)
axs[1,1].set_xticks([-60, -40, -20, 0])
axs[1,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[1,1].tick_params(axis='both', which='major', labelsize=12)
axs[1,1].set_facecolor('lightgrey')
axs[1,1].grid()

CS = axs[1,2].contourf(lon_osnap, depth_osnap, ratio_lambda_temp_3_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[1,2].contourf(lon_osnap, depth_osnap, sig_lambda_temp_pos_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].contourf(lon_osnap, depth_osnap, sig_lambda_temp_neg_3_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[1,2].set_title('f) $R^T_{\mathrm{RES}}$ (E$_3$/E$_1$)', fontsize=16)
axs[1,2].set_ylim(depth_osnap[-1], 0)
axs[1,2].set_xlim(-60, 0)
axs[1,2].set_xticks([-60, -40, -20, 0])
axs[1,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[1,2].tick_params(axis='both', which='major', labelsize=12)
axs[1,2].set_facecolor('lightgrey')
axs[1,2].grid()
fig.colorbar(CS)

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))
axs[2,0].contourf(lon_osnap, depth_osnap, ratio_var_temp_2_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,0].contourf(lon_osnap, depth_osnap, p_var_temp_sig_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,0].set_ylabel('Depth [m]', fontsize=14)
axs[2,0].set_title('g) $R^T_{\mathrm{VAR}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,0].set_ylim(depth_osnap[-1], 0)
axs[2,0].set_xlim(-60, 0)
axs[2,0].set_xticks([-60, -40, -20, 0])
axs[2,0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[2,0].tick_params(axis='both', which='major', labelsize=12)
axs[2,0].set_facecolor('lightgrey')
axs[2,0].grid()

axs[2,1].contourf(lon_osnap, depth_osnap, ratio_var_temp_2_osnap, cmap=cmap)
axs[2,1].contourf(lon_osnap, depth_osnap, ratio_ac_temp_2_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,1].contourf(lon_osnap, depth_osnap, sig_ac_temp_pos_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].contourf(lon_osnap, depth_osnap, sig_ac_temp_neg_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,1].set_title('h) $R^T_{\mathrm{AC}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,1].set_ylim(depth_osnap[-1], 0)
axs[2,1].set_xlim(-60, 0)
axs[2,1].set_xticks([-60, -40, -20, 0])
axs[2,1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[2,1].tick_params(axis='both', which='major', labelsize=12)
axs[2,1].set_facecolor('lightgrey')
axs[2,1].grid()

CS = axs[2,2].contourf(lon_osnap, depth_osnap, ratio_lambda_temp_2_osnap, cmap='seismic', extend='both', norm=divnorm, levels=np.arange(0,5,0.25))
axs[2,2].contourf(lon_osnap, depth_osnap, sig_lambda_temp_pos_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].contourf(lon_osnap, depth_osnap, sig_lambda_temp_neg_2_osnap, colors='none', linewidth=1.0, hatches='xx')
axs[2,2].set_title('i) $R^T_{\mathrm{RES}}$ (E$_3$/E$_2$)', fontsize=16)
axs[2,2].set_ylim(depth_osnap[-1], 0)
axs[2,2].set_xlim(-60, 0)
axs[2,2].set_xticks([-60, -40, -20, 0])
axs[2,2].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$'])
axs[2,2].tick_params(axis='both', which='major', labelsize=12)
axs[2,2].set_facecolor('lightgrey')
axs[2,2].grid()
fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'mean_ratio_VAR_AC1_LAMBDA_temp_branch12vsbranch23_OSNAP_lineardetrending_single_estimator.png', dpi=400)

plt.show()

#%% Determining how many significant gridcells there are per ratio

number_grid_cells_samba = ratio_var_salt_samba.count()
number_grid_cells_rapid = ratio_var_salt_rapid.count() #longitude must run from -60 to -20, but its the same when you take the whole array. its only bathymetry anyway
number_grid_cells_osnap = ratio_var_salt_osnap.count()

index_var_salt_samba = p_var_salt_sig_samba.count()/number_grid_cells_samba
index_var_temp_samba = p_var_temp_sig_samba.count()/number_grid_cells_samba

index_var_salt_samba_2 = p_var_salt_sig_2_samba.count()/number_grid_cells_samba
index_var_temp_samba_2 = p_var_temp_sig_2_samba.count()/number_grid_cells_samba

index_var_salt_samba_3 = p_var_salt_sig_3_samba.count()/number_grid_cells_samba
index_var_temp_samba_3 = p_var_temp_sig_3_samba.count()/number_grid_cells_samba

index_ac_salt_samba = sig_ac_salt_pos_samba.count()/number_grid_cells_samba
index_ac_temp_samba = sig_ac_temp_pos_samba.count()/number_grid_cells_samba

index_ac_salt_samba_2 = sig_ac_salt_pos_2_samba.count()/number_grid_cells_samba
index_ac_temp_samba_2 = sig_ac_temp_pos_2_samba.count()/number_grid_cells_samba

index_ac_salt_samba_3 = sig_ac_salt_pos_3_samba.count()/number_grid_cells_samba
index_ac_temp_samba_3 = sig_ac_temp_pos_3_samba.count()/number_grid_cells_samba

index_lambda_salt_samba = sig_lambda_salt_neg_samba.count()/number_grid_cells_samba
index_lambda_temp_samba = sig_lambda_temp_neg_samba.count()/number_grid_cells_samba

index_lambda_salt_samba_2 = sig_lambda_salt_neg_2_samba.count()/number_grid_cells_samba
index_lambda_temp_samba_2 = sig_lambda_temp_neg_2_samba.count()/number_grid_cells_samba

index_lambda_salt_samba_3 = sig_lambda_salt_neg_3_samba.count()/number_grid_cells_samba
index_lambda_temp_samba_3 = sig_lambda_temp_neg_3_samba.count()/number_grid_cells_samba

#%%

index_var_salt_rapid = p_var_salt_sig_rapid.count()/number_grid_cells_rapid
index_var_temp_rapid = p_var_temp_sig_rapid.count()/number_grid_cells_rapid

index_var_salt_rapid_2 = p_var_salt_sig_2_rapid.count()/number_grid_cells_rapid
index_var_temp_rapid_2 = p_var_temp_sig_2_rapid.count()/number_grid_cells_rapid

index_var_salt_rapid_3 = p_var_salt_sig_3_rapid.count()/number_grid_cells_rapid
index_var_temp_rapid_3 = p_var_temp_sig_3_rapid.count()/number_grid_cells_rapid

index_ac_salt_rapid = sig_ac_salt_pos_rapid.count()/number_grid_cells_rapid
index_ac_temp_rapid = sig_ac_temp_pos_rapid.count()/number_grid_cells_rapid

index_ac_salt_rapid_2 = sig_ac_salt_pos_2_rapid.count()/number_grid_cells_rapid
index_ac_temp_rapid_2 = sig_ac_temp_pos_2_rapid.count()/number_grid_cells_rapid

index_ac_salt_rapid_3 = sig_ac_salt_pos_3_rapid.count()/number_grid_cells_rapid
index_ac_temp_rapid_3 = sig_ac_temp_pos_3_rapid.count()/number_grid_cells_rapid

index_lambda_salt_rapid = sig_lambda_salt_neg_rapid.count()/number_grid_cells_rapid
index_lambda_temp_rapid = sig_lambda_temp_neg_rapid.count()/number_grid_cells_rapid

index_lambda_salt_rapid_2 = sig_lambda_salt_neg_2_rapid.count()/number_grid_cells_rapid
index_lambda_temp_rapid_2 = sig_lambda_temp_neg_2_rapid.count()/number_grid_cells_rapid

index_lambda_salt_rapid_3 = sig_lambda_salt_neg_3_rapid.count()/number_grid_cells_rapid
index_lambda_temp_rapid_3 = sig_lambda_temp_neg_3_rapid.count()/number_grid_cells_rapid

#%%

index_var_salt_osnap = p_var_salt_sig_osnap.count()/number_grid_cells_osnap
index_var_temp_osnap = p_var_temp_sig_osnap.count()/number_grid_cells_osnap

index_var_salt_osnap_2 = p_var_salt_sig_2_osnap.count()/number_grid_cells_osnap
index_var_temp_osnap_2 = p_var_temp_sig_2_osnap.count()/number_grid_cells_osnap

index_var_salt_osnap_3 = p_var_salt_sig_3_osnap.count()/number_grid_cells_osnap
index_var_temp_osnap_3 = p_var_temp_sig_3_osnap.count()/number_grid_cells_osnap

index_ac_salt_osnap = sig_ac_salt_pos_osnap.count()/number_grid_cells_osnap
index_ac_temp_osnap = sig_ac_temp_pos_osnap.count()/number_grid_cells_osnap

index_ac_salt_osnap_2 = sig_ac_salt_pos_2_osnap.count()/number_grid_cells_osnap
index_ac_temp_osnap_2 = sig_ac_temp_pos_2_osnap.count()/number_grid_cells_osnap

index_ac_salt_osnap_3 = sig_ac_salt_pos_3_osnap.count()/number_grid_cells_osnap
index_ac_temp_osnap_3 = sig_ac_temp_pos_3_osnap.count()/number_grid_cells_osnap

index_lambda_salt_osnap = sig_lambda_salt_neg_osnap.count()/number_grid_cells_osnap
index_lambda_temp_osnap = sig_lambda_temp_neg_osnap.count()/number_grid_cells_osnap

index_lambda_salt_osnap_2 = sig_lambda_salt_neg_2_osnap.count()/number_grid_cells_osnap
index_lambda_temp_osnap_2 = sig_lambda_temp_neg_2_osnap.count()/number_grid_cells_osnap

index_lambda_salt_osnap_3 = sig_lambda_salt_neg_3_osnap.count()/number_grid_cells_osnap
index_lambda_temp_osnap_3 = sig_lambda_temp_neg_3_osnap.count()/number_grid_cells_osnap


#%% Now normalize per area

#dz and dx of grid points along SAMBA
fh_grid_samba = netcdf.Dataset(directory + 'SAMBA_lon_lat_depth_area.nc','r')

depth_grid 	= fh_grid_samba.variables['depth'][:]
lon_grid	= fh_grid_samba.variables['lon'][:]
dx_samba = fh_grid_samba.variables['dx'][:]
dz_samba = fh_grid_samba.variables['dz'][:]

area_grids_samba = ma.masked_all((len(depth_grid), len(lon_grid)))

for depth_i in range(len(depth_grid)):
    for lon_i in range(len(lon_grid)):
        area_grids_samba[depth_i, lon_i] = dx_samba[lon_i]*dz_samba[depth_i]

if np.shape(p_var_salt_sig_samba) != np.shape(area_grids_samba):
    print('Shapes are not the same!!')
    
total_area_samba = np.nansum(area_grids_samba[~p_var_salt_samba.mask])

area_index_var_salt_samba = np.nansum(area_grids_samba[~p_var_salt_sig_samba.mask])/total_area_samba
area_index_var_temp_samba = np.nansum(area_grids_samba[~p_var_temp_sig_samba.mask])/total_area_samba

area_index_var_salt_samba_2 = np.nansum(area_grids_samba[~p_var_salt_sig_2_samba.mask])/total_area_samba
area_index_var_temp_samba_2 = np.nansum(area_grids_samba[~p_var_temp_sig_2_samba.mask])/total_area_samba

area_index_var_salt_samba_3 = np.nansum(area_grids_samba[~p_var_salt_sig_3_samba.mask])/total_area_samba
area_index_var_temp_samba_3 = np.nansum(area_grids_samba[~p_var_temp_sig_3_samba.mask])/total_area_samba

area_index_ac_salt_samba = np.nansum(area_grids_samba[~sig_ac_salt_pos_samba.mask])/total_area_samba
area_index_ac_temp_samba = np.nansum(area_grids_samba[~sig_ac_temp_pos_samba.mask])/total_area_samba

area_index_ac_salt_samba_2 = np.nansum(area_grids_samba[~sig_ac_salt_pos_2_samba.mask])/total_area_samba
area_index_ac_temp_samba_2 = np.nansum(area_grids_samba[~sig_ac_temp_pos_2_samba.mask])/total_area_samba

area_index_ac_salt_samba_3 = np.nansum(area_grids_samba[~sig_ac_salt_pos_3_samba.mask])/total_area_samba
area_index_ac_temp_samba_3 = np.nansum(area_grids_samba[~sig_ac_temp_pos_3_samba.mask])/total_area_samba

area_index_lambda_salt_samba = np.nansum(area_grids_samba[~sig_lambda_salt_neg_samba.mask])/total_area_samba
area_index_lambda_temp_samba = np.nansum(area_grids_samba[~sig_lambda_temp_neg_samba.mask])/total_area_samba

area_index_lambda_salt_samba_2 = np.nansum(area_grids_samba[~sig_lambda_salt_neg_2_samba.mask])/total_area_samba
area_index_lambda_temp_samba_2 = np.nansum(area_grids_samba[~sig_lambda_temp_neg_2_samba.mask])/total_area_samba

area_index_lambda_salt_samba_3 = np.nansum(area_grids_samba[~sig_lambda_salt_neg_3_samba.mask])/total_area_samba
area_index_lambda_temp_samba_3 = np.nansum(area_grids_samba[~sig_lambda_temp_neg_3_samba.mask])/total_area_samba

#%% 

#dz and dx of grid points along RAPID
fh_grid_rapid = netcdf.Dataset(directory + 'RAPID_lon_lat_depth_area.nc','r')

depth_grid 	= fh_grid_rapid.variables['depth'][:]
lon_grid	= fh_grid_rapid.variables['lon'][:]
dx_rapid = fh_grid_rapid.variables['dx'][:]
dz_rapid = fh_grid_rapid.variables['dz'][:]

area_grids_rapid = ma.masked_all((len(depth_grid), len(lon_grid)))

for depth_i in range(len(depth_grid)):
    for lon_i in range(len(lon_grid)):
        area_grids_rapid[depth_i, lon_i] = dx_rapid[lon_i]*dz_rapid[depth_i]

if np.shape(p_var_salt_sig_rapid[:,0:len(lon_grid)]) != np.shape(area_grids_rapid):
    print('Shapes are not the same!!')
    
total_area_rapid = np.nansum(area_grids_rapid[~ratio_ac_salt_rapid[:,0:len(lon_grid)].mask])

area_index_var_salt_rapid = np.nansum(area_grids_rapid[~p_var_salt_sig_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_var_temp_rapid = np.nansum(area_grids_rapid[~p_var_temp_sig_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_var_salt_rapid_2 = np.nansum(area_grids_rapid[~p_var_salt_sig_2_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_var_temp_rapid_2 = np.nansum(area_grids_rapid[~p_var_temp_sig_2_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_var_salt_rapid_3 = np.nansum(area_grids_rapid[~p_var_salt_sig_3_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_var_temp_rapid_3 = np.nansum(area_grids_rapid[~p_var_temp_sig_3_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_ac_salt_rapid = np.nansum(area_grids_rapid[~sig_ac_salt_pos_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_ac_temp_rapid = np.nansum(area_grids_rapid[~sig_ac_temp_pos_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_ac_salt_rapid_2 = np.nansum(area_grids_rapid[~sig_ac_salt_pos_2_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_ac_temp_rapid_2 = np.nansum(area_grids_rapid[~sig_ac_temp_pos_2_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_ac_salt_rapid_3 = np.nansum(area_grids_rapid[~sig_ac_salt_pos_3_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_ac_temp_rapid_3 = np.nansum(area_grids_rapid[~sig_ac_temp_pos_3_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_lambda_salt_rapid = np.nansum(area_grids_rapid[~sig_lambda_salt_neg_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_lambda_temp_rapid = np.nansum(area_grids_rapid[~sig_lambda_temp_neg_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_lambda_salt_rapid_2 = np.nansum(area_grids_rapid[~sig_lambda_salt_neg_2_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_lambda_temp_rapid_2 = np.nansum(area_grids_rapid[~sig_lambda_temp_neg_2_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

area_index_lambda_salt_rapid_3 = np.nansum(area_grids_rapid[~sig_lambda_salt_neg_3_rapid[:,0:len(lon_grid)].mask])/total_area_rapid
area_index_lambda_temp_rapid_3 = np.nansum(area_grids_rapid[~sig_lambda_temp_neg_3_rapid[:,0:len(lon_grid)].mask])/total_area_rapid

#%%

#dz and dx of grid points along OSNAP
fh_grid_osnap = netcdf.Dataset(directory + 'OSNAP_lat_lon_depth_area_closest_to.nc','r')

depth_grid 	= fh_grid_osnap.variables['depth'][:]
lon_grid	= fh_grid_osnap.variables['lon'][:]
dx_osnap = fh_grid_osnap.variables['dx'][:]
dz_osnap = fh_grid_osnap.variables['dz'][:]

area_grids_osnap = ma.masked_all((len(depth_grid), len(lon_grid)))

for depth_i in range(len(depth_grid)):
    for lon_i in range(len(lon_grid)):
        area_grids_osnap[depth_i, lon_i] = dx_osnap[lon_i]*dz_osnap[depth_i]

if np.shape(p_var_salt_sig_osnap) != np.shape(area_grids_osnap):
    print('Shapes are not the same!!')
    
total_area_osnap = np.nansum(area_grids_osnap[~p_var_salt_osnap.mask])

area_index_var_salt_osnap = np.nansum(area_grids_osnap[~p_var_salt_sig_osnap.mask])/total_area_osnap
area_index_var_temp_osnap = np.nansum(area_grids_osnap[~p_var_temp_sig_osnap.mask])/total_area_osnap

area_index_var_salt_osnap_2 = np.nansum(area_grids_osnap[~p_var_salt_sig_2_osnap.mask])/total_area_osnap
area_index_var_temp_osnap_2 = np.nansum(area_grids_osnap[~p_var_temp_sig_2_osnap.mask])/total_area_osnap

area_index_var_salt_osnap_3 = np.nansum(area_grids_osnap[~p_var_salt_sig_3_osnap.mask])/total_area_osnap
area_index_var_temp_osnap_3 = np.nansum(area_grids_osnap[~p_var_temp_sig_3_osnap.mask])/total_area_osnap

area_index_ac_salt_osnap = np.nansum(area_grids_osnap[~sig_ac_salt_pos_osnap.mask])/total_area_osnap
area_index_ac_temp_osnap = np.nansum(area_grids_osnap[~sig_ac_temp_pos_osnap.mask])/total_area_osnap

area_index_ac_salt_osnap_2 = np.nansum(area_grids_osnap[~sig_ac_salt_pos_2_osnap.mask])/total_area_osnap
area_index_ac_temp_osnap_2 = np.nansum(area_grids_osnap[~sig_ac_temp_pos_2_osnap.mask])/total_area_osnap

area_index_ac_salt_osnap_3 = np.nansum(area_grids_osnap[~sig_ac_salt_pos_3_osnap.mask])/total_area_osnap
area_index_ac_temp_osnap_3 = np.nansum(area_grids_osnap[~sig_ac_temp_pos_3_osnap.mask])/total_area_osnap

area_index_lambda_salt_osnap = np.nansum(area_grids_osnap[~sig_lambda_salt_neg_osnap.mask])/total_area_osnap
area_index_lambda_temp_osnap = np.nansum(area_grids_osnap[~sig_lambda_temp_neg_osnap.mask])/total_area_osnap

area_index_lambda_salt_osnap_2 = np.nansum(area_grids_osnap[~sig_lambda_salt_neg_2_osnap.mask])/total_area_osnap
area_index_lambda_temp_osnap_2 = np.nansum(area_grids_osnap[~sig_lambda_temp_neg_2_osnap.mask])/total_area_osnap

area_index_lambda_salt_osnap_3 = np.nansum(area_grids_osnap[~sig_lambda_salt_neg_3_osnap.mask])/total_area_osnap
area_index_lambda_temp_osnap_3 = np.nansum(area_grids_osnap[~sig_lambda_temp_neg_3_osnap.mask])/total_area_osnap





