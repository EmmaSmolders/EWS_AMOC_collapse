#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 09:19:06 2025

@author: 6008399

Locations of estimated tipping times CESM QE

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
from scipy.io import loadmat

#Making pathway to folder with all data
directory = r'/../../Data/EWS/'
directory_figures = r'/Users/6008399/Documents/PhD/Cesm_collapse/Figures/'

#%% Ratio of variance using single estimator (E3/E1)

var_salt_3 = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_SALT_month_1-12_34S_branch_1650_600.nc','r')
var_temp_3 = netcdf.Dataset(directory + 'Atlantic_single_estimator_variance_TEMP_month_1-12_34S_branch_1650_600.nc','r')
      
ratio_var_salt  = var_salt_3.variables['SIGMA_SALT'][:]
ratio_var_temp  = var_temp_3.variables['SIGMA_TEMP'][:]

p_var_salt  = var_salt_3.variables['p_value_LEVENE'][:]
p_var_temp  = var_temp_3.variables['p_value_LEVENE'][:]

lat = var_salt_3.variables['lat'][:]
lon = var_salt_3.variables['lon'][:]
depth = var_salt_3.variables['depth'][:]

p_var_salt_sig = ma.masked_all((len(depth), len(lon)))
p_var_temp_sig = ma.masked_all((len(depth), len(lon)))

#Create mask where the p-value is smaller than 0.01
for i in range(len(depth)):
    for j in range(len(lon)):
        if p_var_salt[i,j] < 0.05 and ratio_var_salt[i,j] > 1:
            #print('hoi')
            p_var_salt_sig[i,j] = 1
        else:
            continue
        
        if p_var_temp[i,j] < 0.05 and ratio_var_temp[i,j] > 1:
            #print('hoi')
            p_var_temp_sig[i,j] = 1
        else:
            continue

#%% Ratio of AC1 using single estimator

ac_salt_3   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_SALT_month_1-12_34S_branch_1650_600.nc','r')
ac_temp_3   = netcdf.Dataset(directory + 'Atlantic_single_estimator_AC1_TEMP_month_1-12_34S_branch_1650_600.nc','r')
       
ratio_ac_salt     = ac_salt_3.variables['RATIO_AC1_SALT'][:]
ratio_ac_temp     = ac_temp_3.variables['RATIO_AC1_TEMP'][:]

sig_ac_salt_pos   = ac_salt_3.variables['SIG_RATIO_AC1_SALT_pos'][:]
sig_ac_salt_neg   = ac_salt_3.variables['SIG_RATIO_AC1_SALT_neg'][:]

sig_ac_temp_pos   = ac_temp_3.variables['SIG_RATIO_AC1_TEMP_pos'][:]
sig_ac_temp_neg   = ac_temp_3.variables['SIG_RATIO_AC1_TEMP_neg'][:]

#%% Ratio of lambda using single estimator

lambda_salt_3   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_34S_branch_1650_600.nc','r')
lambda_temp_3   = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_34S_branch_1650_600.nc','r')

ratio_lambda_salt     = lambda_salt_3.variables['RATIO_LAMBDA_SALT'][:]
ratio_lambda_temp     = lambda_temp_3.variables['RATIO_LAMBDA_TEMP'][:]

sig_lambda_salt_pos  = lambda_salt_3.variables['SIG_RATIO_LAMBDA_SALT_pos'][:]
sig_lambda_salt_neg  = lambda_salt_3.variables['SIG_RATIO_LAMBDA_SALT_neg'][:]

sig_lambda_temp_pos  = lambda_temp_3.variables['SIG_RATIO_LAMBDA_TEMP_pos'][:]
sig_lambda_temp_neg  = lambda_temp_3.variables['SIG_RATIO_LAMBDA_TEMP_neg'][:]

#%%

#dz and dx of grid points along SAMBA
fh_grid_samba = netcdf.Dataset(directory + 'SAMBA_lon_lat_depth_area.nc','r')

depth_grid 	= fh_grid_samba.variables['depth'][:]
lon_grid	= fh_grid_samba.variables['lon'][:]
dx_samba = fh_grid_samba.variables['dx'][:]
dz_samba = fh_grid_samba.variables['dz'][:]

fh = netcdf.Dataset(directory + 'lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','r')
fh_temp = netcdf.Dataset(directory + 'lambda_TEMP_transient_1500_endtimes_1700_1702_spacing_2_window_70_SAMBA_all_gridpoints.nc','r')

depth             = fh.variables['depth'][:]          #Depth [m]
lon               = fh.variables['lon'][0:len(lon_grid)]          #Longitude [degE]
lambda_salt       = fh.variables['lambda_values_1700'][:,:,0:len(lon_grid)] 
time              = fh.variables['time_1700'][:]

fh.close()

if lon[0] != lon_grid[0]:
    print('Longitude ranges are not equal! Change this!!')

#Change point times (from Matlab procedure)
data = loadmat(directory + 'changepoint_time_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
time_change = data['change_point_time_all_start']
time_change = time_change[:,:,0:len(lon_grid)]
time_change = ma.masked_invalid(time_change)

#Times at zero (from Matlab procedure)
data_zero = loadmat(directory + 'time_at_zero_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
time_zero = data_zero['time_at_zero_all_start']
time_zero = time_zero[:,:,0:len(lon_grid)]
time_zero = ma.masked_invalid(time_zero)

data_zero_temp = loadmat(directory + 'time_at_zero_all_SAMBA_TEMP_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
time_zero_temp = data_zero_temp['time_at_zero_all_start']
time_zero_temp = time_zero_temp[:,:,0:len(lon_grid)]
time_zero_temp = ma.masked_invalid(time_zero_temp)

#Range of CPend
end_year_change = np.linspace(1580, 1645, 1646-1580) 

#Now the negative slopes are indicated by 9e10 ('inf'), but this is not handy for the median determination. So make it smaller
time_zero[time_zero == 9e10] = 9e6
time_zero_temp[time_zero_temp == 9e10] = 9e6

#%%

cmap = colors.ListedColormap(['white'])

#divnorm = mcolors.TwoSlopeNorm(vmin=1600, vcenter=800000, vmax=9000000)

fig, axs = plt.subplots(1, 2, figsize=(14, 6))
axs[0].contourf(lon, depth, ratio_var_salt, cmap=cmap)
axs[0].contourf(lon, depth, time_zero[0,0,:,:], cmap='seismic', levels=np.linspace(1600,15000000,3))
axs[0].contourf(lon, depth, sig_lambda_salt_neg, colors='none', linewidth=1.0, hatches='xx')
axs[0].set_ylabel('Depth [m]', fontsize = 14)
#axs[0,0].set_xlabel('Longitude $^\circ$E', fontsize = 12)
axs[0].set_title('a) Tipping time regions CP$_{end}$ = 1580 (RES$^S$ at 34$^\circ$S)', fontsize = 16)
axs[0].set_ylim(depth[-1], 0)
axs[0].set_xlim(-60,20)
axs[0].set_xticks([-60, -40, -20, 0, 20])
axs[0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0].tick_params(axis='both', which='major', labelsize=12)
axs[0].set_facecolor('lightgrey')
axs[0].grid()

axs[1].contourf(lon, depth, ratio_var_salt, cmap=cmap)
CS = axs[1].contourf(lon, depth, time_zero[0,-1,:,:], cmap='seismic', levels=np.linspace(1600,15000000,3))
axs[1].contourf(lon, depth, sig_lambda_salt_neg, colors='none', linewidth=1.0, hatches='xx')
#axs[1].set_ylabel('Depth [m]', fontsize = 14)
#axs[0,0].set_xlabel('Longitude $^\circ$E', fontsize = 12)
axs[1].set_title('b) Tipping time regions CP$_{end}$ = 1645 (RES$^S$ at 34$^\circ$S)', fontsize = 16)
axs[1].set_ylim(depth[-1], 0)
axs[1].set_xlim(-60,20)
axs[1].set_xticks([-60, -40, -20, 0, 20])
axs[1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1].tick_params(axis='both', which='major', labelsize=12)
axs[1].set_facecolor('lightgrey')
axs[1].grid()

#fig.colorbar(CS)

plt.tight_layout()
plt.savefig(directory_figures +'tipping_time_regions_LAMBDA_salt_CPbegin_CPend.png', dpi=400)
plt.show()

#%%

cmap = colors.ListedColormap(['white'])

fig, axs = plt.subplots(1, 2, figsize=(14, 6))
axs[0].contourf(lon, depth, ratio_var_temp, cmap=cmap)
axs[0].contourf(lon, depth, time_zero_temp[0,0,:,:], cmap='seismic', levels=np.linspace(1600,15000000,3))
axs[0].contourf(lon, depth, sig_lambda_temp_neg, colors='none', linewidth=1.0, hatches='xx')
axs[0].set_ylabel('Depth [m]', fontsize = 14)
#axs[0,0].set_xlabel('Longitude $^\circ$E', fontsize = 12)
axs[0].set_title('a) Tipping time regions CP$_{end}$ = 1580 (RES$^T$ at 34$^\circ$S)', fontsize = 16)
axs[0].set_ylim(depth[-1], 0)
axs[0].set_xlim(-60,20)
axs[0].set_xticks([-60, -40, -20, 0, 20])
axs[0].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[0].tick_params(axis='both', which='major', labelsize=12)
axs[0].set_facecolor('lightgrey')
axs[0].grid()

axs[1].contourf(lon, depth, ratio_var_temp, cmap=cmap)
CS = axs[1].contourf(lon, depth, time_zero_temp[0,-1,:,:], levels=np.linspace(1600,15000000,3))
axs[1].contourf(lon, depth, sig_lambda_temp_neg, colors='none', linewidth=1.0, hatches='xx')
#axs[1].set_ylabel('Depth [m]', fontsize = 14)
#axs[0,0].set_xlabel('Longitude $^\circ$E', fontsize = 12)
axs[1].set_title('b) Tipping time regions CP$_{end}$ = 1645 (RES$^T$ at 34$^\circ$S)', fontsize = 16)
axs[1].set_ylim(depth[-1], 0)
axs[1].set_xlim(-60,20)
axs[1].set_xticks([-60, -40, -20, 0, 20])
axs[1].set_xticklabels(['60$^{\circ}$W', '40$^{\circ}$W', '20$^{\circ}$W', '0$^{\circ}$', '20$^{\circ}$E'])
axs[1].tick_params(axis='both', which='major', labelsize=12)
axs[1].set_facecolor('lightgrey')
axs[1].grid()

plt.tight_layout()
plt.savefig(directory_figures +'tipping_time_regions_LAMBDA_temp_CPbegin_CPend.png', dpi=400)
plt.show()

