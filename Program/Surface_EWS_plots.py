#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 19:25:13 2024

@author: 6008399

Surface plots of significant regions of the single estimator ratio of EWS of the QE simulation 

Figure 3 of paper, and Figures S6-S8

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

#%% Read in data (window 70)

#%% Restoring rate

#Salinity (E2/E1)
lambda_salt1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12.nc','r')
lambda_salt2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_105m.nc','r')
lambda_salt3 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_209m.nc','r')
lambda_salt4 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_408m.nc','r')
lambda_salt5 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_707m.nc','r')
lambda_salt6 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_1106m.nc','r')
lambda_salt7 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_1968m.nc','r')
lambda_salt8 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_2889m.nc','r')
lambda_salt9 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_SALT_month_1-12_4375m.nc','r')

depth = lambda_salt1.variables['depth'][:]
lon = lambda_salt1.variables['lon'][:]
lat = lambda_salt1.variables['lat'][:]

ratio_lambda_salt1  = lambda_salt1.variables['RATIO_LAMBDA_SALT'][0:1,:,:]
ratio_lambda_salt2  = lambda_salt2.variables['RATIO_LAMBDA_SALT'][10:11,:,:]
ratio_lambda_salt3  = lambda_salt3.variables['RATIO_LAMBDA_SALT'][20:21,:,:]
ratio_lambda_salt4  = lambda_salt4.variables['RATIO_LAMBDA_SALT'][30:31,:,:]
ratio_lambda_salt5  = lambda_salt5.variables['RATIO_LAMBDA_SALT'][36:37,:,:]
ratio_lambda_salt6  = lambda_salt6.variables['RATIO_LAMBDA_SALT'][40:41,:,:]
ratio_lambda_salt7  = lambda_salt7.variables['RATIO_LAMBDA_SALT'][45:46,:,:]
ratio_lambda_salt8  = lambda_salt8.variables['RATIO_LAMBDA_SALT'][49:50,:,:]
ratio_lambda_salt9  = lambda_salt9.variables['RATIO_LAMBDA_SALT'][55:56,:,:]

sig_lambda_salt_neg1  = lambda_salt1.variables['SIG_RATIO_LAMBDA_SALT_neg'][0:1,:,:]
sig_lambda_salt_neg2  = lambda_salt2.variables['SIG_RATIO_LAMBDA_SALT_neg'][10:11,:,:]
sig_lambda_salt_neg3  = lambda_salt3.variables['SIG_RATIO_LAMBDA_SALT_neg'][20:21,:,:]
sig_lambda_salt_neg4  = lambda_salt4.variables['SIG_RATIO_LAMBDA_SALT_neg'][30:31,:,:]
sig_lambda_salt_neg5  = lambda_salt5.variables['SIG_RATIO_LAMBDA_SALT_neg'][36:37,:,:]
sig_lambda_salt_neg6  = lambda_salt6.variables['SIG_RATIO_LAMBDA_SALT_neg'][40:41,:,:]
sig_lambda_salt_neg7  = lambda_salt7.variables['SIG_RATIO_LAMBDA_SALT_neg'][45:46,:,:]
sig_lambda_salt_neg8  = lambda_salt8.variables['SIG_RATIO_LAMBDA_SALT_neg'][49:50,:,:]
sig_lambda_salt_neg9  = lambda_salt9.variables['SIG_RATIO_LAMBDA_SALT_neg'][55:56,:,:]

ratio_lambda_salt = np.ma.concatenate((ratio_lambda_salt1, ratio_lambda_salt2, ratio_lambda_salt3, ratio_lambda_salt4, ratio_lambda_salt5, ratio_lambda_salt6, ratio_lambda_salt7, ratio_lambda_salt8, ratio_lambda_salt9))
sig_lambda_salt_neg = np.ma.concatenate((sig_lambda_salt_neg1, sig_lambda_salt_neg2, sig_lambda_salt_neg3, sig_lambda_salt_neg4, sig_lambda_salt_neg5, sig_lambda_salt_neg6, sig_lambda_salt_neg7, sig_lambda_salt_neg8, sig_lambda_salt_neg9))

#Salinity (E3/E1)
lambda_salt1_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_5m.nc','r')
lambda_salt2_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_105m.nc','r')
lambda_salt3_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_209m.nc','r')
lambda_salt4_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_408m.nc','r')
lambda_salt5_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_707m.nc','r')
lambda_salt6_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_1106m.nc','r')
lambda_salt7_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_1968m.nc','r')
lambda_salt8_E3_E1= netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_2889m.nc','r')
lambda_salt9_E3_E1= netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_600_depth_4375m.nc','r')

ratio_lambda_salt1_E3_E1  = lambda_salt1_E3_E1.variables['RATIO_LAMBDA_SALT'][0:1,:,:]
ratio_lambda_salt2_E3_E1  = lambda_salt2_E3_E1.variables['RATIO_LAMBDA_SALT'][10:11,:,:]
ratio_lambda_salt3_E3_E1  = lambda_salt3_E3_E1.variables['RATIO_LAMBDA_SALT'][20:21,:,:]
ratio_lambda_salt4_E3_E1  = lambda_salt4_E3_E1.variables['RATIO_LAMBDA_SALT'][30:31,:,:]
ratio_lambda_salt5_E3_E1  = lambda_salt5_E3_E1.variables['RATIO_LAMBDA_SALT'][36:37,:,:]
ratio_lambda_salt6_E3_E1  = lambda_salt6_E3_E1.variables['RATIO_LAMBDA_SALT'][40:41,:,:]
ratio_lambda_salt7_E3_E1  = lambda_salt7_E3_E1.variables['RATIO_LAMBDA_SALT'][45:46,:,:]
ratio_lambda_salt8_E3_E1  = lambda_salt8_E3_E1.variables['RATIO_LAMBDA_SALT'][49:50,:,:]
ratio_lambda_salt9_E3_E1  = lambda_salt9_E3_E1.variables['RATIO_LAMBDA_SALT'][55:56,:,:]

sig_lambda_salt_neg1_E3_E1  = lambda_salt1_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][0:1,:,:]
sig_lambda_salt_neg2_E3_E1  = lambda_salt2_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][10:11,:,:]
sig_lambda_salt_neg3_E3_E1  = lambda_salt3_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][20:21,:,:]
sig_lambda_salt_neg4_E3_E1  = lambda_salt4_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][30:31,:,:]
sig_lambda_salt_neg5_E3_E1  = lambda_salt5_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][36:37,:,:]
sig_lambda_salt_neg6_E3_E1  = lambda_salt6_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][40:41,:,:]
sig_lambda_salt_neg7_E3_E1  = lambda_salt7_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][45:46,:,:]
sig_lambda_salt_neg8_E3_E1  = lambda_salt8_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][49:50,:,:]
sig_lambda_salt_neg9_E3_E1  = lambda_salt9_E3_E1.variables['SIG_RATIO_LAMBDA_SALT_neg'][55:56,:,:]

ratio_lambda_salt_E3_E1 = np.ma.concatenate((ratio_lambda_salt1_E3_E1, ratio_lambda_salt2_E3_E1, ratio_lambda_salt3_E3_E1, ratio_lambda_salt4_E3_E1, ratio_lambda_salt5_E3_E1, ratio_lambda_salt6_E3_E1, ratio_lambda_salt7_E3_E1, ratio_lambda_salt8_E3_E1, ratio_lambda_salt9_E3_E1))
sig_lambda_salt_neg_E3_E1 = np.ma.concatenate((sig_lambda_salt_neg1_E3_E1, sig_lambda_salt_neg2_E3_E1, sig_lambda_salt_neg3_E3_E1, sig_lambda_salt_neg4_E3_E1, sig_lambda_salt_neg5_E3_E1, sig_lambda_salt_neg6_E3_E1, sig_lambda_salt_neg7_E3_E1, sig_lambda_salt_neg8_E3_E1, sig_lambda_salt_neg9_E3_E1))

#Temperature (E2/E1)
lambda_temp1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12.nc','r')
lambda_temp2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_105m.nc','r')
lambda_temp3 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_209m.nc','r')
lambda_temp4 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_408m.nc','r')
lambda_temp5 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_707m.nc','r')
lambda_temp6 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_1106m.nc','r')
lambda_temp7 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_1968m.nc','r')
lambda_temp8 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_2889m.nc','r')
lambda_temp9 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_4375m.nc','r')

ratio_lambda_temp1  = lambda_temp1.variables['RATIO_LAMBDA_TEMP'][0:1,:,:]
ratio_lambda_temp2  = lambda_temp2.variables['RATIO_LAMBDA_TEMP'][10:11,:,:]
ratio_lambda_temp3  = lambda_temp3.variables['RATIO_LAMBDA_TEMP'][20:21,:,:]
ratio_lambda_temp4  = lambda_temp4.variables['RATIO_LAMBDA_TEMP'][30:31,:,:]
ratio_lambda_temp5  = lambda_temp5.variables['RATIO_LAMBDA_TEMP'][36:37,:,:]
ratio_lambda_temp6  = lambda_temp6.variables['RATIO_LAMBDA_TEMP'][40:41,:,:]
ratio_lambda_temp7  = lambda_temp7.variables['RATIO_LAMBDA_TEMP'][45:46,:,:]
ratio_lambda_temp8  = lambda_temp8.variables['RATIO_LAMBDA_TEMP'][49:50,:,:]
ratio_lambda_temp9  = lambda_temp9.variables['RATIO_LAMBDA_TEMP'][55:56,:,:]

sig_lambda_temp_neg1  = lambda_temp1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][0:1,:,:]
sig_lambda_temp_neg2  = lambda_temp2.variables['SIG_RATIO_LAMBDA_TEMP_neg'][10:11,:,:]
sig_lambda_temp_neg3  = lambda_temp3.variables['SIG_RATIO_LAMBDA_TEMP_neg'][20:21,:,:]
sig_lambda_temp_neg4  = lambda_temp4.variables['SIG_RATIO_LAMBDA_TEMP_neg'][30:31,:,:]
sig_lambda_temp_neg5  = lambda_temp5.variables['SIG_RATIO_LAMBDA_TEMP_neg'][36:37,:,:]
sig_lambda_temp_neg6  = lambda_temp6.variables['SIG_RATIO_LAMBDA_TEMP_neg'][40:41,:,:]
sig_lambda_temp_neg7  = lambda_temp7.variables['SIG_RATIO_LAMBDA_TEMP_neg'][45:46,:,:]
sig_lambda_temp_neg8  = lambda_temp8.variables['SIG_RATIO_LAMBDA_TEMP_neg'][49:50,:,:]
sig_lambda_temp_neg9  = lambda_temp9.variables['SIG_RATIO_LAMBDA_TEMP_neg'][55:56,:,:]

ratio_lambda_temp = np.ma.concatenate((ratio_lambda_temp1, ratio_lambda_temp2, ratio_lambda_temp3, ratio_lambda_temp4, ratio_lambda_temp5, ratio_lambda_temp6, ratio_lambda_temp7, ratio_lambda_temp8, ratio_lambda_temp9))
sig_lambda_temp_neg = np.ma.concatenate((sig_lambda_temp_neg1, sig_lambda_temp_neg2, sig_lambda_temp_neg3, sig_lambda_temp_neg4, sig_lambda_temp_neg5, sig_lambda_temp_neg6, sig_lambda_temp_neg7, sig_lambda_temp_neg8, sig_lambda_temp_neg9))

#Temperature (E3/E1)
lambda_temp1_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_5m.nc','r')
lambda_temp2_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_105m.nc','r')
lambda_temp3_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_209m.nc','r')
lambda_temp4_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_408m.nc','r')
lambda_temp5_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_707m.nc','r')
lambda_temp6_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_1106m.nc','r')
lambda_temp7_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_1968m.nc','r')
lambda_temp8_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_600_depth_2889m.nc','r')
lambda_temp9_E3_E1 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA1_TEMP_month_1-12_4375m.nc','r')

ratio_lambda_temp1_E3_E1  = lambda_temp1_E3_E1.variables['RATIO_LAMBDA_TEMP'][0:1,:,:]
ratio_lambda_temp2_E3_E1  = lambda_temp2_E3_E1.variables['RATIO_LAMBDA_TEMP'][10:11,:,:]
ratio_lambda_temp3_E3_E1  = lambda_temp3_E3_E1.variables['RATIO_LAMBDA_TEMP'][20:21,:,:]
ratio_lambda_temp4_E3_E1  = lambda_temp4_E3_E1.variables['RATIO_LAMBDA_TEMP'][30:31,:,:]
ratio_lambda_temp5_E3_E1  = lambda_temp5_E3_E1.variables['RATIO_LAMBDA_TEMP'][36:37,:,:]
ratio_lambda_temp6_E3_E1  = lambda_temp6_E3_E1.variables['RATIO_LAMBDA_TEMP'][40:41,:,:]
ratio_lambda_temp7_E3_E1  = lambda_temp7_E3_E1.variables['RATIO_LAMBDA_TEMP'][45:46,:,:]
ratio_lambda_temp8_E3_E1  = lambda_temp8_E3_E1.variables['RATIO_LAMBDA_TEMP'][49:50,:,:]
ratio_lambda_temp9_E3_E1  = lambda_temp9_E3_E1.variables['RATIO_LAMBDA_TEMP'][55:56,:,:]

sig_lambda_temp_neg1_E3_E1  = lambda_temp1_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][0:1,:,:]
sig_lambda_temp_neg2_E3_E1  = lambda_temp2_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][10:11,:,:]
sig_lambda_temp_neg3_E3_E1  = lambda_temp3_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][20:21,:,:]
sig_lambda_temp_neg4_E3_E1  = lambda_temp4_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][30:31,:,:]
sig_lambda_temp_neg5_E3_E1  = lambda_temp5_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][36:37,:,:]
sig_lambda_temp_neg6_E3_E1  = lambda_temp6_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][40:41,:,:]
sig_lambda_temp_neg7_E3_E1  = lambda_temp7_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][45:46,:,:]
sig_lambda_temp_neg8_E3_E1  = lambda_temp8_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][49:50,:,:]
sig_lambda_temp_neg9_E3_E1  = lambda_temp9_E3_E1.variables['SIG_RATIO_LAMBDA_TEMP_neg'][55:56,:,:]

ratio_lambda_temp_E3_E1 = np.ma.concatenate((ratio_lambda_temp1_E3_E1, ratio_lambda_temp2_E3_E1, ratio_lambda_temp3_E3_E1, ratio_lambda_temp4_E3_E1, ratio_lambda_temp5_E3_E1, ratio_lambda_temp6_E3_E1, ratio_lambda_temp7_E3_E1, ratio_lambda_temp8_E3_E1, ratio_lambda_temp9_E3_E1))
sig_lambda_temp_neg_E3_E1 = np.ma.concatenate((sig_lambda_temp_neg1_E3_E1, sig_lambda_temp_neg2_E3_E1, sig_lambda_temp_neg3_E3_E1, sig_lambda_temp_neg4_E3_E1, sig_lambda_temp_neg5_E3_E1, sig_lambda_temp_neg6_E3_E1, sig_lambda_temp_neg7_E3_E1, sig_lambda_temp_neg8_E3_E1, sig_lambda_temp_neg9_E3_E1))


#%% Plot (Figure 3 of paper)

depth_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8]

fig, axs = plt.subplots(3,3,subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,10.5))

#axs is a 2D array of 'GeoAxes'. Flatten it into a 1D array
axs=axs.flatten()

for i, depth_i in enumerate(depth_idx):
    
    print(depth_i)
    print(i)

    data_EWS = sig_lambda_salt_neg_E3_E1[depth_i,:,:]
    
    print(data_EWS.count()/np.size(sig_lambda_temp_neg[0,:,:]))
    
    #Contour plot
    cs2 = axs[i].contourf(lon, lat, data_EWS, transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))

    axs[i].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
    axs[i].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
    axs[i].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
    #Title each subplot with the depth
    axs[i].set_title(''+str(int(depth[depth_idx_1[i]]))+'m')
    
    #Draw coastlines for each subplot
    axs[i].coastlines()
    
    axs[i].set_xlim([-80,20])
    
    # Longitude labels
    axs[i].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axs[i].xaxis.set_major_formatter(lon_formatter)

    # Latitude labels
    axs[i].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axs[i].yaxis.set_major_formatter(lat_formatter)

# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.25, top=0.93, left=0.05, right=0.99,
                    wspace=0.02, hspace=0.3)

# Add a colorbar axis at the bottom of the graph
#cbar_ax = fig.add_axes([0.15, 0.2, 0.7, 0.02])

# Draw the colorbar
#cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal')

# Add a big title at the top
plt.suptitle('Significant $R^S_\mathrm{RES}$ regions (E2/E1)')
#plt.savefig(directory_figures + 'Atlantic_LAMBDA_regions_window_70_sig_095_TEMP_single_estimator_branch_1650_600.pdf')
plt.show()


#%% Difference in surface plots for using different branches

lambda_salt1_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_1500_depth_5m.nc','r')
lambda_salt3_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_1500_depth_209m.nc','r')
lambda_salt6_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_1500_depth_1106m.nc','r')
lambda_salt7_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_SALT_month_1-12_surface_branch_1650_1500_depth_1968m.nc','r')

sig_lambda_salt_neg1_E3_E2  = lambda_salt1_E3_E2.variables['SIG_RATIO_LAMBDA_SALT_neg'][0:1,:,:]
sig_lambda_salt_neg3_E3_E2  = lambda_salt3_E3_E2.variables['SIG_RATIO_LAMBDA_SALT_neg'][20:21,:,:]
sig_lambda_salt_neg6_E3_E2  = lambda_salt6_E3_E2.variables['SIG_RATIO_LAMBDA_SALT_neg'][40:41,:,:]
sig_lambda_salt_neg7_E3_E2  = lambda_salt7_E3_E2.variables['SIG_RATIO_LAMBDA_SALT_neg'][45:46,:,:]

fig, axs = plt.subplots(3,3,subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,10.5))

axs[0,0].contourf(lon, lat, sig_lambda_salt_neg1[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[0,0].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[0,0].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[0,0].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[0,0].set_title('E2/E1 (5m)')
axs[0,0].coastlines()
axs[0,0].set_xlim([-80,20])
    
axs[0,0].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,0].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[0,0].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,0].yaxis.set_major_formatter(lat_formatter)

axs[0,1].contourf(lon, lat, sig_lambda_salt_neg1_E3_E1[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[0,1].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[0,1].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[0,1].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[0,1].set_title('E3/E1 (5m)')
axs[0,1].coastlines()
axs[0,1].set_xlim([-80,20])
    
axs[0,1].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,1].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[0,1].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,1].yaxis.set_major_formatter(lat_formatter)

axs[0,2].contourf(lon, lat, sig_lambda_salt_neg1_E3_E2[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[0,2].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[0,2].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[0,2].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[0,2].set_title('E3/E2 (5m)')
axs[0,2].coastlines()
axs[0,2].set_xlim([-80,20])
    
axs[0,2].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,2].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[0,2].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,2].yaxis.set_major_formatter(lat_formatter)

axs[1,0].contourf(lon, lat, sig_lambda_salt_neg3[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[1,0].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[1,0].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[1,0].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[1,0].set_title('E2/E1 (209m)')
axs[1,0].coastlines()
axs[1,0].set_xlim([-80,20])
    
axs[1,0].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,0].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[1,0].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,0].yaxis.set_major_formatter(lat_formatter)

axs[1,1].contourf(lon, lat, sig_lambda_salt_neg3_E3_E1[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[1,1].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[1,1].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[1,1].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[1,1].set_title('E3/E1 (209m)')
axs[1,1].coastlines()
axs[1,1].set_xlim([-80,20])
    
axs[1,1].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,1].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[1,1].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,1].yaxis.set_major_formatter(lat_formatter)

axs[1,2].contourf(lon, lat, sig_lambda_salt_neg3_E3_E2[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[1,2].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[1,2].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[1,2].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[1,2].set_title('E3/E2 (209m)')
axs[1,2].coastlines()
axs[1,2].set_xlim([-80,20])
    
axs[1,2].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,2].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[1,2].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,2].yaxis.set_major_formatter(lat_formatter)

axs[2,0].contourf(lon, lat, sig_lambda_salt_neg7[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[2,0].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[2,0].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[2,0].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[2,0].set_title('E2/E1 (1967m)')
axs[2,0].coastlines()
axs[2,0].set_xlim([-80,20])
    
axs[2,0].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[2,0].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[2,0].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[2,0].yaxis.set_major_formatter(lat_formatter)

axs[2,1].contourf(lon, lat, sig_lambda_salt_neg7_E3_E1[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[2,1].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[2,1].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[2,1].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[2,1].set_title('E3/E1 (1967m)')
axs[2,1].coastlines()
axs[2,1].set_xlim([-80,20])
    
axs[2,1].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[2,1].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[2,1].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[2,1].yaxis.set_major_formatter(lat_formatter)

axs[2,2].contourf(lon, lat, sig_lambda_salt_neg7_E3_E2[0,:,:], transform=ccrs.PlateCarree(), colors=['green'])#, levels=np.linspace(0,1,2))
axs[2,2].axhline(y=60, color='k', linestyle='--', linewidth=0.7)
axs[2,2].axhline(y=26, color='k', linestyle='--', linewidth=0.7)
axs[2,2].axhline(y=-34.5, color='k', linestyle='--', linewidth=0.7)
    
axs[2,2].set_title('E3/E2 (1967m)')
axs[2,2].coastlines()
axs[2,2].set_xlim([-80,20])
    
axs[2,2].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[2,2].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[2,2].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[2,2].yaxis.set_major_formatter(lat_formatter)

plt.suptitle('Significant R$_{RES}^S$ regions')

# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.25, top=0.93, left=0.05, right=0.99,
                    wspace=0.02, hspace=0.3)

plt.savefig(directory_figures + 'Atlantic_DIFF_LAMBDA_regions_window_70_sig_095_SALT_single_estimator_branches_E1_E2_E3.pdf')

#%% temperature

lambda_temp1_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_1500_depth_5m.nc','r')
lambda_temp3_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_1500_depth_209m.nc','r')
lambda_temp6_E3_E2 = netcdf.Dataset(directory + 'Atlantic_single_estimator_LAMBDA_TEMP_month_1-12_surface_branch_1650_1500_depth_1106m.nc','r')

sig_lambda_temp_neg1_E3_E2  = lambda_temp1_E3_E2.variables['SIG_RATIO_LAMBDA_TEMP_neg'][0:1,:,:]
sig_lambda_temp_neg3_E3_E2  = lambda_temp3_E3_E2.variables['SIG_RATIO_LAMBDA_TEMP_neg'][20:21,:,:]
sig_lambda_temp_neg6_E3_E2  = lambda_temp6_E3_E2.variables['SIG_RATIO_LAMBDA_TEMP_neg'][40:41,:,:]

fig, axs = plt.subplots(3,3,subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,10.5))

axs[0,0].contourf(lon, lat, sig_lambda_temp_neg1[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[0,0].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[0,0].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[0,0].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[0,0].set_title('E2/E1 (5m)')
axs[0,0].coastlines()
axs[0,0].set_xlim([-80,20])
    
axs[0,0].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,0].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[0,0].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,0].yaxis.set_major_formatter(lat_formatter)

axs[0,1].contourf(lon, lat, sig_lambda_temp_neg1_E3_E1[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[0,1].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[0,1].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[0,1].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[0,1].set_title('E3/E1 (5m)')
axs[0,1].coastlines()
axs[0,1].set_xlim([-80,20])
    
axs[0,1].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,1].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[0,1].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,1].yaxis.set_major_formatter(lat_formatter)

axs[0,2].contourf(lon, lat, sig_lambda_temp_neg1_E3_E2[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[0,2].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[0,2].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[0,2].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[0,2].set_title('E3/E2 (5m)')
axs[0,2].coastlines()
axs[0,2].set_xlim([-80,20])
    
axs[0,2].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,2].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[0,2].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,2].yaxis.set_major_formatter(lat_formatter)

axs[1,0].contourf(lon, lat, sig_lambda_temp_neg3[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[1,0].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[1,0].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[1,0].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[1,0].set_title('E2/E1 (209m)')
axs[1,0].coastlines()
axs[1,0].set_xlim([-80,20])
    
axs[1,0].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,0].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[1,0].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,0].yaxis.set_major_formatter(lat_formatter)

axs[1,1].contourf(lon, lat, sig_lambda_temp_neg3_E3_E1[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[1,1].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[1,1].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[1,1].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[1,1].set_title('E3/E1 (209m)')
axs[1,1].coastlines()
axs[1,1].set_xlim([-80,20])
    
axs[1,1].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,1].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[1,1].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,1].yaxis.set_major_formatter(lat_formatter)

axs[1,2].contourf(lon, lat, sig_lambda_temp_neg3_E3_E2[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[1,2].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[1,2].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[1,2].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[1,2].set_title('E3/E2 (209m)')
axs[1,2].coastlines()
axs[1,2].set_xlim([-80,20])
    
axs[1,2].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,2].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[1,2].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,2].yaxis.set_major_formatter(lat_formatter)

axs[2,0].contourf(lon, lat, sig_lambda_temp_neg6[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[2,0].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[2,0].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[2,0].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[2,0].set_title('E2/E1 (1106m)')
axs[2,0].coastlines()
axs[2,0].set_xlim([-80,20])
    
axs[2,0].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[2,0].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[2,0].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[2,0].yaxis.set_major_formatter(lat_formatter)

axs[2,1].contourf(lon, lat, sig_lambda_temp_neg6_E3_E1[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[2,1].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[2,1].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[2,1].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[2,1].set_title('E3/E1 (1106m)')
axs[2,1].coastlines()
axs[2,1].set_xlim([-80,20])
    
axs[2,1].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[2,1].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[2,1].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[2,1].yaxis.set_major_formatter(lat_formatter)

axs[2,2].contourf(lon, lat, sig_lambda_temp_neg6_E3_E2[0,:,:], transform=ccrs.PlateCarree(), colors=['blue'])#, levels=np.linspace(0,1,2))
axs[2,2].axhline(y=60, color='r', linestyle='--', linewidth=0.7)
axs[2,2].axhline(y=26, color='r', linestyle='--', linewidth=0.7)
axs[2,2].axhline(y=-34.5, color='r', linestyle='--', linewidth=0.7)
    
axs[2,2].set_title('E3/E2 (1106m)')
axs[2,2].coastlines()
axs[2,2].set_xlim([-80,20])
    
axs[2,2].set_xticks(np.arange(-70,20,40), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[2,2].xaxis.set_major_formatter(lon_formatter)

# Latitude labels
axs[2,2].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[2,2].yaxis.set_major_formatter(lat_formatter)

plt.suptitle('Significant R$_{RES}^T$ regions')

# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.25, top=0.93, left=0.05, right=0.99,
                    wspace=0.02, hspace=0.3)

plt.savefig(directory_figures + 'Atlantic_DIFF_LAMBDA_regions_window_70_sig_095_TEMP_single_estimator_branches_E1_E2_E3.pdf')

