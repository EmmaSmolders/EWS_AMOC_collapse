#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 12:39:23 2024

@author: 6008399

Changepoint analysis for QE using all gridpoints 

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
from scipy.optimize import fsolve
from numpy.polynomial.polynomial import Polynomial
from scipy.io import loadmat
from statsmodels.graphics.gofplots import qqplot

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/Cesm_collapse/netcdf/'
directory_figures = r'/Users/6008399/Documents/PhD/Cesm_collapse/Figures/'

#%% Read in data 

#dz and dx of grid points along SAMBA
fh_grid_samba = netcdf.Dataset(directory + 'SAMBA_lon_lat_depth_area.nc','r')

depth_grid 	= fh_grid_samba.variables['depth'][:]
lon_grid	= fh_grid_samba.variables['lon'][:]
dx_samba = fh_grid_samba.variables['dx'][:]
dz_samba = fh_grid_samba.variables['dz'][:]

fh = netcdf.Dataset(directory + 'lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','r')

depth             = fh.variables['depth'][:]          #Depth [m]
lon               = fh.variables['lon'][0:len(lon_grid)]          #Longitude [degE]
lambda_salt       = fh.variables['lambda_values_1700'][:,:,0:len(lon_grid)] 
lambda_salt       = -lambda_salt #Restoring rate is negative in Boers but should be positive
time              = fh.variables['time_1700'][:]

fh.close()

if lon[0] != lon_grid[0]:
    print('Longitude ranges are not equal! Change this!!')

#Change point times (from Matlab procedure)
data = loadmat(directory + 'changepoint_time_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
time_change = data['change_point_time_all_start']
time_change = time_change[:,:,:,0:len(lon_grid)]
time_change = ma.masked_invalid(time_change)

#Times at zero (from Matlab procedure)
data_zero = loadmat(directory + 'time_at_zero_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
time_zero = data_zero['time_at_zero_all_start']
time_zero = time_zero[:,:,:,0:len(lon_grid)]
time_zero = ma.masked_invalid(time_zero)

fh_r = loadmat(directory + 'r_value_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
r_value = fh_r['r_value_all_start']
r_value = r_value[:,:,:,0:len(lon_grid)]
r_value = ma.masked_invalid(r_value)

fh_trend = loadmat(directory + 'trend_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
trend = fh_trend['trend_all_start']
trend = trend[:,:,:,0:len(lon_grid)]
trend = ma.masked_invalid(trend)

fh_base = loadmat(directory + 'base_all_SAMBA_SALT_QE_1500_1700_window_70_change_1580_1645_all_CP_pertimeseries_NOSELECTIONBIAS.mat')
base = fh_base['base_all_start']
base = base[:,:,:,0:len(lon_grid)]
base = ma.masked_invalid(base)

#Range of CPend
end_year_change = np.linspace(1580, 1645, 1646-1580) 

if len(end_year_change) != len(time_zero[0,:,0,0]):
    print('CPend ranges are not equal! Change this!!')

#%% Determine p-value and r-value for all changepoint times

window = 70

p = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
p_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
#trend = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
trend_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
#base = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
base_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
depth_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
lon_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
#r_value = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
r_value_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
time_zero_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
time_zero_positive = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))
time_change_good = ma.masked_all((10, len(end_year_change), len(depth), len(lon)))

for cp_i in range(10):
    for i in range(len(end_year_change)):
        print(i)
    
        if len(end_year_change) > 1:
            for depth_i in range(len(depth)):
                for lon_i in range(len(lon)):
                
                    if time_change[cp_i, i, depth_i, lon_i] is ma.masked:
                        continue
            
                    else:
                        
                        time_idx = np.where(time == time_change[cp_i,i,depth_i,lon_i])

                        data_test = lambda_salt[time_idx[0][0]:-(window//2),int(depth_i), int(lon_i)]
                        time_test = time[time_idx[0][0]:-(window//2)]

                        #Determine the shuffled Fourier Kendall test (following Boers 2021)
                        #p[cp_i,i,depth_i,lon_i], trend[cp_i,i,depth_i, lon_i], base[cp_i,i,depth_i,lon_i], r_value[cp_i,i,depth_i,lon_i]	= Kendall_tau_test(data_test, 1000, time_test)
    
                    #Only keep fits with p < 0.05 (monotonically increasing) and r > 0.7 (goodness of fit)
                    #if p[cp_i,i,depth_i,lon_i] < 0.05 and r_value[cp_i,i,depth_i,lon_i] > 0.7:
                    if r_value[cp_i,i,depth_i,lon_i] < 0:
                        time_zero_positive[cp_i,i,depth_i,lon_i] = time_zero[cp_i,i,depth_i,lon_i]
                        
                    if r_value[cp_i,i,depth_i,lon_i] < -0.7 and time_zero[cp_i,i,depth_i,lon_i] > 1665:
                        r_value_good[cp_i,i,depth_i,lon_i] = r_value[cp_i,i,depth_i,lon_i]
                        p_good[cp_i,i,depth_i,lon_i] = p[cp_i,i,depth_i,lon_i]
                        trend_good[cp_i,i,depth_i,lon_i] = trend[cp_i,i,depth_i,lon_i]
                        base_good[cp_i,i,depth_i,lon_i] = base[cp_i,i,depth_i,lon_i]
                        time_zero_good[cp_i,i,depth_i,lon_i] = time_zero[cp_i,i,depth_i,lon_i]
                        time_change_good[cp_i,i,depth_i,lon_i] = time_change[cp_i,i,depth_i,lon_i]

  
#%%

plt.figure(figsize=(8,4))
plt.hist(time_zero_good.flatten(), bins=100, alpha=0.5, color='b', density=True)
#ci_lower, ci_upper = np.percentile(time_zero_good.compressed(), [2.5, 97.5])
#plt.axvline(x=ci_lower, color='r', linestyle='dashed', linewidth=2, label='95% confidence interval')
#plt.axvline(x=ci_upper, color='r', linestyle='dashed', linewidth=2)
plt.grid()
plt.xlabel('Model years')
plt.ylabel('Frequency')
plt.title('Estimated tipping times (salinity) for all fits (R>0)')

plt.figure(figsize=(8,4))
plt.hist(time_zero.flatten(), bins=100, alpha=0.5, color='b', density=True)
plt.grid()
plt.xlabel('Model years')
plt.ylabel('Frequency')
plt.title('Estimated tipping times (salinity) for all fits')

fig, axs = plt.subplots(1, 1, figsize=(8,6))
#ci_lower, ci_upper = np.nanpercentile(time_zero.compressed(), [2.5, 97.5])
axs.hist(r_value_good.flatten(), bins=100, alpha=0.5, color='b', density=False)
axs.set_title('R values (salinity, window 70 years)')
axs.set_ylabel('Frequency')
axs.set_xlabel('R-value')
axs.legend()
axs.grid()

#%% 

area_grids = ma.masked_all((10,len(depth), len(lon)))

for depth_i in range(len(depth)):
    for lon_i in range(len(lon)):
        area_grids[:,depth_i, lon_i] = dx_samba[lon_i]*dz_samba[depth_i]

area_grids	= area_grids / area_grids.max()

#%% Normalize to give same weight to all points (for time_zero_good, so with some R restriction)

number_bins = 100 

normalized_hist = ma.masked_all((len(time_zero_good[0,:,1,1]), number_bins))
bins = ma.masked_all((len(time_zero_good[0,:,1,1]), number_bins + 1))
#final_hist = ma.masked_all(len(time_change[0]))

for i in range(len(time_zero_good[0,:,1,1])):
    
    # Create a mask and use the non-mask values only for the histogram
    data = time_zero_good[:,i,:,:].flatten()

    normalized_hist[i,:], bins[i,:] = np.histogram(data, bins=number_bins, range=(np.nanmin(time_zero_good), np.nanmax(time_zero_good)), density=True, weights=area_grids.flatten())
    print(np.sum(normalized_hist[i,:] * np.diff(bins[i,:]))) #should be very close to 1

#Take sum of normalized histograms to get final histogram (taking the average gives exactly the same shape)    
final_hist = np.nansum(normalized_hist, axis=0)

plt.figure()
plt.title('Summed histogram of normalized PDFs salinity')
plt.stairs(final_hist, bins[-1,:], fill=True)
#plt.xlim(1650, 3000)
#plt.ylim(0,0.1)
plt.grid()

#%%

indices = [0, 10, 20, 30, 40, 50, 63, 75, 100]
#indices = [27, 28, 29, 30, 31, 32, 33, 34, 35]

fig, axs = plt.subplots(3, 3, figsize=(12, 10))

for i, index in enumerate(indices):
    print(index)
    # Flatten the 2D matrix into a 1D array
    time_at_zero_all_flat = time_zero_good[:,index, :, :].flatten()

    # Create a subplot
    ax = axs[i // 3, i % 3]
    ax.hist(time_at_zero_all_flat, bins=50)
    #ax.set_title(f'({np.count_nonzero(~np.isnan(time_zero_good[index, :, :]))} points, window {window}y, CP < y{end_year_change[index]})')
    ax.set_title(f'({time_zero_good[:,index, :, :].count()} points, window {window}y, CP < y{end_year_change[index]})')
    #ax.set_xlim([2019, 2100])

    # Calculate the 95% confidence interval
    ci_lower = np.percentile(time_at_zero_all_flat, 2.5)
    ci_upper = np.percentile(time_at_zero_all_flat, 97.5)

    # Add vertical dashed red lines to the histogram
    #ax.axvline(ci_lower, color='r', linestyle='--')
    #ax.axvline(ci_upper, color='r', linestyle='--')

plt.tight_layout()
#plt.savefig(directory_figures + 'PDFs_salinity_samba_transient_range_CP_restriction_r_07_p_005_RCP85_window10.pdf')
plt.show()

#%%

fig, axs = plt.subplots(3, 3, figsize=(12, 10))

for i, index in enumerate(indices):
    print(index)
    # Flatten the 2D matrix into a 1D array
    time_at_zero_all_flat = time_zero[:,index, :, :].flatten()

    # Create a subplot
    ax = axs[i // 3, i % 3]
    ax.hist(time_at_zero_all_flat, bins=50)
    #ax.set_title(f'({np.count_nonzero(~np.isnan(time_zero_good[index, :, :]))} points, window {window}y, CP < y{end_year_change[index]})')
    ax.set_title(f'({time_zero[:,index, :, :].count()} points, window {window}y, CP < y{end_year_change[index]})')
    #ax.set_xlim([2019, 2100])

    # Calculate the 95% confidence interval
    ci_lower = np.percentile(time_at_zero_all_flat, 2.5)
    ci_upper = np.percentile(time_at_zero_all_flat, 97.5)

    # Add vertical dashed red lines to the histogram
    #ax.axvline(ci_lower, color='r', linestyle='--')
    #ax.axvline(ci_upper, color='r', linestyle='--')

plt.tight_layout()
#plt.savefig(directory_figures + 'PDFs_salinity_samba_transient_range_CP_restriction_r_07_p_005_RCP85_window10.pdf')
plt.show()

#%% Plot fits 

from itertools import chain

closest_indices = np.argsort(np.abs(time_zero_good[0,-1,:,:] - 1759))
closest_indices = closest_indices[12::]

# Assuming time_zero_good is your 2D numpy array
diff = np.abs(time_zero_good[0,-1,:,:] - 1759)  # Find the absolute difference between array values and 2050

# Initialize lists to store the best and worst indices
best_indices = []
worst_indices = []

# Get the 9 best (closest) and worst (furthest) indices
for _ in range(12):
    best_index = np.nanargmin(diff)
    worst_index = np.nanargmax(diff)
    best_indices.append(best_index)
    worst_indices.append(worst_index)
    
    # Set the current best and worst values to NaN so they are not selected again in the next iteration
    diff.flat[best_index] = np.nan
    diff.flat[worst_index] = np.nan

# Convert the indices back to 2D
best_indices_2d = np.unravel_index(best_indices, time_zero_good[0,-1,:,:].shape)
worst_indices_2d = np.unravel_index(worst_indices, time_zero_good[0,-1,:,:].shape)

fig, axs = plt.subplots(3, 3, figsize=(18, 8))

for i in range(0,9):
#for i in chain(range(0, 4), range(6, 10)):
    depth_i = int(worst_indices_2d[0][i]) 
    lon_i = int(worst_indices_2d[1][i])
    
    data_point = lambda_salt[:, depth_i, lon_i]
    time_point = time
    
    trend_plot = trend[0,-1, depth_i, lon_i]
    base_plot = base[0,-1,  depth_i, lon_i]
    time_zero_plot = int(time_zero_good[0,-1,  depth_i, lon_i])
    time_change_plot = int(time_change_good[0,-1,  depth_i, lon_i])
    time_fit = np.linspace(time_change_plot, time_zero_plot, int(time_zero_plot-time_change_plot))
    
    ax = axs[i//3, i%3]
    ax.plot(time_point, data_point, color='black', linewidth=1.5, label='data')
    ax.plot(time_fit, trend_plot*time_fit + base_plot, color='red', linewidth=1.5, label='fit')
    ax.axvline(x=time_change_plot, color='blue', linewidth=1.5, label='CP yr '+str(time_change_plot)+'')
    ax.legend(loc='upper right')
    ax.set_xlabel('Model year')
    ax.set_ylabel('Restoring rate')
    ax.grid(True)
    ax.set_title('(('+str(int(depth[depth_i]))+'m, '+str(int(lon[lon_i]))+'°E), tipping time '+str(int(time_zero_plot))+'')

plt.tight_layout()
#plt.suptitle(f'Restoring rate temperature SAMBA, end year {end_year}, MinThreshold {minthreshold}')
#plt.savefig(directory_figures + 'LAMBDA_QE_fit_to_zero_SAMBA_SALT_cesm_w70_minthreshold1_r_07_9worstfits_1758.pdf')
plt.show()

#%% Determine mean and 90% for final histogram over different ranges of end years (start with 10 years mean)

bins = bins[0,:] #all bins are the same

ci_lower_bin_center = ma.masked_all(len(end_year_change))
ci_upper_bin_center = ma.masked_all(len(end_year_change))
median_bin_center = ma.masked_all(len(end_year_change))
mean_bin_center = ma.masked_all(len(end_year_change))
final_hist = ma.masked_all((len(end_year_change), number_bins))

for end_i in range(1,len(end_year_change)):
    print(end_i)
    #final_hist[end_i,:] = final_hist[end_i,:] = np.nanmean(normalized_hist[0:end_i, :], axis = 0)
    final_hist[end_i,:] = normalized_hist[end_i,:]
    
    # Calculate the cumulative sum of the histogram
    cdf = np.cumsum(final_hist[end_i,:])
    cdf /= cdf[-1]
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    
    # Find the bin where the CDF first exceeds 0.10, 0.5, and 0.90
    ci_lower_bin_index = np.where(cdf >= 0.10)[0][0]
    median_bin_index = np.where(cdf >= 0.5)[0][0]
    ci_upper_bin_index = np.where(cdf >= 0.90)[0][0]
    
    ci_lower_bin_center[end_i] = bin_centers[ci_lower_bin_index]
    median_bin_center[end_i] = bin_centers[median_bin_index]
    ci_upper_bin_center[end_i] = bin_centers[ci_upper_bin_index]
    mean_bin_center[end_i] = np.nansum(bin_centers * final_hist[end_i,:]) / np.sum(final_hist[end_i,:])
    
plt.figure()
plt.title('Summed histogram of normalized PDFs salinity QE')
plt.stairs(final_hist[30,:], bins, fill=True)
plt.grid()
#plt.xlim(2020,2100)

fig, axs = plt.subplots(4,1, figsize=(8,6))

axs[0].plot(end_year_change[:], mean_bin_center)
axs[0].set_title('Mean of estimated tipping times')
axs[0].set_xlabel('Change point end year restriction')
axs[0].set_ylabel('Year')
axs[0].grid('True')

axs[1].plot(end_year_change[:], median_bin_center)
axs[1].set_title('Median of estimated tipping times')
axs[1].set_xlabel('Change point end year restriction')
axs[1].set_ylabel('Year')
axs[1].grid('True')

axs[2].plot(end_year_change[:], ci_lower_bin_center)
axs[2].set_title('Lower 90% confidence interval of estimated tipping times')
axs[2].set_xlabel('Change point end year restriction')
axs[2].set_ylabel('Year')
axs[2].grid('True')

axs[3].plot(end_year_change[:], ci_upper_bin_center)
axs[3].set_title('Upper 90% confidence interval of estimated tipping times')
axs[3].set_xlabel('Change point end year restriction')
axs[3].set_ylabel('Year')
axs[3].grid('True')

#plt.suptitle('Salinity along SAMBA (ORAS5)', y=1.02, fontsize=20)
plt.tight_layout()
#plt.savefig(directory_figures + 'mean_median_CI_salinity_samba_RCP85_range_CP_restriction.pdf')

#%% Average mean, 10 and 90CI for the different CPends

print(f'{np.mean(mean_bin_center)}% mean estimated tipping time for all CPends')
print(f'{np.mean(ci_lower_bin_center)}% mean 10\% CI for all CPends')
print(f'{np.mean(ci_upper_bin_center)}% mean 90\% CI for all CPends')
print(f'{np.median(median_bin_center)}% median 90\% CI for all CPends')

fig, ax = plt.subplots(1,2, figsize=(10,4))

ax[0].plot(end_year_change[:], median_bin_center, 'k--', label='Median')
ax[0].plot(end_year_change[:], mean_bin_center, color='black', label='Mean')
ax[0].fill_between(end_year_change[:], ci_lower_bin_center, ci_upper_bin_center, color='green', alpha=0.2, label='10-90% CI')
ax[0].plot(end_year_change[:], ci_lower_bin_center, color='green', linewidth = 2, alpha=0.8)
ax[0].plot(end_year_change[:], ci_upper_bin_center, color='green', linewidth = 2, alpha=0.8)
ax[0].axhline(y=1758, color='red', linestyle='-', label='CESM tipping point')
#ax[0].fill_between(end_year_change_1500[10::], 1741, 1775, color='red', alpha=0.4)
ax[0].set_title('c) $\\tau_e$ using $\mathrm{RES}^S$ at $34^\circ$S (R<-0.7)')
ax[0].set_xlabel('CP$_\mathrm{end}$ [model year]')
ax[0].set_ylabel('$\\tau_e$ [model year]')
ax[0].grid('True')
ax[0].legend(loc=2)
ax[0].set_ylim(1650, 2200)

ax[1].stairs(final_hist[-1,:], bins, fill=True, color='green', alpha=0.8)
ax[1].grid('True')
ax[1].set_title('d) Average PDF of $\\tau_e$ using $\mathrm{RES}^S$ at $34^\circ$S (R<-0.7)')
ax[1].set_xlabel('$\\tau_e$ [model year]')
ax[1].axvline(1758, color='red')
ax[1].set_xlim(1650, 1955)
ax[1].set_ylim(0,0.02)

plt.tight_layout()
#plt.suptitle('Temperature along SAMBA (quasi-equilibrium)', y=1.02, fontsize=20)
plt.savefig(directory_figures + 'mean_median_CI_salinity_samba_QE_range_CP_restriction_allgridpoints_areaweighted_window70_1500_1700_ALLCP_R07.pdf', bbox_inches='tight')

#%% Normalize to give same weight to all points (for time_zero_positive, so with some R>0)

number_bins = 1000

normalized_hist_positive = ma.masked_all((len(time_zero_positive[0,:,1,1]), number_bins))
bins_positive  = ma.masked_all((len(time_zero_positive[0,:,1,1]), number_bins + 1))
#final_hist = ma.masked_all(len(time_change[0]))

for i in range(len(time_zero_positive[0,:,1,1])):
    
    # Create a mask and use the non-mask values only for the histogram
    #mask = time_zero_good[i,:,:].mask
    
    #Take PDF, because you want every CPend to weigh the same, ongeacht van het aantal punten 
    #normalized_hist[i,:], bins[i,:] = np.histogram(time_zero_good[i,~mask], bins=number_bins, range=(np.min(time_zero_good), np.max(time_zero_good)), density=True, weights=area_grids[~mask])
    data = time_zero_positive[:,i,:,:].flatten()

    normalized_hist_positive[i,:], bins_positive[i,:] = np.histogram(data, bins=number_bins, range=(np.nanmin(time_zero_positive), np.nanmax(time_zero_positive)), density=True, weights=area_grids.flatten())
    print(np.sum(normalized_hist_positive[i,:] * np.diff(bins_positive[i,:]))) #should be very close to 1

#Take sum of normalized histograms to get final histogram (taking the average gives exactly the same shape)    
final_hist_positive = np.nansum(normalized_hist_positive, axis=0)

plt.figure()
plt.title('Summed histogram of normalized PDFs salinity')
plt.stairs(final_hist_positive, bins_positive[-1,:], fill=True)
#plt.xlim(1650, 3000)
#plt.ylim(0,0.1)
plt.grid()

#%% Determine mean and 90% for final histogram over different ranges of end years (start with 10 years mean)

bins_positive = bins_positive[0,:] #all bins are the same

ci_lower_bin_center_positive  = ma.masked_all(len(end_year_change))
ci_upper_bin_center_positive  = ma.masked_all(len(end_year_change))
median_bin_center_positive  = ma.masked_all(len(end_year_change))
mean_bin_center_positive  = ma.masked_all(len(end_year_change))
final_hist_positive  = ma.masked_all((len(end_year_change), number_bins))

for end_i in range(1, len(end_year_change)):
    print(end_i)
    #final_hist_positive[end_i,:] = final_hist_positive[end_i,:] = np.nanmean(normalized_hist_positive[0:end_i, :], axis = 0)
    final_hist_positive[end_i,:] = normalized_hist_positive[end_i, :]
    
    # Calculate the cumulative sum of the histogram
    cdf = np.cumsum(final_hist_positive[end_i,:])
    cdf /= cdf[-1]
    bin_centers = 0.5 * (bins_positive[:-1] + bins_positive[1:])
    
    # Find the bin where the CDF first exceeds 0.10, 0.5, and 0.90
    ci_lower_bin_index = np.where(cdf >= 0.10)[0][0]
    median_bin_index = np.where(cdf >= 0.5)[0][0]
    ci_upper_bin_index = np.where(cdf >= 0.90)[0][0]
    
    ci_lower_bin_center_positive[end_i] = bin_centers[ci_lower_bin_index]
    median_bin_center_positive[end_i] = bin_centers[median_bin_index]
    ci_upper_bin_center_positive[end_i] = bin_centers[ci_upper_bin_index]
    mean_bin_center_positive[end_i] = np.nansum(bin_centers * final_hist_positive[end_i,:]) / np.sum(final_hist_positive[end_i,:])
    
plt.figure()
plt.title('Summed histogram of normalized PDFs salinity QE')
plt.stairs(final_hist_positive[-1,:], bins_positive, fill=True)
plt.grid()
#plt.xlim(2020,2100)

fig, axs = plt.subplots(4,1, figsize=(8,6))

axs[0].plot(end_year_change[:], mean_bin_center_positive)
axs[0].set_title('Mean of estimated tipping times')
axs[0].set_xlabel('Change point end year restriction')
axs[0].set_ylabel('Year')
axs[0].grid('True')

axs[1].plot(end_year_change[:], median_bin_center_positive)
axs[1].set_title('Median of estimated tipping times')
axs[1].set_xlabel('Change point end year restriction')
axs[1].set_ylabel('Year')
axs[1].grid('True')

axs[2].plot(end_year_change[:], ci_lower_bin_center_positive)
axs[2].set_title('Lower 90% confidence interval of estimated tipping times')
axs[2].set_xlabel('Change point end year restriction')
axs[2].set_ylabel('Year')
axs[2].grid('True')

axs[3].plot(end_year_change[:], ci_upper_bin_center_positive)
axs[3].set_title('Upper 90% confidence interval of estimated tipping times')
axs[3].set_xlabel('Change point end year restriction')
axs[3].set_ylabel('Year')
axs[3].grid('True')

#plt.suptitle('Salinity along SAMBA (ORAS5)', y=1.02, fontsize=20)
plt.tight_layout()
#plt.savefig(directory_figures + 'mean_median_CI_salinity_samba_RCP85_range_CP_restriction.pdf')

#%% Average mean, 10 and 90CI for the different CPends

print(f'{np.mean(mean_bin_center)}% mean estimated tipping time for all CPends')
print(f'{np.mean(ci_lower_bin_center)}% mean 10\% CI for all CPends')
print(f'{np.mean(ci_upper_bin_center)}% mean 90\% CI for all CPends')
print(f'{np.median(median_bin_center)}% median 90\% CI for all CPends')

fig, ax = plt.subplots(1,2, figsize=(10,4))

ax[0].plot(end_year_change[:], median_bin_center_positive, 'k--', label='Median')
ax[0].plot(end_year_change[:], mean_bin_center_positive, color='black', label='Mean')
ax[0].fill_between(end_year_change[:], ci_lower_bin_center_positive, ci_upper_bin_center_positive, color='green', alpha=0.2, label='10-90% CI')
ax[0].plot(end_year_change[:], ci_lower_bin_center_positive, color='green', linewidth = 2, alpha=0.8)
ax[0].plot(end_year_change[:], ci_upper_bin_center_positive, color='green', linewidth = 2, alpha=0.8)
ax[0].axhline(y=1758, color='red', linestyle='-', label='CESM tipping point')
#ax[0].fill_between(end_year_change_1500[10::], 1741, 1775, color='red', alpha=0.4)
ax[0].set_title('a) $\\tau_e$ using $\mathrm{RES}^S$ at $34^\circ$S (R<0)')
ax[0].set_xlabel('CP$_\mathrm{end}$ [model year]')
ax[0].set_ylabel('$\\tau_e$ [model year]')
ax[0].grid('True')
ax[0].legend(loc=2)
ax[0].set_ylim(1650, 2200)

ax[1].stairs(final_hist_positive[-1,:], bins_positive, fill=True, color='green', alpha=0.8)
ax[1].grid('True')
ax[1].set_title('b) Average PDF of $\\tau_e$ using $\mathrm{RES}^S$ at $34^\circ$S (R<0)')
ax[1].set_xlabel('$\\tau_e$ [model year]')
ax[1].axvline(1758, color='red')
ax[1].set_xlim(1600, 5000)
ax[1].set_ylim(0,0.02)

plt.tight_layout()
#plt.suptitle('Temperature along SAMBA (quasi-equilibrium)', y=1.02, fontsize=20)
plt.savefig(directory_figures + 'mean_median_CI_salinity_samba_QE_range_CP_restriction_allgridpoints_areaweighted_window70_1500_1700_ALLCP_R0.pdf', bbox_inches='tight')

#%% Normalize to give same weight to all points (for time_zero, so without any restriction)

number_bins = 5000000

normalized_hist_ALL = ma.masked_all((len(time_zero[0,:,1,1]), number_bins))
bins_ALL = ma.masked_all((len(time_zero[0,:,1,1]), number_bins + 1))
#final_hist = ma.masked_all(len(time_change[0]))

 
for i in range(len(time_zero[0,:,0,0])):
    data = time_zero[:,i,:,:].flatten()

    normalized_hist_ALL[i,:], bins_ALL[i,:] = np.histogram(data, bins=number_bins, range=(np.nanmin(time_zero), np.nanmax(time_zero)), density=True, weights=area_grids.flatten())
    print(np.sum(normalized_hist_ALL[i,:] * np.diff(bins_ALL[i,:]))) #should be very close to 1

#Take sum of normalized histograms to get final histogram (taking the average gives exactly the same shape)    
final_hist_ALL = np.nansum(normalized_hist_ALL, axis=0)

plt.figure()
plt.title('Summed histogram of normalized PDFs salinity')
plt.stairs(final_hist_ALL, bins_ALL[-1,:], fill=True)
#plt.xlim(1650, 3000)
#plt.ylim(0,0.1)
plt.grid()

#%% Determine mean and 90% for final histogram over different ranges of end years (start with 10 years mean)

bins = bins_ALL[0,:] #all bins are the same

ci_lower_bin_center_ALL = ma.masked_all(len(end_year_change))
ci_upper_bin_center_ALL = ma.masked_all(len(end_year_change))
median_bin_center_ALL = ma.masked_all(len(end_year_change))
mean_bin_center_ALL = ma.masked_all(len(end_year_change))
final_hist_ALL = ma.masked_all((len(end_year_change), number_bins))

for end_i in range(1,len(end_year_change)):
    #final_hist_ALL[end_i,:] = np.nanmean(normalized_hist_ALL[0:end_i, :], axis = 0)
    final_hist_ALL[end_i,:] = normalized_hist_ALL[end_i,:]
    
    # Calculate the cumulative sum of the histogram
    cdf = np.cumsum(final_hist_ALL[end_i,:])
    cdf /= cdf[-1]
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    
    # Find the bin where the CDF first exceeds 0.10, 0.5, and 0.90
    ci_lower_bin_index = np.where(cdf >= 0.10)[0][0]
    median_bin_index = np.where(cdf >= 0.5)[0][0]
    ci_upper_bin_index = np.where(cdf >= 0.90)[0][0]
    
    ci_lower_bin_center_ALL[end_i] = bin_centers[ci_lower_bin_index]
    median_bin_center_ALL[end_i] = bin_centers[median_bin_index]
    ci_upper_bin_center_ALL[end_i] = bin_centers[ci_upper_bin_index]
    mean_bin_center_ALL[end_i] = np.nansum(bin_centers * final_hist_ALL[end_i,:]) / np.sum(final_hist_ALL[end_i,:])
  
plt.figure()
plt.title('Summed histogram of normalized PDFs salinity QE')
plt.stairs(final_hist_ALL[10,:], bins)#, fill=True)
plt.grid()
#plt.xlim(2020,2100)

fig, axs = plt.subplots(4,1, figsize=(8,6))

axs[0].plot(end_year_change, mean_bin_center_ALL)
axs[0].set_title('Mean of estimated tipping times')
axs[0].set_xlabel('Change point end year restriction')
axs[0].set_ylabel('Year')
axs[0].grid('True')

axs[1].plot(end_year_change, median_bin_center_ALL)
axs[1].set_title('Median of estimated tipping times')
axs[1].set_xlabel('Change point end year restriction')
axs[1].set_ylabel('Year')
axs[1].grid('True')

axs[2].plot(end_year_change, ci_lower_bin_center_ALL)
axs[2].set_title('Lower 90% confidence interval of estimated tipping times')
axs[2].set_xlabel('Change point end year restriction')
axs[2].set_ylabel('Year')
axs[2].grid('True')

axs[3].plot(end_year_change, ci_upper_bin_center_ALL)
axs[3].set_title('Upper 90% confidence interval of estimated tipping times')
axs[3].set_xlabel('Change point end year restriction')
axs[3].set_ylabel('Year')
axs[3].grid('True')

#plt.suptitle('Salinity along SAMBA (ORAS5)', y=1.02, fontsize=20)
plt.tight_layout()
#plt.savefig(directory_figures + 'mean_median_CI_salinity_samba_RCP85_range_CP_restriction.pdf')

plt.figure()
plt.plot(end_year_change[:], ci_upper_bin_center_ALL)
plt.plot(end_year_change[:], ci_lower_bin_center_ALL)
plt.plot(end_year_change[:], mean_bin_center_ALL)
plt.plot(end_year_change[:], median_bin_center_ALL, '--')
plt.legend()


#%%

median_time_zero = np.zeros(len(time_zero[0,:,0,0]))
percentage_tipping_time = np.zeros(len(time_zero[0,:,0,0]))
percentage_greater_than_1e10 = np.zeros(len(time_zero[0,:,0,0]))
total_entries = time_zero[0,0,:,:].size
area_tipping_time = np.zeros(len(time_zero[0,:,0,0]))
area_greater_than_1e10 = np.zeros(len(time_zero[0,:,0,0]))
total_area = np.nansum(area_grids[0, ~lambda_salt[0,:,:].mask])

for i in range(len(time_zero[0,:,0,0])):
    data = time_zero[:,i,:,:].flatten()
    median_time_zero[i] = np.nanmedian(time_zero[:,i,:,:].data)
    percentage_tipping_time[i] = (np.count_nonzero(~np.isnan(time_zero[0,i,:,:])) / total_entries) *100
    percentage_greater_than_1e10[i] = ((np.nansum(time_zero[0,i,:,:] > 1e6) / total_entries)) * 100
    area_tipping_time[i] = (area_grids[0,~np.isnan(time_zero[0,i,:,:])].sum() / total_area) * 100
    area_greater_than_1e10[i] = (area_grids[0,time_zero[0,i,:,:]>1e6].sum() / total_area)  *100


#%%

time_crop_1			  = 1800
factor_time_crop_1	  = 20
median_bin_center_ALL_cropped = np.zeros(len((median_bin_center_ALL)))

for i in range(len(median_bin_center_ALL)):
    if median_bin_center_ALL[i] > time_crop_1:
        median_bin_center_ALL_cropped[i]   = ((median_bin_center_ALL[i] - time_crop_1) / factor_time_crop_1) + time_crop_1
    else:
        median_bin_center_ALL_cropped[i] = median_bin_center_ALL[i]
        
# Create a figure and axis
fig, ax = plt.subplots(2,1, figsize=(8,8))

# Plot percentage_tipping_time on the primary y-axis
color = 'tab:blue'

ax[1].set_ylabel('Percentage area grid cells [%]', fontsize=14)
ax[1].plot(area_tipping_time, color=color, label='All grid cells with CP', linewidth=2)
ax[1].plot(area_greater_than_1e10, '--', color=color, label='R > 0', linewidth=2)
ax[1].tick_params(axis='x', labelsize=14, width=2, length=6, labelcolor='black')
ax[1].tick_params(axis='y', labelsize=14, width=2, length=6)
#ax[1].set_ylim(0,2)
ax[1].set_xlim(0, len(end_year_change))
ax[1].grid()
ax[1].set_title('b) Percentage grid cells used in $\\tau_e$  estimation', fontsize=16)

# Manually set final tick labels
ax[1].set_xticklabels(['1580', '1590', '1600', '1610', '1620', '1630', '1640'])

color = 'tab:red'
ax[0].set_ylabel('Median $\\tau_e$ [model year]', fontsize=14)
ax[0].plot(np.linspace(0,len(end_year_change), np.size(median_bin_center_ALL)), median_bin_center_ALL_cropped, color=color, label='All grid cells with CP', linewidth=2)
ax[0].plot(np.linspace(0,len(end_year_change), np.size(median_bin_center_positive)), median_bin_center_positive, 'o', color=color, label='R<0', linewidth=2)
ax[0].plot(np.linspace(0,len(end_year_change), np.size(median_bin_center)), median_bin_center, '--', color=color, label='R<-0.7', linewidth=2)
#ax2.plot(np.linspace(10,100, np.size(median_bin_center)), median_bin_center_ALL, color=color, label='All slopes', linewidth=2)
ax[0].tick_params(axis='both', labelsize=14, width=2, length=6)
ax[0].set_ylim(1700, 1901)
ax[0].set_xlim(0, len(end_year_change))
ax[1].set_xlabel('CP$_{end}$ [model year]', color='black', fontsize=14)
ax[0].grid()
ax[0].set_xticklabels(['1580', '1590', '1600', '1610', '1620', '1630', '1640'])
ax[0].set_title('a) Median $\\tau_e$ using $\mathrm{RES}^S$ at $34^\circ$S', fontsize=16)

# Rescale certain ticks
labels = ax[0].get_yticks()
for i, val in enumerate(labels):
    if val > time_crop_1:
        labels[i] = ((val - time_crop_1) * factor_time_crop_1) + time_crop_1
labels = labels.astype(int)
ax[0].set_yticklabels(labels)

# Manually set final tick labels
#ax[1].set_yticklabels(['1700', '1750', '1800', '1850', '1900'])

ax[0].hlines(y=1758, xmin=0, xmax=len(end_year_change), color='black', linewidth=2)

# Combine both legends
#lines1, labels1 = ax1.get_legend_handles_labels()
#lines2, labels2 = ax[0,1].get_legend_handles_labels()
ax[0].legend(loc='upper left', fontsize=12)
ax[1].legend(loc='upper left', fontsize=12)

plt.tight_layout()

#plt.suptitle('$\\tau_e$ using $\mathrm{RES}^S$ at $34^\circ$S', fontsize=16)
plt.savefig(directory_figures + 'Median_percentage_tipping_time_ALL_positive_R07_salinity_samba_QE_range_CP_restriction_allgridpoints_areaweighted_window70_1500_1700_ALLCP.pdf', bbox_inches='tight')


plt.show()
