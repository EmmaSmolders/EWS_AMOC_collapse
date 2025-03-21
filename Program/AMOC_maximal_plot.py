#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 16:06:58 2024

@author: 6008399

Time series of AMOC maximum at certain latitude of CESM transient and branched simulations (Figure 1 op paper)

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
directory	= '../../../Data/AMOC/'

#%% maximal AMOC strength under 500m

fh_transient = netcdf.Dataset(directory + 'AMOC_max_year_0-2200.nc','r')
fh_branch1   = netcdf.Dataset(directory + 'AMOC_max_year_600-1200_branch_600.nc','r')
fh_branch2   = netcdf.Dataset(directory + 'AMOC_max_year_1500-2100_branch_1500.nc','r')
fh_branch3   = netcdf.Dataset(directory + 'AMOC_max_year_2900-3051_branch_1650.nc','r')

time_transient      = fh_transient.variables['time'][:]     #time in model years
time_branch1        = fh_branch1.variables['time'][:]       #time in model years
time_branch2        = fh_branch2.variables['time'][:]       #time in model years
time_branch3        = fh_branch3.variables['time'][:]       #time in model years

lat                 = fh_transient.variables['lat'][:]      #latitude [degN]

depth_max_transient = fh_transient.variables['depth'][:]    #depth of AMOC maximum [m]
depth_max_branch1   = fh_branch1.variables['depth'][:]      #depth of AMOC maximum [m]
depth_max_branch2   = fh_branch2.variables['depth'][:]      #depth of AMOC maximum [m]
depth_max_branch3   = fh_branch3.variables['depth'][:]      #depth of AMOC maximum [m]

AMOC_max_transient  = fh_transient.variables['AMOC_max'][:] #maximum of AMOC [Sv]
AMOC_max_branch1    = fh_branch1.variables['AMOC_max'][:]   #maximum of AMOC [Sv]
AMOC_max_branch2    = fh_branch2.variables['AMOC_max'][:]   #maximum of AMOC [Sv]
AMOC_max_branch3    = fh_branch3.variables['AMOC_max'][:]   #maximum of AMOC [Sv]

fh_transient.close()
fh_branch1.close()
fh_branch2.close()
fh_branch3.close()

#%% AMOC structure last 50 years of branch 1 and 2

fh = netcdf.Dataset(directory + 'AMOC_structure_year_1050-1100_branch_600.nc', 'r')

depth	= fh.variables['depth'][:] 	 / 1000.0	
lat	    = fh.variables['lat'][:] 		
AMOC_1	= np.mean(fh.variables['AMOC'][:], axis = 0)	

fh.close()

fh = netcdf.Dataset(directory + 'AMOC_structure_year_1950-2000_branch_1500.nc', 'r')

AMOC_2	= np.mean(fh.variables['AMOC'][:], axis = 0)	

fh.close()

#%% FOV 

fh_fov_60N_transient = netcdf.Dataset(directory + 'FOV_60N_year_0_2200.nc','r')
fh_fov_26N_transient = netcdf.Dataset(directory + 'FOV_26N_year_0_2200.nc','r')
fh_fov_34S_transient = netcdf.Dataset(directory + 'FOV_34S_year_0_2200.nc','r')

fh_fov_60N_branch1 = netcdf.Dataset(directory + 'FOV_60N_year_600_1100_branch_600.nc','r')
fh_fov_26N_branch1 = netcdf.Dataset(directory + 'FOV_26N_year_600_1100_branch_600.nc','r')
fh_fov_34S_branch1 = netcdf.Dataset(directory + 'FOV_34S_year_600_1100_branch_600.nc','r')

fh_fov_60N_branch2 = netcdf.Dataset(directory + 'FOV_60N_year_1500_2000_branch_1500.nc','r')
fh_fov_26N_branch2 = netcdf.Dataset(directory + 'FOV_26N_year_1500_2000_branch_1500.nc','r')
fh_fov_34S_branch2 = netcdf.Dataset(directory + 'FOV_34S_year_1500_2000_branch_1500.nc','r')

fh_fov_60N_branch3 = netcdf.Dataset(directory + 'FOV_60N_year_2900_3051_branch_1650.nc','r')
fh_fov_26N_branch3 = netcdf.Dataset(directory + 'FOV_26N_year_2900_3051_branch_1650.nc','r')
fh_fov_34S_branch3 = netcdf.Dataset(directory + 'FOV_34S_year_2900_3051_branch_1650.nc','r')

FOV_60N_transient  = fh_fov_60N_transient.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_26N_transient  = fh_fov_26N_transient.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_34S_transient  = fh_fov_34S_transient.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])

FOV_60N_branch1  = fh_fov_60N_branch1.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_26N_branch1  = fh_fov_26N_branch1.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_34S_branch1  = fh_fov_34S_branch1.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])

FOV_60N_branch2  = fh_fov_60N_branch2.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_26N_branch2  = fh_fov_26N_branch2.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_34S_branch2  = fh_fov_34S_branch2.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])

FOV_60N_branch3  = fh_fov_60N_branch3.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_26N_branch3  = fh_fov_26N_branch3.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])
FOV_34S_branch3  = fh_fov_34S_branch3.variables['F_OV'][:] #Freshwater transport (overturning component [Sv])


#%% Figure 1 of revised paper (as function of Fh and including E3)

fig, axs = plt.subplots(2, 3, figsize=(12, 6))

axs[0,0].plot(time_transient*0.0003, AMOC_max_transient[:, lat_idx_samba], color = 'black', label='Transient', linewidth=0.5)

axs[0,0].plot(time_branch1[0]*0.0003, np.mean(AMOC_max_branch1[350:500, lat_idx_samba]), 'o', color = 'orchid', label='E1', markersize=10)
axs[0,0].vlines(x = time_branch1[0]*0.0003, ymin=np.min(AMOC_max_branch1[350:500, lat_idx_samba]), ymax= np.max(AMOC_max_branch1[350:500, lat_idx_samba]), color='orchid', linewidth=2)
axs[0,0].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.min(AMOC_max_branch1[350:500, lat_idx_samba]), color='orchid', linewidth=2)
axs[0,0].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.max(AMOC_max_branch1[350:500, lat_idx_samba]), color='orchid', linewidth=2)

axs[0,0].plot(time_branch2[0]*0.0003, np.mean(AMOC_max_branch2[350:500, lat_idx_samba]), 'o', color = 'red', label='E2', markersize=10)
axs[0,0].vlines(x = time_branch2[0]*0.0003, ymin=np.min(AMOC_max_branch2[350:500, lat_idx_samba]), ymax= np.max(AMOC_max_branch2[350:500, lat_idx_samba]), color='red', linewidth=2)
axs[0,0].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.min(AMOC_max_branch2[350:500, lat_idx_samba]), color='red', linewidth=2)
axs[0,0].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.max(AMOC_max_branch2[350:500, lat_idx_samba]), color='red', linewidth=2)

axs[0,0].plot(time_branch3[0]*0.0003 - 1250*0.0003, np.mean(AMOC_max_branch3[:, lat_idx_samba]), 'o', color = 'darkturquoise', label='E3', markersize=10)
axs[0,0].vlines(x = time_branch3[0]*0.0003 - 1250*0.0003, ymin=np.min(AMOC_max_branch3[:, lat_idx_samba]), ymax= np.max(AMOC_max_branch3[:, lat_idx_samba]), color='darkturquoise', linewidth=2)
axs[0,0].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.min(AMOC_max_branch3[:, lat_idx_samba]), color='darkturquoise', linewidth=2)
axs[0,0].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.max(AMOC_max_branch3[:, lat_idx_samba]), color='darkturquoise', linewidth=2)

axs[0,0].set_ylim(-1,21)
axs[0,0].set_xlim(0,0.66)
axs[0,0].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

axs[0,0].set_ylabel('Volume transport [Sv]')
#axs[0,0].set_xlabel('Time [model years]')
axs[0,0].set_title('a) AMOC strength at 34$^\circ$S')
axs[0,0].grid()

graph_1         = axs[0,0].plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Transient')
graph_2         = axs[0,0].plot([-100, -100], [-100, -100], 'o', color ='orchid', label = 'E1')
graph_3         = axs[0,0].plot([-100, -100], [-100, -100], 'o', color ='red',  label = 'E2')
graph_4         = axs[0,0].plot([-100, -100], [-100, -100], 'o', color = 'darkturquoise', label = 'E3')

graphs          = graph_1 + graph_2 + graph_3 + graph_4
legend_labels   = [l.get_label() for l in graphs]
legend_1        = axs[0,0].legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)


ax2 	= fig.add_axes([0.75, 0.145, 0.212, 0.150])

CS	 = ax2.contourf(lat, depth, AMOC_2 - AMOC_1, levels = np.arange(-4, 4.01, 0.5), extend = 'both', cmap = 'RdBu_r')
cbar	= colorbar(CS, ticks = np.arange(-4, 4.01, 2))
cbar.set_label('AMOC diff (Sv)', fontsize = 9)

ax2.set_xlim(-30, 62)
ax2.set_ylim(5, 0)
ax2.set_yticks([0, 2, 4])
ax2.set_ylabel('Depth (km)', fontsize = 8)	

ax2.set_xticks(np.arange(-20, 60.1, 20))
ax2.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'], fontsize = 8)

CS_1	= ax2.contour(lat, depth, AMOC_1, levels = [14], colors = 'k', linewidths = 1)
CS_1	= ax2.contour(lat, depth, AMOC_1, levels = [12], colors = 'r', linewidths = 1)
CS_1	= ax2.contour(lat, depth, AMOC_1, levels = [9], colors = 'b', linewidths = 1)
CS_1	= ax2.contour(lat, depth, AMOC_1, levels = [-1], colors = 'k', linewidths = 1)

ax2.set_title('E2 minus E1', fontsize = 9)

ax3 	= fig.add_axes([0.01, 0.11, 0.25, 0.25], projection = ccrs.Orthographic(-30, 10))

ax3.coastlines(resolution='110m')
ax3.gridlines()
ax3.add_feature(cfeature.LAND, zorder=10)
ax3.set_global()


lon1     = np.arange(0, 361)
lat1     = np.arange(-90, 91)
field   = np.ones((len(lat1), len(lon1))) * -0.35
CS      = ax3.contourf(lon1, lat1, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon2     = np.arange(-100, -5)
lat2     = np.arange(20, 43)
field   = np.ones((len(lat2), len(lon2))) * 0.35
CS      = ax3.contourf(lon2, lat2, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon3     = np.arange(-100, 3)
lat3     = np.arange(42, 51)
field   = np.ones((len(lat3), len(lon3))) * 0.35
CS      = ax3.contourf(lon3, lat3, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

ax3.text(320, 38, '$+F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())
ax3.text(340, -10, '$-F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())

x_1	= np.arange(-65, 7.1, 0.1)
y_1	= np.zeros(len(x_1)) + 60.0
y_2	= np.arange(58, 62.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(58, 62.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax3.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

x_1	= np.arange(-81, -9.99, 0.1)
y_1	= np.zeros(len(x_1)) + 26.0
y_2	= np.arange(24, 28.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(24, 28.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax3.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

x_1	= np.arange(-60, 20.01, 0.1)
y_1	= np.zeros(len(x_1)) - 34
y_2	= np.arange(-37, -30.99, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(-37, -30.99, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax3.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

axs[0,1].plot(time_transient*0.0003, AMOC_max_transient[:, lat_idx_rapid], color = 'black', label='Transient', linewidth=0.5)

axs[0,1].plot(time_branch1[0]*0.0003, np.mean(AMOC_max_branch1[350:500, lat_idx_rapid]), 'o', color = 'orchid', label='E1', markersize=10)
axs[0,1].vlines(x = time_branch1[0]*0.0003, ymin=np.min(AMOC_max_branch1[350:500, lat_idx_rapid]), ymax= np.max(AMOC_max_branch1[350:500, lat_idx_rapid]), color='orchid', linewidth=2)
axs[0,1].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.min(AMOC_max_branch1[350:500, lat_idx_rapid]), color='orchid', linewidth=2)
axs[0,1].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.max(AMOC_max_branch1[350:500, lat_idx_rapid]), color='orchid', linewidth=2)

axs[0,1].plot(time_branch2[0]*0.0003, np.mean(AMOC_max_branch2[350:500, lat_idx_rapid]), 'o', color = 'red', label='E2', markersize=10)
axs[0,1].vlines(x = time_branch2[0]*0.0003, ymin=np.min(AMOC_max_branch2[350:500, lat_idx_rapid]), ymax= np.max(AMOC_max_branch2[350:500, lat_idx_rapid]), color='red', linewidth=2)
axs[0,1].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.min(AMOC_max_branch2[350:500, lat_idx_rapid]), color='red', linewidth=2)
axs[0,1].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.max(AMOC_max_branch2[350:500, lat_idx_rapid]), color='red', linewidth=2)

axs[0,1].plot(time_branch3[0]*0.0003 - 1250*0.0003, np.mean(AMOC_max_branch3[:, lat_idx_rapid]), 'o', color = 'darkturquoise', label='E3', markersize=10)
axs[0,1].vlines(x = time_branch3[0]*0.0003 - 1250*0.0003, ymin=np.min(AMOC_max_branch3[:, lat_idx_rapid]), ymax= np.max(AMOC_max_branch3[:, lat_idx_rapid]), color='darkturquoise', linewidth=2)
axs[0,1].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.min(AMOC_max_branch3[:, lat_idx_rapid]), color='darkturquoise', linewidth=2)
axs[0,1].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.max(AMOC_max_branch3[:, lat_idx_rapid]), color='darkturquoise', linewidth=2)

axs[0,1].set_ylim(-1,21)
axs[0,1].set_xlim(0,0.66)
axs[0,1].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
#axs[0,1].set_xlabel('Time [model years]') 
axs[0,1].set_title('b) AMOC strength at 26$^\circ$N')
axs[0,1].grid()
#axs[1].legend()

axs[0,2].plot(time_transient*0.0003, AMOC_max_transient[:, lat_idx_60N], color = 'black', label='Transient', linewidth=0.5)

axs[0,2].plot(time_branch1[0]*0.0003, np.mean(AMOC_max_branch1[350:500, lat_idx_60N]), 'o', color = 'orchid', label='E1', markersize=10)
axs[0,2].vlines(x = time_branch1[0]*0.0003, ymin=np.min(AMOC_max_branch1[350:500, lat_idx_60N]), ymax= np.max(AMOC_max_branch1[350:500, lat_idx_60N]), color='orchid', linewidth=2)
axs[0,2].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.min(AMOC_max_branch1[350:500, lat_idx_60N]), color='orchid', linewidth=2)
axs[0,2].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.max(AMOC_max_branch1[350:500, lat_idx_60N]), color='orchid', linewidth=2)

axs[0,2].plot(time_branch2[0]*0.0003, np.mean(AMOC_max_branch2[350:500, lat_idx_60N]), 'o', color = 'red', label='E2', markersize=10)
axs[0,2].vlines(x = time_branch2[0]*0.0003, ymin=np.min(AMOC_max_branch2[350:500, lat_idx_60N]), ymax= np.max(AMOC_max_branch2[350:500, lat_idx_60N]), color='red', linewidth=2)
axs[0,2].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.min(AMOC_max_branch2[350:500, lat_idx_60N]), color='red', linewidth=2)
axs[0,2].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.max(AMOC_max_branch2[350:500, lat_idx_60N]), color='red', linewidth=2)

axs[0,2].plot(time_branch3[0]*0.0003 - 1250*0.0003, np.mean(AMOC_max_branch3[:, lat_idx_60N]), 'o', color = 'darkturquoise', label='E3', markersize=10)
axs[0,2].vlines(x = time_branch3[0]*0.0003 - 1250*0.0003, ymin=np.min(AMOC_max_branch3[:, lat_idx_60N]), ymax= np.max(AMOC_max_branch3[:, lat_idx_60N]), color='darkturquoise', linewidth=2)
axs[0,2].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.min(AMOC_max_branch3[:, lat_idx_60N]), color='darkturquoise', linewidth=2)
axs[0,2].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.max(AMOC_max_branch3[:, lat_idx_60N]), color='darkturquoise', linewidth=2)

axs[0,2].set_ylim(-1,21)
axs[0,2].set_xlim(0,0.66)
axs[0,2].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

axs[0,2].set_title('c) AMOC strength at 60$^\circ$N')
axs[0,2].grid()
#axs[2].legend()

axs[1,0].plot(time_transient*0.0003, FOV_34S_transient, color = 'black', label='Transient', linewidth=0.5)

axs[1,0].plot(time_branch1[0]*0.0003, np.mean(FOV_34S_branch1[350:500]), 'o', color = 'orchid', label='E1', markersize=10)
axs[1,0].vlines(x = time_branch1[0]*0.0003, ymin=np.min(FOV_34S_branch1[350:500]), ymax= np.max(FOV_34S_branch1[350:500]), color='orchid', linewidth=2)
axs[1,0].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.min(FOV_34S_branch1[350:500]), color='orchid', linewidth=2)
axs[1,0].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.max(FOV_34S_branch1[350:500]), color='orchid', linewidth=2)

axs[1,0].plot(time_branch2[0]*0.0003, np.mean(FOV_34S_branch2[350:500]), 'o', color = 'red', label='E2', markersize=10)
axs[1,0].vlines(x = time_branch2[0]*0.0003, ymin=np.min(FOV_34S_branch2[350:500]), ymax= np.max(FOV_34S_branch2[350:500]), color='red', linewidth=2)
axs[1,0].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.min(FOV_34S_branch2[350:500]), color='red', linewidth=2)
axs[1,0].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.max(FOV_34S_branch2[350:500]), color='red', linewidth=2)

axs[1,0].plot(time_branch3[0]*0.0003 - 1250*0.0003, np.mean(FOV_34S_branch3), 'o', color = 'darkturquoise', label='E3', markersize=10)
axs[1,0].vlines(x = time_branch3[0]*0.0003 - 1250*0.0003, ymin=np.min(FOV_34S_branch3), ymax= np.max(FOV_34S_branch3), color='darkturquoise', linewidth=2)
axs[1,0].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.min(FOV_34S_branch3), color='darkturquoise', linewidth=2)
axs[1,0].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.max(FOV_34S_branch3), color='darkturquoise', linewidth=2)

axs[1,0].set_xlim(0,0.66)
axs[1,0].set_ylim(-0.75,0.3)
axs[1,0].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])


axs[1,0].set_ylabel('Freshwater transport [Sv]')
axs[1,0].set_xlabel('Freshwater flux forcing F$_H$ [Sv]')
axs[1,0].set_title('d) $F_{\mathrm{ov}}$ at 34$^\circ$S')
axs[1,0].grid()
#axs[1,0].legend(loc=3)

axs[1,1].plot(time_transient*0.0003, FOV_26N_transient, color = 'black', label='Transient', linewidth=0.5)

axs[1,1].plot(time_branch1[0]*0.0003, np.mean(FOV_26N_branch1[350:500]), 'o', color = 'orchid', label='E1', markersize=10)
axs[1,1].vlines(x = time_branch1[0]*0.0003, ymin=np.min(FOV_26N_branch1[350:500]), ymax= np.max(FOV_26N_branch1[350:500]), color='orchid', linewidth=2)
axs[1,1].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.min(FOV_26N_branch1[350:500]), color='orchid', linewidth=2)
axs[1,1].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.max(FOV_26N_branch1[350:500]), color='orchid', linewidth=2)

axs[1,1].plot(time_branch2[0]*0.0003, np.mean(FOV_26N_branch2[350:500]), 'o', color = 'red', label='E2', markersize=10)
axs[1,1].vlines(x = time_branch2[0]*0.0003, ymin=np.min(FOV_26N_branch2[350:500]), ymax= np.max(FOV_26N_branch2[350:500]), color='red', linewidth=2)
axs[1,1].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.min(FOV_26N_branch2[350:500]), color='red', linewidth=2)
axs[1,1].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.max(FOV_26N_branch2[350:500]), color='red', linewidth=2)

axs[1,1].plot(time_branch3[0]*0.0003 - 1250*0.0003, np.mean(FOV_26N_branch3), 'o', color = 'darkturquoise', label='E3', markersize=10)
axs[1,1].vlines(x = time_branch3[0]*0.0003 - 1250*0.0003, ymin=np.min(FOV_26N_branch3), ymax= np.max(FOV_26N_branch3), color='darkturquoise', linewidth=2)
axs[1,1].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.min(FOV_26N_branch3), color='darkturquoise', linewidth=2)
axs[1,1].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.max(FOV_26N_branch3), color='darkturquoise', linewidth=2)

axs[1,1].set_xlim(0,0.66)
axs[1,1].set_ylim(-0.75,0.3)
axs[1,1].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
#axs[1,1].set_ylabel('$F_{OV}$ [Sv]')
axs[1,1].set_xlabel('Freshwater flux forcing F$_H$ [Sv]')
axs[1,1].set_title('e) $F_{\mathrm{ov}}$ at 26$^\circ$N')
axs[1,1].set_ylim(-0.75,0.3)
axs[1,1].grid()

axs[1,2].plot(time_transient*0.0003, FOV_60N_transient, color = 'black', label='Transient', linewidth=0.5)

axs[1,2].plot(time_branch1[0]*0.0003, np.mean(FOV_60N_branch1[350:500]), 'o', color = 'orchid', label='E1', markersize=10)
axs[1,2].vlines(x = time_branch1[0]*0.0003, ymin=np.min(FOV_60N_branch1[350:500]), ymax= np.max(FOV_60N_branch1[350:500]), color='orchid', linewidth=2)
axs[1,2].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.min(FOV_60N_branch1[350:500]), color='orchid', linewidth=2)
axs[1,2].hlines(xmin = time_branch1[0]*0.0003 - 20*0.0003, xmax = time_branch1[0]*0.0003 + 20*0.0003, y=np.max(FOV_60N_branch1[350:500]), color='orchid', linewidth=2)

axs[1,2].plot(time_branch2[0]*0.0003, np.mean(FOV_60N_branch2[350:500]), 'o', color = 'red', label='E2', markersize=10)
axs[1,2].vlines(x = time_branch2[0]*0.0003, ymin=np.min(FOV_60N_branch2[350:500]), ymax= np.max(FOV_60N_branch2[350:500]), color='red', linewidth=2)
axs[1,2].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.min(FOV_60N_branch2[350:500]), color='red', linewidth=2)
axs[1,2].hlines(xmin = time_branch2[0]*0.0003 - 20*0.0003, xmax = time_branch2[0]*0.0003 + 20*0.0003, y=np.max(FOV_60N_branch2[350:500]), color='red', linewidth=2)

axs[1,2].plot(time_branch3[0]*0.0003 - 1250*0.0003, np.mean(FOV_60N_branch3), 'o', color = 'darkturquoise', label='E3', markersize=10)
axs[1,2].vlines(x = time_branch3[0]*0.0003 - 1250*0.0003, ymin=np.min(FOV_60N_branch3), ymax= np.max(FOV_60N_branch3), color='darkturquoise', linewidth=2)
axs[1,2].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.min(FOV_60N_branch3), color='darkturquoise', linewidth=2)
axs[1,2].hlines(xmin = (time_branch3[0] - 1250 - 20)*0.0003, xmax = (time_branch3[0] - 1250 + 20)*0.0003, y=np.max(FOV_60N_branch3), color='darkturquoise', linewidth=2)

axs[1,2].set_xlim(0,0.66)
axs[1,2].set_ylim(-0.75,0.3)
axs[1,2].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
axs[1,2].set_xlabel('Freshwater flux forcing F$_H$ [Sv]')
axs[1,2].set_title('f) $F_{\mathrm{ov}}$ at 60$^\circ$N')
axs[1,2].set_ylim(-0.75,0.3)
axs[1,2].grid()

plt.tight_layout()
plt.show()
