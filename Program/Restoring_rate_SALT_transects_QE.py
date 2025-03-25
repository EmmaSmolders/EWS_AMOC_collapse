#Script determines the restoring rate of a variable in the quasi-equilibrium simulation 

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

#Making pathway to folder with all data
directory_transient	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/TEMP_SALT_ATLANTIC_transects/'
directory		= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_lambda/'
directory_output	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/'

#-----------------------------------------------------------------------------------------

#Functions
def run_fit_a_ar1(x, w):
    
  n = x.shape[0] #n is length of the input array x
  xs = np.zeros_like(x) #xs will store the autoregressive coefficients (same shape as x)
  xs_err = np.zeros_like(x)
  
  #Setting first (window/2) and last (window/2) entries of array to nan (because of sliding window)
  for i in range(w // 2):
     xs[i] = np.nan

  for i in range(n - w // 2, n):
     xs[i] = np.nan
     
  #Iteration over entries in between the nans at both ends in which the autoregressive modeling is performed
  for i in range(w // 2, n - w // 2):
     xw = x[i - w // 2 : i + w // 2 + 1] #extracts a window of size w centred around current index i
     xw = xw - xw.mean() #centers the window by subtraction of mean

     p0, p1 = np.polyfit(np.arange(xw.shape[0]), xw, 1) #linear fit to centered window
     xw = xw - p0 * np.arange(xw.shape[0]) - p1 #removes the linear trend from the centered window

     dxw = xw[1:] - xw[:-1] #Calculate difference of detrended window? dx/dt?

     xw = sm.add_constant(xw) #add constant to detrended window

     model = sm.GLSAR(dxw, xw[:-1], rho=1) #Create a GLSAR (Generalized Least Squares AutoRegresive) model with first differences as the dependent variable and the detrended window as the independent varaible. rho=1 specifies an AR model of order 1
     results = model.iterative_fit(maxiter=10) #Iteratively fits the GLSAR model with a maximum of 10 iterations

     a = results.params[1] #retrieves the AR coefficient from the fitted model

     xs[i] = a
     xs_err[i] = results.bse[1] #standard error

  return xs

#-------------------------------------------------------------------------------------------------------------
#Determine restoring rate of transient simulation
fh = netcdf.Dataset(directory_transient + 'TEMP_SALT_34S_year_0-1750.nc','r')
#fh = netcdf.Dataset(directory_transient + 'TEMP_SALT_26N_year_0-1750.nc','r')
#fh = netcdf.Dataset(directory_transient + 'TEMP_SALT_60N_year_1500-1750.nc', 'r')

#Get variables (yearly averaged)
time  = fh.variables['time'][:]     #time [year]
dens  = fh.variables['PD'][:]       #potential density [kg/m^3]
salt  = fh.variables['SALT'][:]     #salinity [g/kg]
temp  = fh.variables['TEMP'][:]     #temperature [deg C]
depth = fh.variables['depth'][:]    #depth [m]
lon   = fh.variables['lon'][:]      
lat   = fh.variables['lat'][:]

fh.close()

#-------------------------------------------------------------------------------------------------------------

print(time)

#Sliding window (years)
window = 70 

#Variables you want to use
data = salt

# Choose the range of end times for which you want to determine lambda
end_times = range(1700, 1703, 2)
#end_times = range(200, 203, 2)

# Empty dictionary to store the lambda values for each end time
lambda_values = {}

#time_begin = 1500
time_begin = 0

plt.figure()
plt.plot(data[time_begin:end_times[0], 0,0])
plt.show()

if ma.isMaskedArray(data[0:end_times[0],0,0]):
    print("Array is a masked array.")
else:
    print("Array is not a masked array.")

# Determine lambda in gridpoints for each end time
for time_end in end_times:
    print(time_end)
    # Empty array for the current end time
    lambda_transect = ma.masked_all((len(time[time_begin:time_end]), len(depth), len(lon)))

    for depth_i in range(len(depth)):
        print(depth_i)
        for lon_i in range(len(lon)):
            print(lon_i)

	    #Skip if land
            if data[time_begin:time_end,depth_i,lon_i] is ma.masked:
                continue

            if ma.count(data[time_begin:time_end,depth_i,lon_i]) == 0:
                continue

            print(data[time_begin:time_end,depth_i,lon_i])
            print(ma.count(data[time_begin:time_end,depth_i,lon_i]))	    

            lambda_transect[:,depth_i,lon_i] = run_fit_a_ar1(data[time_begin:time_end,depth_i,lon_i], window)

    # Store the lambda values for the current end time in the dictionary
    lambda_values[time_end] = lambda_transect

#-------------------------------------------------------------------------------------------------------------
#Store data
with netcdf.Dataset(directory_output + 'lambda_SALT_transient_1500_endtimes_1700_1702_spacing_2_window_'+str(window)+'_SAMBA_all_gridpoints.nc', 'w') as f:
    # Create the dimensions
    f.createDimension('depth', len(depth))
    f.createDimension('lon', len(lon))
    f.createDimension('end_time', len(end_times))

    # Create the variables
    depths_var = f.createVariable('depth', float, ('depth',))
    lons_var = f.createVariable('lon', float, ('lon',))
    end_times_var = f.createVariable('end_time', int, ('end_time',))

    # Assign the data to the variables
    depths_var[:] = depth
    lons_var[:] = lon
    end_times_var[:] = list(end_times)

    # Create a new dimension and variable for each end_time
    for i, end_time in enumerate(end_times):
        time_dim = f.createDimension(f'time_{end_time}', len(time[time_begin:end_time]))
        time_var = f.createVariable(f'time_{end_time}', float, (f'time_{end_time}',))
        time_var[:] = time[time_begin:end_time]

        lambda_values_var = f.createVariable(f'lambda_values_{end_time}', float, (f'time_{end_time}', 'depth', 'lon'))
        lambda_values_var[:, :, :] = lambda_values[end_time]







