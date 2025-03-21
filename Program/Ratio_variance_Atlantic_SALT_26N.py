#Script which determines single estimator of variance of salinity branch 1 and 2 using model years 350 to 500

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
#import statsmodels.api as sm

#Making pathway to folder with all data
directory_data  = '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_surface/'
directory	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_variance/'
directory_data_1650  	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/TEMP_SALT_ATLANTIC_transects/'

def find_nearest(array, value):
    #array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def AtlanticRegion(lon, lat, data):

	#lon_1 = -61 #osnap
	#lon_2 = 21  #osnap
	
	lon_1 = -80 #rapid
	lon_2 = -12 #rapid

	lon_1_idx = find_nearest(lon[0,:], lon_1)
	lon_2_idx = find_nearest(lon[0,:], lon_2)

	lon = lon[:,lon_1_idx - 1:lon_2_idx + 1]
	lat = lat[:,lon_1_idx - 1:lon_2_idx + 1]
	data = data[:,:,:,lon_1_idx - 1:lon_2_idx + 1]

	return lon, lat, data

def ReadinData(filename, transect = False):

	fh = netcdf.Dataset(filename, 'r')

	depth   	= fh.variables['depth'][:] 	 		#Depth (m)
	lon 		= fh.variables['lon'][:]			#Longitude
	lat 		= fh.variables['lat'][:]			#Latitude 
	time		= fh.variables['time'][:]			#Starting year of window (time = time_year[:len(time_year) - window +1])
	salt		= fh.variables['SALT'][:,:,:]		 	#Salinity 

	fh.close()
	
	if transect == False:
		lon, lat, salt = AtlanticRegion(lon, lat, salt)

	return lon, lat, depth, time, salt

def TrendRemover(time, data):
	"""Removes linear trend"""

	time_fit	= np.arange(1, len(time) + 1)

	#Determine the detrended time series
	trend, base 	= polyfit(time_fit, data, 1)
	data_res	= data - ((trend * time_fit) + base)

	return data_res

#----------------------------------------------------------------------------------------------------------------

#Yearly averaged data from month 1 to 12
month_start	= 1
month_end	= 12

#FILL HERE IN WHICH BRANCHES YOU WANT TO USE (E1 = 600, E2 = 1500, E3 = 1650)
branch1 = 1500
branch2 = 1650

#Which index?
rapid_idx = 270 #(26N)

#Read in data of branch 1 and 2 (yearly averged, only last 150 years). I do it in this way for 26N, because the data is not on Snellius anymore and I haven't saved the data per transect (only for 34S and osnap)
if branch1 == 600:
	lon, lat, depth, time1, salt1		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y600_950-1000.nc', False)
	lon, lat, depth, time2, salt2		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y600_1001-1050.nc', False)
	lon, lat, depth, time3, salt3		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y600_1051-1100.nc', False)
	
	# Concatenate arrays and delete
	time_branch1 = np.concatenate((time1, time2, time3[0:-1]), axis=0) #Take only 150 years, instead of 151 (even number for Fourier_surrogates function)
	salt_branch1 = np.concatenate((salt1, salt2, salt3[0:-1]), axis=0)
	
	#Select transect
	lat1 = lat[rapid_idx, :]
	lon1 = lon[rapid_idx, :]
	print(lat1)
	print(lon1)
	salt_branch1 = salt_branch1[:,:,rapid_idx,:]

	del(salt1, salt2, salt3)

if branch1 == 1500:
	lon, lat, depth, time4, salt4		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1850-1899.nc', False)
	lon, lat, depth, time5, salt5		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1900-1950.nc', False)
	lon, lat, depth, time6, salt6		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1951-2000.nc', False)

	time_branch1 = np.concatenate((time4, time5, time6[0:-1]), axis=0)
	salt_branch1 = np.concatenate((salt4, salt5, salt6[0:-1]), axis=0)
	
	#Select transect
	lat1 = lat[rapid_idx, :]
	lon1 = lon[rapid_idx, :]
	print(lat1)
	print(lon1)
	salt_branch1 = salt_branch1[:,:,rapid_idx,:]

	del(salt4, salt5, salt6)
	
if branch2 == 1650:
	lon, lat, depth, time_branch2, salt_branch2		= ReadinData(directory_data_1650 + 'TEMP_SALT_26N_year_branched_off_y1650_2900-3050.nc', True)
	print(lat)
	print(lon)

if not np.allclose(lat, lat1, atol=1e-10): #apparently there can be a differnce up to e-14?
    print('Latitude array is not the same')
    #sys.exit()
	
if not np.array_equal(lon, lon1):
    print('Latitude array is not the same')
    #sys.exit()
	
if np.shape(salt_branch2) != np.shape(salt_branch1):
	print('Length of array is not the same, and should be equal to 150')
	#sys.exit()
	
print(lat)
print(lon)
print(np.shape(salt_branch2))
print(np.shape(salt_branch1))

#sys.exit()

#Empty arrays
ratio_var_salt = ma.masked_all((len(depth), len(lon)))
var_salt1 = ma.masked_all((len(depth), len(lon)))
var_salt2 = ma.masked_all((len(depth), len(lon)))
stat = ma.masked_all((len(depth), len(lon)))
p = ma.masked_all((len(depth), len(lon)))

#Determine single estimator (VAR_E2/VAR_E1) for each grid point
for depth_i in range(len(depth)):
	print(depth_i)
	for lon_i in range(len(lon)):

		if salt_branch2[0, depth_i, lon_i] is ma.masked:
				#If masked elements (=land), skip
				continue
		
		#Variance
		var_salt1[depth_i, lon_i] = np.var(TrendRemover(time_branch1, salt_branch1[:, depth_i, lon_i]))
		var_salt2[depth_i, lon_i] = np.var(TrendRemover(time_branch2, salt_branch2[:, depth_i, lon_i]))

		#Levene's test to test statistical difference between E1 and E2
		stat[depth_i, lon_i], p[depth_i, lon_i] = stats.levene(TrendRemover(time_branch1, salt_branch1[:, depth_i, lon_i]), TrendRemover(time_branch2, salt_branch2[:, depth_i, lon_i]))
		
		#Single estimator
		ratio_var_salt[depth_i, lon_i] = var_salt2[depth_i, lon_i]/var_salt1[depth_i, lon_i]

	#-----------------------------------------------------------------------------------------

	fh = netcdf.Dataset(directory+'Atlantic_single_estimator_variance_SALT_month_'+str(month_start)+'-'+str(month_end)+'_26N_branch_'+str(branch2)+'_'+str(branch1)+'.nc', 'w')

	fh.createDimension('depth', len(depth))
	fh.createDimension('nlat', len(lat))
	fh.createDimension('nlon', len(lon))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('nlat', float, ('nlat'), zlib=True)
	#fh.createVariable('nlon', float, ('nlon'), zlib=True)
	fh.createVariable('lat', float, ('nlat'), zlib=True)
	fh.createVariable('lon', float, ('nlat'), zlib=True)
	fh.createVariable('SIGMA_SALT', float, ('depth', 'nlat'), zlib=True)
	fh.createVariable('VAR_SALT_'+str(branch1)+'', float, ('depth', 'nlat'), zlib=True)
	fh.createVariable('VAR_SALT_'+str(branch2)+'', float, ('depth', 'nlat'), zlib=True)
	fh.createVariable('STAT_LEVENE', float, ('depth', 'nlat'), zlib=True)
	fh.createVariable('p_value_LEVENE', float, ('depth', 'nlat'), zlib=True)

	fh.variables['depth'].long_name          		= 'Depth of midpoint'
	fh.variables['lon'].long_name            		= 'Array of longitudes'
	fh.variables['lat'].long_name            		= 'Array of latitudes'
	fh.variables['SIGMA_SALT'].long_name     		= 'Single estimator of variance salinity'
	fh.variables['VAR_SALT_'+str(branch1)+''].long_name   	= 'Variance of branch '+str(branch1)+' (model year 350-500)'
	fh.variables['VAR_SALT_'+str(branch2)+''].long_name  	= 'Variance of branch '+str(branch2)+' (model year 350-500)'
	fh.variables['STAT_LEVENE'].long_name    		= 'The test statistic of the Levenes test'
	fh.variables['p_value_LEVENE'].long_name 		= 'The p-value of the Levenes test'

	fh.variables['depth'].units             = 'm'
	fh.variables['lon'].units               = 'Degrees E'
	fh.variables['lat'].units               = 'Degrees N'
	fh.variables['SIGMA_SALT'].units        = '-'

	#Writing data to correct variable       
	fh.variables['nlat'][:]                 	= np.arange(len(lat))
	fh.variables['lon'][:]                  	= lon
	fh.variables['lat'][:]                  	= lat
	fh.variables['depth'][:]                	= depth
	fh.variables['SIGMA_SALT'][:]      		= ratio_var_salt
	fh.variables['VAR_SALT_'+str(branch1)+''][:]    = var_salt1
	fh.variables['VAR_SALT_'+str(branch2)+''][:]    = var_salt2
	fh.variables['STAT_LEVENE'][:]      		= stat
	fh.variables['p_value_LEVENE'][:]      		= p

	fh.close()	








	
