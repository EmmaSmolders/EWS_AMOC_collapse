#Script which determines single estimator of variance of salinity branch 1 and 2 using model years 350 to 500

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats

#Making pathway to folder with all data
directory_data  = '/home/smolders/CESM_Collapse/Data/CESM/Ocean/TEMP_SALT_ATLANTIC_transects/'
directory	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_variance/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	#Get the t-grid
	depth   	= fh.variables['depth'][:] 	 		#Depth (m)
	lon 		= fh.variables['lon'][:]			#Longitude
	lat 		= fh.variables['lat'][:]			#Latitude 
	time		= fh.variables['time'][:]			#Starting year 
	salt		= fh.variables['SALT'][:,:,:]		 	#Salinity 

	fh.close()

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

#CHOOSE TRANSECT (OSNAP OR SAMBA)
transect = 'OSNAP'
#transect = '34S'

#Read in yearly averaged data
if branch1 == 600:
	#SAMBA
	if transect == '34S':	
		lon, lat, depth, time1, salt1		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y600_600-1200.nc')
		time1 = time1[350:500]	
		salt1 = salt1[350:500,:,:]
	
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time1, salt1		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_950-1100_branch600.nc')
		time1 = time1[0:150]
		salt1 = salt1[0:150,:,:]
	
if branch1 == 1500:
	#SAMBA
	if transect == '34S':
		lon, lat, depth, time1, salt1		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y1500_1500-2100.nc')
		time1 = time1[350:500]	
		salt1 = salt1[350:500,:,:]
		
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time1, salt1		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_1850-2000_branch1500.nc')
		time1 = time1[0:150]
		salt1 = salt1[0:150,:,:]
	
if branch2 == 1650:
	#SAMBA
	if transect == '34S':
		lon, lat, depth, time2, salt2		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y1650_2900-3050.nc')
		
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time2, salt2		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_2900-3050_branch_1650.nc')

if len(time1) != len(time2):
	print('Time lengths are not the same!')
	sys.exit()

#Empty arrays
ratio_var_salt = ma.masked_all((len(depth), len(lon)))
var_salt1 = ma.masked_all((len(depth), len(lon)))
var_salt2 = ma.masked_all((len(depth), len(lon)))
stat = ma.masked_all((len(depth), len(lon)))
p = ma.masked_all((len(depth), len(lon)))
#sigma_mean_salt = ma.masked_all((len(depth), len(lat), len(lon[0])))
#ratio_sig_salt = ma.masked_all((len(depth), len(lat), len(lon[0])))

#Determine single estimator (VAR_E2/VAR_E1) for each grid point
for depth_i in range(len(depth)):
	print(depth_i)
	for lon_i in range(len(lon)):

		if salt1[0, depth_i, lon_i] is ma.masked:
				#If masked elements (=land), skip
				continue
		
		#Variance
		var_salt1[depth_i, lon_i] = np.var(TrendRemover(time1, salt1[:, depth_i, lon_i]))
		var_salt2[depth_i, lon_i] = np.var(TrendRemover(time2, salt2[:, depth_i, lon_i]))

		#Levene's test to test statistical difference between E1 and E2
		stat[depth_i, lon_i], p[depth_i, lon_i] = stats.levene(TrendRemover(time1, salt1[:, depth_i, lon_i]), TrendRemover(time2, salt2[:, depth_i, lon_i]))
		
		#Single estimator
		ratio_var_salt[depth_i, lon_i] = var_salt2[depth_i, lon_i]/var_salt1[depth_i, lon_i]

	#-----------------------------------------------------------------------------------------

	fh = netcdf.Dataset(directory+'Atlantic_single_estimator_variance_SALT_month_'+str(month_start)+'-'+str(month_end)+'_'+str(transect)+'_branch_'+str(branch2)+'_'+str(branch1)+'.nc', 'w')

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








	
