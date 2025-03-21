#Script which determines single estimator of AC1 of temperature branch 1 and 2 using model years 350 to 500

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
import statsmodels.api as sm

#Making pathway to folder with all data
directory_data  = '/home/smolders/CESM_Collapse/Data/CESM/Ocean/TEMP_SALT_ATLANTIC_transects/'
directory	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_variance/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	#Get the t-grid
	depth   	= fh.variables['depth'][:] 	 		#Depth (m)
	lon 		= fh.variables['lon'][:]			#Longitude
	lat 		= fh.variables['lat'][:]			#Latitude 
	time		= fh.variables['time'][:]			#Time [model years]
	temp		= fh.variables['TEMP'][:,:,:]	 		#Temperature

	fh.close()

	return lon, lat, depth, time, temp

def TrendRemover(time, data):
	"""Removes linear trend"""

	time_fit	= np.arange(1, len(time) + 1)

	#Determine the detrended time series
	trend, base 	= polyfit(time_fit, data, 1)
	data_res	= data - ((trend * time_fit) + base)

	return data_res

def Fourrier_surrogates(data, num_surr):
	"""Takes the fourier transorm of the time series and re-shuffles the statistics"""
	data_fourier        = np.fft.rfft(data)
	random_phases       = np.exp(np.random.uniform(0, 2 * np.pi, (num_surr, len(data) // 2 + 1)) * 1.0j)
	data_fourier_new    = data_fourier * random_phases
	data_new            = np.real(np.fft.irfft(data_fourier_new))

	return data_new

#----------------------------------------------------------------------------------------------------------------

#Yearly averaged data from month 1 to 12
month_start	= 1
month_end	= 12
num_surr 	= 2000

#FILL HERE IN WHICH BRANCHES YOU WANT TO USE (E1 = 600, E2 = 1500, E3 = 1650)
branch1 = 1500
branch2 = 1650

#CHOOSE TRANSECT (OSNAP OR SAMBA)
#transect = 'OSNAP'
transect = '34S'

#Read in yearly averaged data
if branch1 == 600:
	#SAMBA
	if transect == '34S':	
		lon, lat, depth, time1, temp1		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y600_600-1200.nc')
		time1 = time1[350:500]	
		temp1 = temp1[350:500,:,:]
	
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time1, temp1		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_950-1100_branch600.nc')
		time1 = time1[0:150]
		temp1 = temp1[0:150,:,:]
	
if branch1 == 1500:
	#SAMBA
	if transect == '34S':
		lon, lat, depth, time1, temp1		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y1500_1500-2100.nc')
		time1 = time1[350:500]	
		temp1 = temp1[350:500,:,:]
		
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time1, temp1		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_1850-2000_branch1500.nc')
		time1 = time1[0:150]
		temp1 = temp1[0:150,:,:]
		
if branch2 == 1500:
	#SAMBA
	if transect == '34S':
		lon, lat, depth, time2, temp2		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y1500_1500-2100.nc')
		time2 = time2[350:500]	
		temp2 = temp2[350:500,:,:]
		
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time2, temp2		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_1850-2000_branch1500.nc')
		time2 = time2[0:150]
		temp2 = temp2[0:150,:,:]
	
if branch2 == 1650:
	#SAMBA
	if transect == '34S':
		lon, lat, depth, time2, temp2		= ReadinData(directory_data + 'TEMP_SALT_34S_year_branched_off_y1650_2900-3050.nc')
		
	#OSNAP
	if transect == 'OSNAP':
		lon, lat, depth, time2, temp2		= ReadinData(directory_data + 'TEMP_SALT_DENS_closest_to_OSNAP_year_2900-3050_branch_1650.nc')

if len(time1) != len(time2):
	print('Time lengths are not the same!')
	sys.exit()

#Empty arrays
ratio_ac_temp_surr = ma.masked_all((num_surr, len(depth), len(lon)))
ratio_ac_temp = ma.masked_all((len(depth), len(lon)))
ac_temp1 = ma.masked_all((num_surr, len(depth), len(lon)))
ac_temp2 = ma.masked_all((num_surr, len(depth), len(lon)))
data1 = ma.masked_all((num_surr, len(time1)))
data2 = ma.masked_all((num_surr, len(time2)))
sig_ratio_ac_temp_pos = ma.masked_all((len(depth), len(lon)))
sig_ratio_ac_temp_neg = ma.masked_all((len(depth), len(lon)))

#Determine single estimator (AC1_E2/AC1_E1) for each grid point
for depth_i in range(len(depth)):
	print(depth_i)
	for lon_i in range(len(lon)):

		if temp1[0, depth_i, lon_i] is ma.masked:
				#If masked elements (=land), skip
				continue

		autocorr_temp1 = sm.tsa.acf(TrendRemover(time1, temp1[:, depth_i, lon_i]))[1]
		autocorr_temp2 = sm.tsa.acf(TrendRemover(time2, temp2[:, depth_i, lon_i]))[1]

		#Only use autocorrelations higher than 0.5 to avoid high ratios due to small numbers
		#if autocorr_temp1 and autocorr_temp2 >= 0.5:
		#	ratio_ac_temp[depth_i, lon_i] = autocorr_temp2/autocorr_temp1
		#else:
		#	continue
		
		#No AC restriction
		ratio_ac_temp[depth_i, lon_i] = autocorr_temp2/autocorr_temp1
		
		#Generate Fourier surrogates of linearly detrended salinity time series to determine surrogate ratios
		data1[:, :] = Fourrier_surrogates(TrendRemover(time1, temp1[:, depth_i, lon_i]), num_surr) #shape (num_surr, time)
		data2[:, :] = Fourrier_surrogates(TrendRemover(time2, temp2[:, depth_i, lon_i]), num_surr)

		#Lag-1 autocorrelation
		for i in range(num_surr):
			ac_temp1[i, depth_i, lon_i] = sm.tsa.acf(data1[i, :])[1] #shape (num_surr, depth, lon)
			ac_temp2[i, depth_i, lon_i] = sm.tsa.acf(data2[i, :])[1]

			#Single estimator, surrogate ratio's
			ratio_ac_temp_surr[i, depth_i, lon_i] = ac_temp2[i, depth_i, lon_i]/ac_temp1[i, depth_i, lon_i]

		#Reject null-hypothesis if probability R>1 is larger than 0.95 (i.e. when less than 5% of the AC distributions overlap)
		sig_ratio_ac_temp1 = ma.masked_where(ratio_ac_temp_surr[:,depth_i,lon_i] <= 1, ratio_ac_temp_surr[:,depth_i,lon_i])
		#print(sig_ratio_ac_temp1)
		sig_ratio_ac_temp1 = np.sum(sig_ratio_ac_temp1 * 0.0 + 1.0, axis = 0) #/ (num_surr * num_surr)
		#print(sig_ratio_ac_temp1)
		sig_ratio_ac_temp_pos[depth_i, lon_i] = ma.masked_where(sig_ratio_ac_temp1 < 0.95*num_surr, sig_ratio_ac_temp1)
		#print(sig_ratio_ac_temp[depth_i, lon_i])
        	
		#Reject null-hypothesis if probability R<1 is larger than 0.95
		sig_ratio_ac_temp1 = ma.masked_where(ratio_ac_temp_surr[:,depth_i,lon_i] >= 1, ratio_ac_temp_surr[:,depth_i,lon_i])
		sig_ratio_ac_temp1 = np.sum(sig_ratio_ac_temp1 * 0.0 + 1.0, axis = 0)# / (num_surr * num_surr)
		sig_ratio_ac_temp_neg[depth_i, lon_i] = ma.masked_where(sig_ratio_ac_temp1 < 0.95*num_surr, sig_ratio_ac_temp1)
	

	#-----------------------------------------------------------------------------------------

	fh = netcdf.Dataset(directory+'Atlantic_single_estimator_AC1_TEMP_month_'+str(month_start)+'-'+str(month_end)+'_'+str(transect)+'_branch_'+str(branch2)+'_'+str(branch1)+'.nc', 'w')

	fh.createDimension('depth', len(depth))
	fh.createDimension('nlat', len(lat))
	fh.createDimension('nlon', len(lon))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('nlat', float, ('nlat'), zlib=True)
	fh.createVariable('nlon', float, ('nlon'), zlib=True)
	fh.createVariable('lat', float, ('nlat'), zlib=True)
	fh.createVariable('lon', float, ('nlat'), zlib=True)
	fh.createVariable('RATIO_AC1_TEMP', float, ('depth', 'nlon'), zlib=True)
	fh.createVariable('SIG_RATIO_AC1_TEMP_pos', float, ('depth', 'nlon'), zlib=True)
	fh.createVariable('SIG_RATIO_AC1_TEMP_neg', float, ('depth', 'nlon'), zlib=True)

	fh.variables['depth'].long_name          		= 'Depth of midpoint'
	fh.variables['lon'].long_name            		= 'Array of longitudes'
	fh.variables['lat'].long_name            		= 'Array of latitudes'
	fh.variables['RATIO_AC1_TEMP'].long_name     		= 'Single estimator of AC1 salinity'
	fh.variables['SIG_RATIO_AC1_TEMP_pos'].long_name   	= '95 CI significance of ratio AC1 salinity R>1'
	fh.variables['SIG_RATIO_AC1_TEMP_neg'].long_name   	= '95 CI significance of ratio AC1 salinity R<1'

	fh.variables['depth'].units             = 'm'
	fh.variables['lon'].units               = 'Degrees E'
	fh.variables['lat'].units               = 'Degrees N'
	fh.variables['RATIO_AC1_TEMP'].units    = '-'

	#Writing data to correct variable       
	fh.variables['nlat'][:]                 	= np.arange(len(lat))
	fh.variables['lon'][:]                  	= lon
	fh.variables['lat'][:]                  	= lat
	fh.variables['depth'][:]                	= depth
	fh.variables['RATIO_AC1_TEMP'][:]      		= ratio_ac_temp
	fh.variables['SIG_RATIO_AC1_TEMP_pos'][:]   	= sig_ratio_ac_temp_pos
	fh.variables['SIG_RATIO_AC1_TEMP_neg'][:]   	= sig_ratio_ac_temp_neg

	fh.close()	








	
