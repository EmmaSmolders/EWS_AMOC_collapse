#Script which determines single estimator of lag-1 autocorrelation of temperature branch 1 and 2 using model years 350 to 500 (for RAPID only)

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
directory_data  	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_surface/'
directory		= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/ATLANTIC_variance/'
directory_data_1650  	= '/home/smolders/CESM_Collapse/Data/CESM/Ocean/TEMP_SALT_ATLANTIC_transects/'

def find_nearest(array, value):
    #array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def AtlanticRegion(lon, lat, data):
	
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
	temp		= fh.variables['TEMP'][:,:,:]		 	#Salinity 

	fh.close()
	
	if transect == False:
		lon, lat, temp = AtlanticRegion(lon, lat, temp)

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

##----------------------------------------------------------------------------------------------------------------

#Yearly averaged data from month 1 to 12
month_start	= 1
month_end	= 12
num_surr 	= 2000

#FILL HERE IN WHICH BRANCHES YOU WANT TO USE (E1 = 600, E2 = 1500, E3 = 1650)
branch1 = 600
branch2 = 1650

#Which index?
rapid_idx = 270 #(26N)

#Read in data of branch 1 and 2 (yearly averged, only last 150 years). I do it in this way for 26N, because the data is not on Snellius anymore and I haven't saved the data per transect (only for 34S)
if branch1 == 600:
	lon, lat, depth, time1, temp1		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y600_950-1000.nc', False)
	lon, lat, depth, time2, temp2		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y600_1001-1050.nc', False)
	lon, lat, depth, time3, temp3		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y600_1051-1100.nc', False)
	
	# Concatenate arrays and delete
	time_branch1 = np.concatenate((time1, time2, time3[0:-1]), axis=0) #Take only 150 years, instead of 151 (even number for Fourier_surrogates function)
	temp_branch1 = np.concatenate((temp1, temp2, temp3[0:-1]), axis=0)
	
	#Select transect
	lat1 = lat[rapid_idx, :]
	lon1 = lon[rapid_idx, :]
	print(lat1)
	print(lon1)
	temp_branch1 = temp_branch1[:,:,rapid_idx,:]

	del(temp1, temp2, temp3)

if branch1 == 1500:
	lon, lat, depth, time4, temp4		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1850-1899.nc', False)
	lon, lat, depth, time5, temp5		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1900-1950.nc', False)
	lon, lat, depth, time6, temp6		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1951-2000.nc', False)

	time_branch1 = np.concatenate((time4, time5, time6[0:-1]), axis=0)
	temp_branch1 = np.concatenate((temp4, temp5, temp6[0:-1]), axis=0)
	
	#Select transect
	lat1 = lat[rapid_idx, :]
	lon1 = lon[rapid_idx, :]
	print(lat1)
	print(lon1)
	temp_branch1 = temp_branch1[:,:,rapid_idx,:]

	del(temp4, temp5, temp6)
	
if branch2 == 1500:
	lon, lat, depth, time7, temp7		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1850-1899.nc', False)
	lon, lat, depth, time8, temp8		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1900-1950.nc', False)
	lon, lat, depth, time9, temp9		= ReadinData(directory_data + 'TEMP_SALT_DENS_surface_year_branched_off_y1500_1951-2000.nc', False)

	time_branch2 = np.concatenate((time7, time8, time9[0:-1]), axis=0)
	temp_branch2 = np.concatenate((temp7, temp8, temp9[0:-1]), axis=0)
	
	#Select transect
	lat = lat[rapid_idx, :]
	lon = lon[rapid_idx, :]
	print(lat)
	print(lon)
	temp_branch2 = temp_branch2[:,:,rapid_idx,:]

	del(temp7, temp8, temp9)
	
if branch2 == 1650:
	lon, lat, depth, time_branch2, temp_branch2		= ReadinData(directory_data_1650 + 'TEMP_SALT_26N_year_branched_off_y1650_2900-3050.nc', True)
	print(lat)
	print(lon)

if not np.allclose(lat, lat1, atol=1e-10): #apparently there can be a differnce up to e-14?
    print('Latitude array is not the same')
    #sys.exit()
	
if not np.array_equal(lon, lon1):
    print('Latitude array is not the same')
    #sys.exit()
	
if np.shape(temp_branch2) != np.shape(temp_branch1):
	print('Length of array is not the same, and should be equal to 150')
	#sys.exit()
	
print(lat)
print(lon)
print(np.shape(temp_branch2))
print(np.shape(temp_branch1))

#Empty arrays
ratio_ac_temp_surr = ma.masked_all((num_surr, len(depth), len(lon)))
ratio_ac_temp = ma.masked_all((len(depth), len(lon)))
ac_temp1 = ma.masked_all((num_surr, len(depth), len(lon)))
ac_temp2 = ma.masked_all((num_surr, len(depth), len(lon)))
data1 = ma.masked_all((num_surr, len(time_branch1), len(depth), len(lon)))
data2 = ma.masked_all((num_surr, len(time_branch2), len(depth), len(lon)))
sig_ratio_ac_temp_pos = ma.masked_all((len(depth), len(lon)))
sig_ratio_ac_temp_neg = ma.masked_all((len(depth), len(lon)))

#Determine single estimator (AC1_E2/AC1_E1) for each grid point
for depth_i in range(len(depth)):
	print(depth_i)
	for lon_i in range(len(lon)):

		if temp_branch2[0, depth_i, lon_i] is ma.masked:
				#If masked elements (=land), skip
				continue

		autocorr_temp1 = sm.tsa.acf(TrendRemover(time_branch1, temp_branch1[:, depth_i, lon_i]))[1]
		autocorr_temp2 = sm.tsa.acf(TrendRemover(time_branch2, temp_branch2[:, depth_i, lon_i]))[1]

		#Only use autocorrelations higher than 0.5 to avoid high ratios due to small numbers
		#if autocorr_temp1 and autocorr_temp2 >= 0.5:
		#	ratio_ac_temp[depth_i, lon_i] = autocorr_temp2/autocorr_temp1
		#else:
		#	continue
		
		#No AC restriction
		ratio_ac_temp[depth_i, lon_i] = autocorr_temp2/autocorr_temp1
		
		#Generate Fourier surrogates of linearly detrended salinity time series to determine surrogate ratios
		data1[:, :, depth_i, lon_i] = Fourrier_surrogates(TrendRemover(time_branch1, temp_branch1[:, depth_i, lon_i]), num_surr) #shape (num_surr, time, depth, lon)
		data2[:, :, depth_i, lon_i] = Fourrier_surrogates(TrendRemover(time_branch2, temp_branch2[:, depth_i, lon_i]), num_surr)

		#Lag-1 autocorrelation
		for i in range(num_surr):
			ac_temp1[i, depth_i, lon_i] = sm.tsa.acf(data1[i, :, depth_i, lon_i])[1] #shape (num_surr, depth, lon)
			ac_temp2[i, depth_i, lon_i] = sm.tsa.acf(data2[i, :, depth_i, lon_i])[1]

			#Single estimator, surrogate ratio's
			ratio_ac_temp_surr[i, depth_i, lon_i] = ac_temp2[i, depth_i, lon_i]/ac_temp1[i, depth_i, lon_i]

		#Reject null-hypothesis if probability R>1 is larger than 0.95
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

	fh = netcdf.Dataset(directory+'Atlantic_single_estimator_AC1_TEMP_month_'+str(month_start)+'-'+str(month_end)+'_26N_branch_'+str(branch2)+'_'+str(branch1)+'.nc', 'w')

	fh.createDimension('depth', len(depth))
	fh.createDimension('nlat', len(lat))
	fh.createDimension('nlon', len(lon))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	#fh.createVariable('nlat', float, ('nlat'), zlib=True)
	fh.createVariable('nlon', float, ('nlon'), zlib=True)
	fh.createVariable('lat', float, ('nlon'), zlib=True)
	fh.createVariable('lon', float, ('nlon'), zlib=True)
	fh.createVariable('RATIO_AC1_TEMP', float, ('depth', 'nlon'), zlib=True)
	fh.createVariable('SIG_RATIO_AC1_TEMP_pos', float, ('depth', 'nlon'), zlib=True)
	fh.createVariable('SIG_RATIO_AC1_TEMP_neg', float, ('depth', 'nlon'), zlib=True)

	fh.variables['depth'].long_name          		= 'Depth of midpoint'
	fh.variables['lon'].long_name            		= 'Array of longitudes'
	fh.variables['lat'].long_name            		= 'Array of latitudes'
	fh.variables['RATIO_AC1_TEMP'].long_name     		= 'Single estimator of AC1 temperature'
	fh.variables['SIG_RATIO_AC1_TEMP_pos'].long_name   	= '95 CI significance of ratio AC1 temperature R>1'
	fh.variables['SIG_RATIO_AC1_TEMP_neg'].long_name   	= '95 CI significance of ratio AC1 temperature R<1'

	fh.variables['depth'].units             = 'm'
	fh.variables['lon'].units               = 'Degrees E'
	fh.variables['lat'].units               = 'Degrees N'
	fh.variables['RATIO_AC1_TEMP'].units    = '-'

	#Writing data to correct variable       
	fh.variables['nlon'][:]                 	= np.arange(len(lon))
	fh.variables['lon'][:]                  	= lon
	fh.variables['lat'][:]                  	= lat
	fh.variables['depth'][:]                	= depth
	fh.variables['RATIO_AC1_TEMP'][:]      		= ratio_ac_temp
	fh.variables['SIG_RATIO_AC1_TEMP_pos'][:]   	= sig_ratio_ac_temp_pos
	fh.variables['SIG_RATIO_AC1_TEMP_neg'][:]   	= sig_ratio_ac_temp_neg

	fh.close()	







