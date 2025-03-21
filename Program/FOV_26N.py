#Program determines the FOV index at 26N

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
#directory_data	= '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.qe_hosing.001/run/'
#directory_data	= '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.qe_hosing_branched_off_y600.001/run/'
#directory_data	= '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.qe_hosing_branched_off_y1500.001/run/'
directory_data = '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.start_500y_after_0.48Sv_run_now_0.495Sv.001/OUTPUT/ocn/hist/monthly/'
directory	= '/home/smolders/CESM_Collapse/Data/CESM/'

def ReadinData(filename, depth_min_index, depth_max_index):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon_u 		= fh.variables['ULONG'][:]						#Longitude
	lat_u 		= fh.variables['ULAT'][:]						#Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 	 / 100.0	#Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 	 / 100.0	#Layer thickness (m)
	dx_u		= fh.variables['DXU'][:]  / 100.0		                        #Zonal grid cell length (m)
	depth_u 	= fh.variables['HU'][:] / 100.0						#Depth at u-grid (m)
	v_vel 		= fh.variables['VVEL'][0, depth_min_index:depth_max_index] / 100.0	#Meridional velocity (m/s)

	#Get the t-grid
	lon_t 		= fh.variables['TLONG'][:]						#Longitude
	lat_t 		= fh.variables['TLAT'][:]						#Latitude 
	dx_t		= fh.variables['DXT'][:]  / 100.0		                  	#Zonal grid cell length (m)
	depth_t 	= fh.variables['HT'][:] / 100.0						#Depth at t-grid (m)
	salt		= fh.variables['SALT'][0, depth_min_index:depth_max_index] 		#Salinity (g / kg)

	fh.close()	

	#Convert the grid to -180 to 180
	lon, lat_u	= ConverterField2D(lon_u, lat_u)
	lon, lat_t	= ConverterField2D(lon_t, lat_t)
	lon, depth_u	= ConverterField2D(lon_u, depth_u)
	lon, depth_t	= ConverterField2D(lon_t, depth_t)
	lon, dx_u	= ConverterField2D(lon_u, dx_u)
	lon, dx_t	= ConverterField2D(lon_t, dx_t)
	lon_u, v_vel	= ConverterField3D(lon_u, v_vel)
	lon_t, salt	= ConverterField3D(lon_t, salt)

	#Section at 26N
	lon_1, lon_2	= 88, 150
	lat_1, lat_2	= 269, 272

	#Now only select a small region near section
	lon_u	= lon_u[lat_1:lat_2, lon_1:lon_2]
	lat_u	= lat_u[lat_1:lat_2, lon_1:lon_2]
	lon_t	= lon_t[lat_1:lat_2, lon_1:lon_2]
	lat_t	= lat_t[lat_1:lat_2, lon_1:lon_2]
	depth_u	= depth_u[lat_1:lat_2, lon_1:lon_2]
	depth_t	= depth_t[lat_1:lat_2, lon_1:lon_2]
	dx_u	= dx_u[lat_1:lat_2, lon_1:lon_2]
	dx_t	= dx_t[lat_1:lat_2, lon_1:lon_2]
	v_vel	= v_vel[:, lat_1:lat_2, lon_1:lon_2]
	salt	= salt[:, lat_1:lat_2, lon_1:lon_2]

	for depth_i in range(len(depth)):
		v_vel[depth_i]	= ma.masked_where(depth_u <= 0.0, v_vel[depth_i])
		salt[depth_i]	= ma.masked_where(depth_t <= 0.0, salt[depth_i])

	#Now only get the relevant index
	lat_u	= lat_u[1]
	lat_t	= lat_t[1:3]
	depth_u	= depth_u[1]
	depth_t	= depth_t[1:3]
	dx_u	= dx_u[1]
	dx_t	= dx_t[1:3]
	v_vel	= v_vel[:, 1]
	salt	= salt[:, 1:3]

	return lon_u, lat_u, lon_t, lat_t, depth, layer, depth_u, depth_t, dx_u, dx_t, v_vel, salt

def ConverterField2D(lon, field):
	"""Shifts field, to -180E to 180E"""
	lon_new		= ma.masked_all(shape(lon))
	field_new	= ma.masked_all(shape(field))

	#Get the corresponding index
	index		= 195

	#Start filling at -180
	lon_new[:, :len(lon[0, index:])] 	= lon[:, index:]
	field_new[:, :len(lon[0, index:])]	= field[:, index:]

	#Fill the remaining part
	lon_new[:, len(lon[0, index:]):] 	= lon[:, :index]
	field_new[:, len(lon[0, index:]):]	= field[:, :index]

	lon_new[lon_new >= 179.9]	= lon_new[lon_new >= 179.9] - 360.0

	return lon_new, field_new

def ConverterField3D(lon, field):
	"""Shifts field, to -180E to 180E"""
	lon_new		= ma.masked_all(shape(lon))
	field_new	= ma.masked_all(shape(field))

	#Get the corresponding index
	index		= 195

	#Start filling at -180E
	lon_new[:, :len(lon[0, index:])] 	= lon[:, index:]
	field_new[:, :, :len(lon[0, index:])]	= field[:, :, index:]

	#Fill the remaining part to 360E
	lon_new[:, len(lon[0, index:]):] 	= lon[:, :index]
	field_new[:, :, len(lon[0, index:]):]	= field[:, :, :index]

	lon_new[lon_new >= 179.9]	= lon_new[lon_new >= 179.9] - 360.0

	return lon_new, field_new
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

year_start	= 2900
year_end	= 3051

#-----------------------------------------------------------------------------------------

files = glob.glob(directory_data+'*pop.h.*.nc')
files.sort()

for i in range(len(files)):
	if len(files[i]) != 218: 
		break

files	= files[:i]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-10:-3]	
	year  = int(date[0:4])
	month = int(date[5:7])

	time[year_i] = year + (month-1) / 12.0

time_start	= (np.abs(time - year_start)).argmin()
time_end	= (np.abs(time - (year_end+1))).argmin()

time		= time[time_start:time_end]
files		= files[time_start:time_end]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-10:-3]	
	year  = int(date[0:4])
	month = int(date[5:7])

	time[year_i] = year + (month-1) / 12.0
#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

depth   	= fh.variables['z_t'][:] / 100	#Depth (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (np.abs(depth_min - depth)).argmin()
depth_max_index 	= (np.abs(depth_max - depth)).argmin() + 1

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon_u, lat_u, lon_t, lat_t, depth, layer, depth_u, depth_t, grid_x_u, grid_x_t, v_vel, salt 	= ReadinData(files[0], depth_min_index, depth_max_index)
layer_field_u					= ma.masked_all((len(depth), len(lon_u[0])))
layer_field_t					= ma.masked_all((len(depth), 2, len(lon_t[0])))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field_u[depth_i]	= layer[depth_i]
	layer_field_u[depth_i]	= ma.masked_array(layer_field_u[depth_i], mask = v_vel[depth_i].mask)
	layer_field_t[depth_i]	= layer[depth_i]
	layer_field_t[depth_i]	= ma.masked_array(layer_field_t[depth_i], mask = salt[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff_u	= np.sum(layer_field_u, axis = 0) - depth_u
	depth_diff_t	= np.sum(layer_field_t, axis = 0) - depth_t

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff_u	= ma.masked_where(depth_diff_u < 0, depth_diff_u)
	depth_diff_u	= depth_diff_u.filled(fill_value = 0.0)
	depth_diff_t	= ma.masked_where(depth_diff_t < 0, depth_diff_t)
	depth_diff_t	= depth_diff_t.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field_u[depth_i]	= layer_field_u[depth_i] - depth_diff_u
	layer_field_t[depth_i]	= layer_field_t[depth_i] - depth_diff_t

#Normalise layer field per layer
layer_field_u_norm  = ma.masked_all(shape(layer_field_u))
grid_x_t_norm	    = ma.masked_all((len(depth), 2, len(lon_t[0])))
lat_weight	    = ma.masked_all((len(depth), 2))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_norm[depth_i]		= layer_field_u[depth_i] / np.sum(layer_field_u[depth_i])
    
	#Normalise the length
	grid_x_t_depth          		= ma.masked_array(grid_x_t, mask = salt[depth_i].mask)

	#Now get the lat weights
	lat_weight[depth_i]			= np.sum(grid_x_t_depth, axis = 1) / np.sum(grid_x_t_depth)

	for lat_i in range(2):
		grid_x_t_norm[depth_i, lat_i]  = grid_x_t_depth[lat_i] / np.sum(grid_x_t_depth[lat_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
time_year		= ma.masked_all(int(len(time)/12))
transport_salt_all	= ma.masked_all(len(time_year))

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] = int(time[year_i*12])
	files_month 	  = files[year_i*12:(year_i+1)*12]

	v_vel	=    ma.masked_all((12, len(depth), 2, len(lon_u[0])))
	salt	=    ma.masked_all((12, len(depth), 2, len(lon_t[0])))

	for file_i in range(len(files_month)):
		lon_u, lat_u, lon_t, lat_t, depth, layer, depth_u, depth_t, grid_x_u, grid_x_t, v_vel_month, salt_month = ReadinData(files_month[file_i], depth_min_index, depth_max_index)

		v_vel[file_i, :, 0]	= v_vel_month
		salt[file_i]		= salt_month

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= ma.masked_array(month_days, salt[:, 0, 0, 50].mask)
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(salt[0, 0]), len(salt[0, 0, 0])))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	v_vel	= np.sum(v_vel * month_days_all, axis = 0)
	salt	= np.sum(salt * month_days_all, axis = 0)
	v_vel	= v_vel[:, 0]

	#------------------------------------------------------------------------------

	#Determine the meridional transport
	transport	= v_vel * layer_field_u * grid_x_u

	#Determine the section averaged velocity (barotropic)
	vel_barotropic	= np.sum(transport) / np.sum(layer_field_u * grid_x_u)

	#Determine the overturning velocity (baroclinic)
	vel_baroclinic	 = v_vel - vel_barotropic

	#Determine the zonal means
	salt_zonal      = np.sum(salt * grid_x_t_norm, axis = 2)  - 35.0
	transport_clin	= np.sum(vel_baroclinic * layer_field_u * grid_x_u, axis = 1)

	#-----------------------------------------------------------------------------------------

	#Take the mean over the two latitudes for the salt
	salt_zonal      = np.sum(salt_zonal * lat_weight, axis = 1)

	#Determine the total salinity transport
	transport_salt_all[year_i]		= (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 

#-----------------------------------------------------------------------------------------

#plot(time_year, transport_salt_all, '-k')
#show()

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_26N_year_'+str(year_start)+'_'+str(year_end)+'_branch_1650.nc', 'w')

fh.createDimension('time', len(time_year))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('F_OV', float, ('time'), zlib=True)

fh.variables['F_OV'].long_name 		= 'Freshwater transport (overturning component)'

fh.variables['time'].units 		= 'Year'
fh.variables['F_OV'].units 		= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time_year
fh.variables['F_OV'][:] 		= transport_salt_all

fh.close()
