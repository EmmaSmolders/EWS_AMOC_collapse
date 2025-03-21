#Program determines the maximum of the AMOC strength (under 500m to exlcude surface currents)

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data (choose directory_data according to input data you want (quasi-equilibrium or equilibrium simulations)
directory_data	= '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.qe_hosing_branched_off_y600.001/run/'
#directory_data	= '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.qe_hosing.001/run/'
#directory_data = '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.start_500y_after_0.48Sv_run_now_0.495Sv.001/OUTPUT/ocn/hist/monthly/'
#directory_data	= '/projects/0/prace_imau/prace_2013081679/cesm1_0_5/b.e10.B1850.f19_g17.qe_hosing_branched_off_y1500.001/run/'
directory	= '/home/smolders/CESM_Collapse/Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	lat 		= fh.variables['lat_aux_grid'][84:]	#Latitude  
	depth   	= fh.variables['moc_z'][33:] / 100	#Depth (m)
	AMOC		= fh.variables['MOC'][0, 1, 0, 33:, 84:]

	fh.close()

	return lat, depth, AMOC

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

year_start	= 600
year_end	= 1200

#-----------------------------------------------------------------------------------------

files = glob.glob(directory_data+'*pop.h.*.nc')
files.sort()

for i in range(len(files)):
	if len(files[i]) != 174: 
		break

files	= files[:i]

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

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat 		  = fh.variables['lat_aux_grid'][84:]	#Latitude  
depth   	= fh.variables['moc_z'][33:] / 100	#Depth (m)
AMOC		  = fh.variables['MOC'][0, 1, 0, 33:, 84:]

fh.close()

#-----------------------------------------------------------------------------------------
time_year	  = ma.masked_all(int(len(time)/12))
AMOC_year	  = ma.masked_all((len(time_year), len(depth), len(lat)))
AMOC_max	  = ma.masked_all((len(time_year), len(lat)))
depth_max 	= ma.masked_all((len(time_year), len(lat)))

for year_i in range(int(np.min(time)), int(np.min(time))+len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i - int(np.min(time))] = year_i
	files_month = glob.glob(directory_data+'b.e10.B1850.f19_g17.qe_hosing_branched_off_y600.001.pop.h.'+str(year_i).zfill(4)+'-*.nc')
	files_month.sort()

	AMOC	=    ma.masked_all((12, len(depth), len(lat)))

	for file_i in range(len(files_month)):
		lat, depth, AMOC_month 	= ReadinData(files_month[file_i])
		AMOC[file_i]		= AMOC_month

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lat)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	AMOC_year[year_i - int(np.min(time))]	= np.sum(AMOC * month_days_all, axis = 0)

	#Now take the maximal value
	for lat_i in range(len(lat)):
		AMOC_max[year_i - int(np.min(time)), lat_i] = np.max(AMOC_year[year_i - int(np.min(time)), :, lat_i])
		max_index = np.argmax(AMOC_year[year_i - int(np.min(time)), :, lat_i])
		depth_max[year_i - int(np.min(time)), lat_i] = depth[max_index]
	        
#-----------------------------------------------------------------------------------------    

print('Data is written to file')
fh = netcdf.Dataset(directory+'/Ocean/AMOC_max_year_'+str(year_start)+'-'+str(year_end)+'_branch_600.nc', 'w')

fh.createDimension('time', len(time_year))
fh.createDimension('lat', len(lat))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('depth', float, ('time', 'lat'), zlib=True)
fh.createVariable('AMOC_max', float, ('time', 'lat'), zlib=True)

fh.variables['AMOC_max'].long_name 		= 'maximum of Atlantic Meridional Overturning Circulation'
fh.variables['depth'].long_name 	  	= 'Depth of maximum AMOC'
fh.variables['lat'].long_name 		  	= 'Array of latitudes'

fh.variables['time'].units 		      	= 'Year'
fh.variables['AMOC_max'].units 		 	= 'Sv'
fh.variables['depth'].units 		    	= 'm'
fh.variables['lat'].units 		     	= 'degrees N'

#Writing data to correct variable	
fh.variables['time'][:]     	  	  	= time_year
fh.variables['depth'][:] 		        = depth_max
fh.variables['lat'][:] 			        = lat
fh.variables['AMOC_max'][:] 		    	= AMOC_max

fh.close()
