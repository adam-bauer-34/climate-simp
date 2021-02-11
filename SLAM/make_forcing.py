import numpy as np
import math as m
import matplotlib.pyplot as plt
import random
import scipy.stats as stat
#from the_model import SLAM
import scipy.signal as signal
#from T_constant_SLAM import T_constant_SLAM
from netCDF4 import Dataset

def make_precip():

	mean_P = 2000 # in mm, cumulative precip in one summer
	std_P  = 500  # in mm, standard deviation of cumulative summer precip

	ind_mean = 10 # in mm, mean value of precipitation event

	Cumulative_rain = np.random.normal(mean_P,std_P)
	Precip_events = []
	i = 0
	while Cumulative_rain > sum(Precip_events):
		Precip_events.append(np.random.gamma(ind_mean,scale=1.0))
		i = i + 1

	return(Precip_events)

def red_forcing():


######### Parameters  ##############################################################################################

	sigma_limit = 900
	sigma_avg   = 800	
	sigma_min   = 200

	rad_anom_norm = 200
	day_length = 15  # IN HOURS

	temp_mean = 292
	temp_anom_norm = 4

	hum_mean = 0.012
	hum_anom_norm = 0.002

	r_temp = 0.6
	r_hum  = 0.6
	r_rad  = 0.6


########################### Setting up some basic cycles  #############################################################

	space_time = (24-day_length)*60/2.

	day_time_shift = np.linspace(space_time,day_length*60 + space_time - 1,day_length*60)
	day_time = np.linspace(1,day_length*60,day_length*60)
	day_MAX  = sigma_limit*np.sin(m.pi*day_time/(day_length*60))
	day_avg  = sigma_avg*np.sin(m.pi*day_time/(day_length*60))
	day_min  = sigma_min*np.sin(m.pi*day_time/(day_length*60))

	pre_day = np.linspace(0,space_time - 1,space_time)
	post_day = np.linspace(day_length*60 + space_time,24*60 - 1, space_time)

	no_rad = np.zeros(len(pre_day))
	
	whole_day = np.hstack((pre_day,day_time_shift,post_day))
	whole_day_max = np.hstack((no_rad,day_MAX,no_rad))
	whole_day_avg = np.hstack((no_rad,day_avg,no_rad))
	whole_day_min = np.hstack((no_rad,day_min,no_rad))

	all_time = 92*24*60

	whole_summer_time = np.linspace(0,92*24*60 - 1,all_time)
	sw_max = np.tile(whole_day_max,92)
	sw_avg = np.tile(whole_day_avg,92)
	sw_min = np.tile(whole_day_min,92)
	
################################## ADDING RED NOISE FORCING ########################################################

	six_hour_len = 92*4

	temp_day_vals = np.linspace(0,3*m.pi/2.,4)
	temp_day = np.sin(temp_day_vals)
	temp_full = np.tile(temp_day,92)
	hum_full  = np.tile(temp_day,92)
	rad_full = np.tile(temp_day,92)

	i = 1
	while i < six_hour_len:

		t_random = random.gauss(0,1)
		q_random = random.gauss(0,1) # FOR UNCORRELATED T AND Q
		rad_random = random.gauss(0,1)

		temp_full[i] = r_temp*temp_full[i-1] + ((1-r_temp**2)**0.5)*t_random
		rad_full[i] = r_rad*rad_full[i-1] + ((1-r_rad**2)**0.5)*rad_random
		hum_full[i] = r_hum*hum_full[i-1] + ((1-r_hum**2)**0.5)*q_random # FOR UNCORRELATED T AND Q

		i = i + 1


	temp_forcing = interp(temp_mean + temp_full*temp_anom_norm,92*60*24)
	q_forcing = interp(hum_mean + hum_full*hum_anom_norm,92*60*24)
	rad_full = rad_full*rad_anom_norm

	rad_full = interp(rad_full,92*60*24)
	precip = np.zeros(92*24*60)
	cf = np.zeros(92*24*60)
	F_solar = np.zeros(92*24*60)

	cloud_count = 0
	cloud_atlas = []  # tells us where clouds are so we can re-correlate P and F

	i = 0
	while i < 92*24*60:

		do_add = sw_avg[i] + rad_full[i]
		max_val = sw_max[i]
		avg_val = sw_avg[i]

		if rad_full[i] < -rad_anom_norm*1.25:  # CLOUD RED NOISE HAS TO BE LARGE FOR CLOUD ATLAS TO GENERATE PRECIP
			cloud_atlas.append(i)
			cloud_count = cloud_count + 1

		F_solar[i] = do_add

		if F_solar[i] < 0:
			F_solar[i] = 0

		if F_solar[i] > sw_max[i]:
			F_solar[i] = sw_max[i]

		if F_solar[i] < sw_min[i]:
			F_solar[i] = sw_min[i]
			
		cf[i] = -rad_full[i]/(sigma_limit*2)
		
		if cf[i] < 0:
			cf[i] = 0
		if cf[i] > 1:
			cf[i] = 1


		i = i + 1

	precip_array = make_precip()

	j = 0
	print(cloud_count)
	while j < len(precip_array):
		my_cloud_dex = random.randint(0,cloud_count)
		precip_dex = cloud_atlas[my_cloud_dex]
		precip[precip_dex] = precip_array[j]
		j = j + 1

	return(F_solar,precip,temp_forcing,q_forcing,cf)


def interp(array,length_new):


	time_old = np.linspace(0,length_new - 1,len(array))
	time_new = np.linspace(0,length_new - 1,length_new)
	new_array = np.interp(time_new,time_old,array)
	return(new_array)


def write_netCDF(nyears):

# Creates an ensemble of forcings for a given number of years

	L = 24*60*92 # minutes in a summer

	rad = np.zeros(shape=(L,nyears))
	precip = np.zeros(shape=(L,nyears))
	t_forcing = np.zeros(shape=(L,nyears))
	q_forcing = np.zeros(shape=(L,nyears))
	cf = np.zeros(shape=(L,nyears))

	i = 0
	while i < nyears:
		one_rad,one_precip,one_t,one_q,one_cf = red_forcing()
		rad[:,i] = one_rad
		precip[:,i] = one_precip
		t_forcing[:,i] = one_t
		q_forcing[:,i] = one_q
		cf[:,i] = one_cf
		i = i + 1

	your_netcdf_name = 'Fresh_Adam'

	my_temp = Dataset(your_netcdf_name+'.nc', 'w', format='NETCDF3_64BIT')

	my_temp.createDimension('time', 24*60*92)
	minutes = my_temp.createVariable('time', 'i4', ('time', ))
	my_temp.createDimension('summer_num',nyears)
	number = my_temp.createVariable('summer_num','i4',('summer_num', ))

	rad_var = my_temp.createVariable('F','f4',('time','summer_num', ))
	tup_var = my_temp.createVariable('Tup','f4',('time','summer_num', ))
	qup_var = my_temp.createVariable('qup','f4',('time','summer_num', ))
	precip_var = my_temp.createVariable('precip','f4',('time','summer_num', ))
	cf_var   = my_temp.createVariable('cf','f4',('time','summer_num', ))

	rad_var[:,:] = rad
	precip_var[:,:] = precip
	cf_var[:,:]   = cf
	tup_var[:,:] = t_forcing
	qup_var[:,:] = q_forcing

	my_temp.close()


