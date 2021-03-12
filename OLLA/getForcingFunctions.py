#!/usr/bin/env python
# coding: utf-8

# # Forcing Generator Functions
# 
# This file contains programs that make the .nc forcing file that will be fed to OLLA. 
# 
# The first function makes the .nc, but requires the latter function to make the forcing on a minute-by-minute scale. In order to make this 
# minute by minute forcing, we need minute by minute precipitation, 
# which is why we have the last function.
# 
# This code largely borrows from make_forcing.py, originally written by Lucas Zeppetello.
# 
# Included functions: getForcingNetCDF, getRedNoise, getPrecip, getInterpolate


import numpy as np 
import random

from netCDF4 import Dataset

random.seed() # Init seed for random number generator. If blank, the seed is the current system time.


# Generate a NetCDF to force the model with


def getForcingNetCDF(tot_years, stats_array, filename):
    
    # This creates a netCDF file with two dimensions: the time (or
    # minutes) and the year (taken from tot_years). This will force
    # OLLA.
    
    mins_in_sum = 92*60*24 # minutes in a summer
    
    # Initialize arrays for forcing variables
    
    F_solar = np.zeros(shape=(mins_in_sum, tot_years))
    precip = np.zeros(shape=(mins_in_sum, tot_years))
    
    #t_upper = np.zeros(shape=(mins_in_sum, tot_years))
    #q_upper = np.zeros(shape=(mins_in_sum, tot_years))
    #cloud_frac = np.zeros(shape=(mins_in_sum, tot_years))
    
    # Loop through initialized arrays and make "by the minute" forcing using 
    # getRedNoise(), below. So to be clear: this loop makes a year's worth of forcing,
    # while getRedNoise() makes forcing by the minute.
    
    for i in range(0, tot_years-1):
        minute_F_solar, minute_precip = getRedNoise(stats_array)
        F_solar[:,i] = minute_F_solar
        precip[:,i] = minute_precip
        
        #t_upper[:,i] = minute_t_upper
        #q_upper[:,i] = minute_q_upper
        #cloud_frac[:,i] = minute_cloud_frac
        
    # Now to write the above forcing information to a .nc
    
    tmp_dataset = Dataset(filename, 'w', format='NETCDF3_64BIT')
    
    # Creating a "dimension" is equivalent to inializing a numpy array of shape
    # 'time'. You just have to add the dimensions to the .nc sequentially, as done below.
    
    tmp_dataset.createDimension('time', mins_in_sum)
    minutes = tmp_dataset.createVariable('time', 'i4', ('time', ))
    
    tmp_dataset.createDimension('summer_number', tot_years)
    summer_number = tmp_dataset.createVariable('summer_number', 'i4', ('summer_number', ))
    
    # Initialize vairables in the .nc that relate to forcing.
    
    F_solar_var = tmp_dataset.createVariable('F_solar', 'f4', ('time', 'summer_number', ))
    precip_var = tmp_dataset.createVariable('precip_syn', 'f4', ('time', 'summer_number', ))
    
    #t_upper_var = tmp_dataset.createVariable('t_upper', 'f4', ('time', 'summer_number', ))
    #q_upper_var = tmp_dataset.createVariable('q_upper', 'f4', ('time', 'summer_number', ))
    #cloud_frac_var = tmp_dataset.createVariable('cloud_frac', 'f4', ('time', 'summer_number', ))
    
    # Fill in values of the forcing variables in the .nc
    
    F_solar_var[:,:] = F_solar 
    precip_var[:,:] = precip
    
    #t_upper_var[:,:] = t_upper
    #q_upper_var[:,:] = q_upper
    #cloud_frac_var[:,:] = cloud_frac
    
    # Voila!
    
    tmp_dataset.close()


# Generate red noise arrays for the forcing NetCDF 


def getRedNoise(stats_array):
    
    ### PARAMETER SETUP ###
    
    # Some parameters of getRedNoise are fed in from the SGP data, and thus are stored in 
    # stats_array.
    
    rad_anom_norm = stats_array[0]
    
    # they create an envelope of allowable sw radiation for the model.
    
    sw_limit = stats_array[1] # highest allowable daily peak in radiation
    sw_average = stats_array[2] # average daily peak in radiation
    sw_min = stats_array[3]  # lowest allowable daily peak in radiation
    
    #temp_mean = stats_array[4]
    #temp_anom_norm = stats_array[5]
    
    #hum_mean = stats_array[6]
    #hum_anom_mean = stats_array[7]
    
    # Set to 0.6 in SLAM, these are correlation coefficients, and must be in (0,1). So for
    # consistency:
    
    r_temp = 0.6
    r_hum = 0.6
    r_rad = 0.7
    
    day_length = 15 # IN HOURS
    
    ### CYCLE INITIALIZATION ###
    
    # The basic cycle is formed by taking a unit of free_time, then a unit of day_length,
    # then another unit of free_time. This reconstructs the full 24 hour day, and will be 
    # a useful structure going forward.
    
    free_time = (24 - day_length)*60/2. # In minutes!
    
    free_time = int(free_time)
    
    # Put the time of the day between, say, 9 and 5 and not from 0 to 8.
    
    day_time_shift = np.linspace(free_time, day_length*60 + free_time - 1, day_length*60)
    day_time = np.linspace(1, day_length*60, day_length*60)
    
    # Create radiation profiles as sin waves
    
    day_max = sw_limit * np.sin(np.pi * day_time * (day_length*60)**(-1))
    day_avg = sw_average * np.sin(np.pi * day_time * (day_length*60)**(-1))
    day_min = sw_min * np.sin(np.pi * day_time * (day_length*60)**(-1))
    
    # make "preday" and "postday," i.e., the time that does not include radiation
    
    pre_day = np.linspace(0, free_time - 1, free_time)
    post_day = np.linspace(day_length*60 + free_time, 24*60 - 1, free_time) # i.e., the time from sundown to the end of the day
    
    # make a 'no radiation' array
    
    no_rad = np.zeros(len(pre_day))
    
    # make an array of times that includes the entire day, from midnight to midnight (the next day)
    
    tot_day = np.hstack((pre_day, day_time_shift, post_day))
    tot_day_sw_max = np.hstack((no_rad, day_max, no_rad))
    tot_day_sw_avg = np.hstack((no_rad, day_avg, no_rad))
    tot_day_sw_min = np.hstack((no_rad, day_min, no_rad))
    
    # now we want to extend the above daily framework to the entire summer!
    
    all_summer_time = 92*24*60 # 92 days in summer, 24 hours in a day, 60 mins in an hour
    
    sw_max_sum = np.tile(tot_day_sw_max, 92) # make 92 days worth of the same radiation
    sw_avg_sum = np.tile(tot_day_sw_avg, 92)
    sw_min_sum = np.tile(tot_day_sw_min, 92)
    
    ### CREATE RED NOISE FORCING ###
    
    # characteristic time scale of correlation is ~6 hours
    
    num_sixhr = 92 * 4 # number of days * number of 6 hour periods per day
    
    day_vals = np.linspace(0, 1.5*np.pi, 4)
    day_trend = np.sin(day_vals)
    
    #temp_summer = np.tile(day_trend, 92) # makes an array with 92 copies of temp_day back to back
    #hum_summer = np.tile(day_trend, 92) 
    
    rad_summer = np.tile(day_trend, 92)
    
    i = 1
    while i < num_sixhr:
        
        #temp_rand = random.gauss(0,1) # sample a gaussian distribution with mean = 0 and std = 1
        #hum_rand = random.gauss(0,1) # UNCORRELATED TEMP AND HUM 
        
        rad_rand = random.gauss(0,1)
        
        # add day by day variability; correlate the data
        # factor of (1-r^2)^0.5 makes the variance of the temp distribution the same as the
        # random distribution's that we chose temp_rand from, i.e., is one
        
        #temp_summer[i] = r_temp * temp_summer[i-1] + ((1 - r_temp**2)**(0.5)) * temp_rand
        #hum_summer[i] = r_hum * hum_summer[i-1] + ((1 - r_hum**2)**(0.5)) * hum_rand
        
        rad_summer[i] = r_rad * rad_summer[i-1] + ((1 - r_rad**2)**(0.5)) * rad_rand
        
        i += 1
    
    #temp_summer_min = getInterpolate(temp_mean + temp_summer * temp_anom_norm, all_summer_time)
    #hum_summer_min = getInterpolate(hum_mean + hum_summer * hum_anom_norm, all_summer_time)
    
    rad_summer_min = getInterpolate(rad_summer * rad_anom_norm, all_summer_time)
    
    # THIS FINISHES TEMPERATURE AND HUMIDITY
    
    ### PRECIP AND CLOUD FORCING ###
    
    cf = np.zeros(all_summer_time)
    F_solar = np.zeros(all_summer_time)
    
    cloud_count = 0
    cloud_atlas = [] # tells us where clouds are so we can correlate P and F, i.e., make sure 
    # there isn't a lot of radiation on rainy days
    
    j = 0
    while j < all_summer_time:
        
        rad_val = sw_avg_sum[j] + rad_summer_min[j]
        max_rad_val = sw_max_sum[j]
        min_rad_val = sw_min_sum[j]
        
        # cloud red noise has to be large for the cloud atlas to generate precip
        if rad_summer_min[j] < -rad_anom_norm * 1.25: 
            cloud_atlas.append(j)
            cloud_count += 1
            
        F_solar[j] = rad_val
        
        if F_solar[j] < 0:
            F_solar[j] = 0
        
        if F_solar[j] > max_rad_val:
            F_solar[j] = max_rad_val
            
        if F_solar[j] < min_rad_val and F_solar[j] > 0:
            F_solar[j] = min_rad_val 
        
        """
        cf[j] = - rad_summer_min[j] * (2 * sw_limit)**(-1)
        
        if cf[j] > 1:
            cf[j] = 1
        
        if cf[j] < 0:
            cf[j] = 0
        """
        
        j += 1
        
    precip = np.zeros(all_summer_time)
    precip_array = getPrecip(stats_array)
    
    k = 0
    
    # make sure precip only happens on days where there are clouds 
    while k < len(precip_array):
        rand_cloud_dex = random.randint(0, cloud_count)
        precip_dex = cloud_atlas[rand_cloud_dex]
        precip[precip_dex] = precip_array[k]
        k += 1
    
    return F_solar, precip


# Simple (linear) interpolation function 


def getInterpolate(array, new_length):
    old_time = np.linspace(0, new_length - 1, len(array))
    new_time = np.linspace(0, new_length - 1, new_length)
    new_array = np.interp(new_time, old_time, array)
    
    return new_array


# Get precipitation for red noise generator


def getPrecip(stats_array):
    
    # Take the relevant statistics from the SGP data to tune our precip
    # events
    
    cum_precip = stats_array[4] # average total rainfall in a summer
    std_precip = stats_array[5]
    
    indiv_event_mean = stats_array[6]
    
    # Generate a normal distribution for the cumulative rain?
    
    cumulative_rain = np.random.gamma(cum_precip, scale=1.0) # sample a dist for tot rain
    print(cumulative_rain)
    
    precip_events = []
    
    i = 0
    
    # Confused about this statement. What does it mean for the sum of the 
    # precipitation events to be less than the *distribution* of 
    # cumulative rain events?
    
    while cumulative_rain > np.sum(precip_events):
        tmp_precip = np.random.gamma(indiv_event_mean, scale=1.0) # why from a gamma dist?
        precip_events.append(tmp_precip)
        i += 1 # neat syntax :) 
    
    #print(precip_events)
    
    return precip_events

