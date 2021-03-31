#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Test version of getPrecipEventMeansXARRAY

def getPrecipEventMeans(precip):
    
    precip_all_array = []

    index_array = getEventIndicies(precip)
    
    N = len(precip)
    event_array = []
    tmp_event = []
    
    for i in range(0, N):
        
        if i in index_array:
            tmp_event.append(float(precip[i]))
            
            if i+1 in index_array:
                continue
            
            else:
                event_array.append(tmp_event)
                tmp_event = []
            
        else:
            continue
            
    N_events = len(event_array)
    mean_array = []
    
    for i in range(0, N_events):
        tmp_mean = np.mean(event_array[i])
        mean_array.append(tmp_mean)
    
            
    return np.mean(mean_array)


# In[ ]:


# Test version of getMeanStdMaxMin

def getTestStats(filearray):
    
    N = len(filearray)
    Tup_all_array = []
    Tup_max_array = []
    
    for i in range(0, N):
        tmp_dataset = xr.open_dataset(filearray[i])
        
        tmp_Tup = tmp_dataset["c_f"].values
        tmp_Tup_max = np.amax(tmp_Tup)
        
        Tup_max_array.append(tmp_Tup_max)
        Tup_all_array.append(tmp_Tup)
    
    Tup_all_array = np.concatenate(Tup_all_array)
    
    Tup_max = np.amax(Tup_max_array)
    
    Tup_mean = np.mean(Tup_all_array)
    Tup_std = np.std(Tup_all_array)
    
    return Tup_mean, Tup_std, Tup_max


# In[ ]:


# Newtonian integrator proof of concept

import numpy as np
import matplotlib.pyplot as plt

y = np.zeros(1000)

x = np.linspace(0, 10, 1000)

for i in range(0, len(y)):
    y[i] = y[i-1] + x[i-1] * .01
    
plt.plot(x, y, color='r')
plt.plot(x, 0.5*x**2, color='b')

# function to calculate the maximum, minimum, mean, and standard deviation of a set of data. antequated after sufficent differences in desired values for precip and radiaiton were discovered.


def getMeanStdMaxMin(filearray, data_var):
    
    N_files = len(filearray) 
    param_all_array = []
    param_max_array = []
    param_cum_array = []
    
    for i in range(0, N_files):
        tmp_dataset = xr.open_dataset(filearray[i]) # Open ith summer's dataset
        
        tmp_param = tmp_dataset[data_var].values # open the "data_var" values for all summers
        
        tmp_param_fixed = getRemoveDefects(tmp_param, data_var) # fix values of below zero
        
        param_all_array.append(tmp_param_fixed) # add values to total array
        param_cum_array.append(np.sum(tmp_param_fixed)) # sum of all values in array
        
        # Find highest peak and lowest peak by splitting total data by the day.
        
        N_days = int(len(tmp_param) * (60*24)**(-1))
        
        #print(N_days)
        
        tmp_param_daysplit = np.split(tmp_param, N_days) # make array of daily arrays
        
        for j in range(0, len(tmp_param_daysplit)):
            
            #print(i,j)

            tmp_daysplit_clean = getRemoveDefects(tmp_param_daysplit[j], data_var) # clean indiv day
            
            if len(tmp_daysplit_clean) > 0: # check to make sure getRemoveDefects didn't clean every value in the array
                tmp_param_max = np.amax(tmp_daysplit_clean) # take peak daily value
                param_max_array.append(tmp_param_max) # append to max array
                
            else: # if getRemoveDefects detected an entire day of bad data, just carry on and don't take the max
                continue 
    
    param_all_array = np.concatenate(param_all_array) # concatenate arrays
    
    param_max = np.amax(param_max_array) # determine max of each peak
    param_min = np.amin(param_max_array) # determine min of each peak
    param_max_mean = np.mean(param_max_array) # average peak value
    
    param_mean = np.mean(param_all_array) # determine mean of all files
    param_std = np.std(param_all_array) # determine standard deviation of all files
    param_cum = np.mean(param_cum_array) # CUMULATIVE sum of all the data ....... GeT uR hEaD oUt Of ThE gUtTeRs
    
    return param_max_mean, param_std, param_max, param_min, param_cum

# calc precipitation stats

def getPrecipStats_OLD(precip_file_array, precip_data_var):
    
    N_files = len(precip_file_array) 
    precip_cum_array = []
    precip_all_array = []
    
    for i in range(0, N_files):
        tmp_dataset = xr.open_dataset(precip_file_array[i]) # Open ith summer's dataset
        
        tmp_precip = tmp_dataset[precip_data_var].values # open the "data_var" values for all summers
        
        tmp_precip_fixed = getRemoveDefects(tmp_precip, precip_data_var) # fix values of below zero
        
        precip_all_array.append(tmp_precip_fixed) # add values to total array
        precip_cum_array.append(np.sum(tmp_precip_fixed)) # sum of all values in array
    
    precip_all_array_con = np.concatenate(precip_all_array) # concatenate arrays
    
    precip_std = np.std(precip_all_array_con) # determine standard deviation of all files
    precip_cum = np.mean(precip_cum_array) # CUMULATIVE sum of all the data ....... GeT uR hEaD oUt Of ThE gUtTeRs
    
    return precip_cum, precip_std

# old red noise function

def getRedNoise(stats_array):
    
    ### PARAMETER SETUP ###
    
    # Some parameters of getRedNoise are fed in from the SGP data, and thus are stored in 
    # stats_array.
    
    #rad_anom_norm = 75.
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
    
    day_length = 13 # IN HOURS
    
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
        rand_cloud_dex = random.randint(0, cloud_count-1)
        precip_dex = cloud_atlas[rand_cloud_dex]
        precip[precip_dex] = precip_array[k]
        k += 1
    
    return F_solar, precip

# Get precipitation for red noise generator, no longer used

def getPrecip(stats_array):
    
    # Take the relevant statistics from the SGP data to tune our precip
    # events
    
    cum_precip = stats_array[4] # average total rainfall in a summer
    std_precip = stats_array[5]
    
    indiv_event_mean = stats_array[6]
    
    # Generate a normal distribution for the cumulative rain?
    
    cumulative_rain = np.random.gamma(cum_precip, scale=1.0) # sample a dist for tot rain
    #print(cumulative_rain)
    
    precip_events = []

    # Confused about this statement. What does it mean for the sum of the 
    # precipitation events to be less than the *distribution* of 
    # cumulative rain events?
    
    while cumulative_rain > np.sum(precip_events):
        tmp_precip = np.random.gamma(indiv_event_mean, scale=.6) # why from a gamma dist?
        if tmp_precip < 0.2:
            continue 
        
        elif tmp_precip > 2.5:
            continue 
            
        else:
            precip_events.append(tmp_precip)
    
    #print(precip_events)
    
    return precip_events