#!/usr/bin/env python
# coding: utf-8

# # getStats
# 
# Function takes a list of files (.nc) and returns statistical quantities of interest to tune our forcing. We need to compute the: 
# 
# - radiation anomaly (standard deviation?)
# - maximum, average, and minimum shortwave radiation limits
# - temperature mean
# - temperature anomaly (standard deviation?)
# - humidity mean 
# - humidity anomaly (standard deviation?)
# - precipitation mean
# - precipitation anomaly (standard deviation?)
# - the mean precipitation in any one event

import numpy as np
import xarray as xr
import glob


# function to get the statistics of the relevant file arrays


def getStats(rad_file_array, rad_data_var, precip_file_array, precip_data_var, threshold):
    
    rad_max, rad_max_mean, rad_min, rad_std = getRadStats(rad_file_array, rad_data_var)
    
    #temp_mean, temp_std, temp_max, temp_min = getMeanStdMaxMin(temp_file_array, "temperature")
    
    #hum_mean, hum_std, hum_max, hum_min = getMeanStdMaxMin(hum_file_array, "humidity")

    precip_cum, precip_std = getPrecipStats(precip_file_array, precip_data_var)
    
    indiv_precip_event_mean_array = getPrecipEventMeansXARRAY_1D(precip_file_array, threshold, precip_data_var) # threshold = precip to consider an event
    
    stats_array = np.zeros(7) # initialize stats array
    
    # [rad_anom_norm, sw_limit, sw_average, sw_min, mean_precip, std_precip, indiv_event_mean] in getForcingFunctions
    
    stats_array[0] = rad_std
    stats_array[1] = rad_max
    stats_array[2] = rad_max_mean
    stats_array[3] = rad_min
    
    #stats_array[4] = temp_mean
    #stats_array[5] = temp_std
    
    #stats_array[6] = hum_mean
    #stats_array[7] = hum_std
    
    stats_array[4] = precip_cum
    stats_array[5] = precip_std
    
    stats_array[6] = np.mean(indiv_precip_event_mean_array)
    
    return stats_array

# calc radiation statistics

def getRadStats(rad_file_array, rad_data_var):
    
    # funciton that returns the max radiation peak value, the average radiation peak value, the minimum radiation peak
    # value, and the standard deviation of the defect corrected, daytime diurnal cycles, i.e., the standard deviation
    # of all data that occurs during the day. 
    
    N_files = len(rad_file_array) 
    rad_max_array = []
    rad_diurn_vals = []
    
    for i in range(0, N_files):
        tmp_dataset = xr.open_dataset(rad_file_array[i]) # Open ith summer's dataset
        
        tmp_rad_raw = tmp_dataset[rad_data_var].values # open the "data_var" values for all summers
        
        tmp_rad_fixed = np.convolve(getRemoveDefects(tmp_rad_raw, rad_data_var), np.ones(60)/60, mode='same') # fix values of below zero
        
        # Find highest peak and lowest peak by splitting total data by the day.
        
        N_days = int(len(tmp_rad_raw) * (60*24)**(-1))
        
        #print(N_days)
        
        tmp_rad_daysplit = np.split(tmp_rad_fixed, N_days) # make array of daily arrays
        
        for j in range(0, N_days):
            
            tmp_rad_day = tmp_rad_daysplit[j]
            
            if len(tmp_rad_day) > 0: # check to make sure getRemoveDefects didn't clean every value in the array
                tmp_rad_max = np.amax(tmp_rad_day) # take peak daily value
                tmp_day_vals1 = tmp_rad_day[660:1440] # first part of diurnal cycle
                tmp_day_vals2 = tmp_rad_day[0:60] # second part of diurnal cycle (due to offset in data)
                rad_max_array.append(tmp_rad_max) # append to max array
                rad_diurn_vals.append(tmp_day_vals1)
                rad_diurn_vals.append(tmp_day_vals2)
                
            else: # if getRemoveDefects detected an entire day of bad data, just carry on and don't take the max
                continue 
    
    rad_max = np.amax(rad_max_array) # determine max of each peak
    rad_min = np.amin(rad_max_array) # determine min of each peak
    rad_max_mean = np.mean(rad_max_array) # average peak value

    rad_diurn_vals_con = np.concatenate(rad_diurn_vals)
    rad_std = np.std(rad_diurn_vals_con) # take standard deviation of diurnal cycle values 
    
    return rad_max, rad_max_mean, rad_min, rad_std

# calc precipitation stats

def getPrecipStats(precip_file_array, precip_data_var):
    
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


# function that removes defects in the radiation data, where radiation values are recorded as -9999 W m^{-2} and less than zero
# as of 3/10, this function now linearly interpolates data with bad values instead of nixing them out of the array.

def getRemoveDefects(array, data_var):
    
    # for radiation, the number should always be positive.
    
    if data_var == "BestEstimate_down_short_hemisp":
    
        good_indices = np.where(array > 0)[0]

        if len(good_indices) == 0: # return empty array if found day of bad data

            return []

        else: # if good day of data with defects, clean it up

            interp_array = np.zeros(len(array))

            for i in range(0, len(array)):

                #print("current index is " + str(i))

                if array[i] <= 0:

                    if i == len(array) - 1: # if the defect value is the last value in the array, just set to the last value, since there are no future data points to interpolate to.
                        interp_array[i] = interp_array[i-1]

                    else:

                        #print("value is " + str(array[i]))

                        good_index = np.where(good_indices > i)[0][0]

                        interp_index = good_indices[good_index]

                        #print("interpolating " + str(i) + " to index number " + str(interp_index))

                        interp_line_slope = (array[interp_index] - interp_array[i-1]) * (float(interp_index) - (float(i) - 1))**(-1)

                        interp_value_i = interp_line_slope * (i - (i - 1)) + interp_array[i-1]

                        #print("interp slope is " + str(interp_line_slope))
                        #print("interped value is " + str(interp_value_i))

                        interp_array[i] = interp_value_i

                else:
                    interp_array[i] = array[i]
                   
            
        return interp_array
    
    # zero is a valid input for precip! 
    
    if data_var == "precip":
        
        good_indices = np.where(array >= 0)[0]

        if len(good_indices) == 0: # return empty array if found day of bad data

            return []

        else: # if good day of data with defects, clean it up

            interp_array = np.zeros(len(array))

            for i in range(0, len(array)):

                #print("current index is " + str(i))

                if array[i] < 0:

                    if i == len(array) - 1:
                        interp_array[i] = interp_array[i-1]

                    else:

                        #print("value is " + str(array[i]))

                        good_index_array = np.where(good_indices > i)[0]
                        
                        if len(good_index_array) == 0: # no higher value that is not bad, happens in file 3, day 35, just set to zero.
                            
                            interp_array[i] = 0
                            
                        else:
                            
                            good_index = good_index_array[0]

                            interp_index = good_indices[good_index]

                            #print("interpolating " + str(i) + " to index number " + str(interp_index))

                            interp_line_slope = (array[interp_index] - interp_array[i-1]) * (float(interp_index) - (float(i) - 1))**(-1)

                            interp_value_i = interp_line_slope * (i - (i - 1)) + interp_array[i-1]

                            #print("interp slope is " + str(interp_line_slope))
                            #print("interped value is " + str(interp_value_i))

                            interp_array[i] = interp_value_i

                else:
                    interp_array[i] = array[i]
                   
            
        return interp_array
    
    else:
        return "ERROR: Data keyword not compatible. Try changing the keyword to either 'precip' for precipitation or 'BestEstimate_down_short_hemisp' for radiation."


# get precip event mean for ONE DIMENSIONAL (ONE SUMMER) DATA

def getPrecipEventMeansXARRAY_1D(filearray, threshold, data_var):
    
    N_files = len(filearray) # number of files
    
    #print(N_files)
    
    precip_all_array = []
    event_array = []
    tmp_event = []
    
    # loop through all files and find the following things:
    
    for i in range(0,N_files):
        
        tmp_dataset = xr.open_dataset(filearray[i]) # open a file's worth of data
        
        tmp_precip_data = tmp_dataset[data_var].values # extract precip values as a numpy array
        
        tmp_precip = np.zeros(tmp_precip_data.shape[0])
        
        if data_var == "precip_syn":
            
            for k in range(0, tmp_precip_data.shape[0]):
                tmp_precip[k] = tmp_precip_data[k][0]
                
        else:
            tmp_precip = tmp_precip_data
    
        #print(tmp_precip)
            
        index_array = np.where(tmp_precip > threshold)[0] # find within this numpy array where an event has occured relative to some threshold

        N_mins = len(tmp_precip)
            
        #print(N_mins)
        
        for j in range(0, N_mins): 
                
            # loop through minute indices. if index j is in index array, then precip event happened at
            # that index. add those precip values to a temporary event array. 

            if j in index_array:
                tmp_event.append(float(tmp_precip[j])) # append precipitation amount for an event

                if j+1 in index_array: # if the next member of the array is also part of the event, keep looping
                    continue

                else:
                    event_array.append(tmp_event) # if the next member is not part of the event, append the event to a list of all events
                    tmp_event = [] # clear temp event array

            else:
                continue 
                    
    
    N_events = len(event_array)
    mean_array = []
    
    for i in range(0, N_events):
        tmp_mean = np.mean(event_array[i])
        mean_array.append(tmp_mean) # take mean of individual events and add them to an array
    
    return mean_array # array of mean of all events


# get half the autocorrelation array that numpy computes 

def getAutocorrelationTrunc(array):
    result = np.correlate(array, array, mode='full')
    return result[int(result.size/2):] # makes t \in (0,infty) instead of t \in (-infty, infty)

# get daily mean and max (can turn on std and skew arrays, but rn they're turned off)

def getDailyMeanMaxArrays(filelist, keyword):
    
    N_files = len(filelist)
            
    daily_mean_array = [] # the mean of each day's radiation
    daily_max_array = [] # the maximum value of the each day's radiation
    #daily_std_array = [] # the standard deviation of each day's diurnal cycle
    #daily_skew_array = [] # the skewness of each day's diurnal cycle
    
    for i in range(0, N_files):
        
        print(i) # what file are we on? 
        
        tmp_ds = xr.open_dataset(filelist[i])
        
        tmp_rad_data = tmp_ds[keyword].values
        
        N_days = int(tmp_rad_data.shape[0] * (60*24)**(-1))
        
        tmp_rad = np.zeros(N_days*60*24) 
        
        if keyword == "F_solar": # keyword check to modify how this function works for synthetic data vs. obs data
            
            for k in range(0, tmp_rad_data.shape[0]): 
                tmp_rad[k] = tmp_rad_data[k][0]
        
        else:
            tmp_rad = tmp_rad_data
        
        #print(tmp_rad)
        
        N_days = int(len(tmp_rad) * (60*24)**(-1))
        
        #print(N_days)
        
        tmp_rad_daysplit = np.split(tmp_rad, N_days)
        
        for j in range(0, N_days):
            
            #print(j)
            
            if keyword == "BestEstimate_down_short_hemisp":
                tmp_rad_day_clean = getRemoveDefects(tmp_rad_daysplit[j], keyword) # day's worth of clean data from obs
            
            else:
                tmp_rad_day_clean = tmp_rad_daysplit[j]
            
            if len(tmp_rad_day_clean) > 0: # if getRemoveDefects didn't remove every data point, do:
                tmp1 = tmp_rad_day_clean[0:60]
                tmp2 = tmp_rad_day_clean[660:1440]
                tmp_diurn = np.hstack((tmp1,tmp2))
                daily_mean_array.append(np.mean(tmp_diurn, axis=0)) # compute mean parameter value for that day
                daily_max_array.append(np.amax(tmp_rad_day_clean, axis=0)) # report maximum param value for that day
                #daily_std_array.append(np.std(tmp_rad_day_clean, axis=0))
                #daily_skew_array.append(sts.skew(tmp_rad_day_clean, axis=0))
                
            else: # if getRemoveDefects found an entire day's worth of bad data, skip that day and continue looping
                continue
        
    #return daily_mean_array, daily_max_array, daily_std_array, daily_skew_array
    return daily_mean_array, daily_max_array


def getDailyMeanMaxArrayssyn(filelist, keyword):
    
    N_files = len(filelist)
            
    daily_mean_array = [] # the mean of each day's radiation
    daily_max_array = [] # the maximum value of the each day's radiation
    #daily_std_array = [] # the standard deviation of each day's diurnal cycle
    #daily_skew_array = [] # the skewness of each day's diurnal cycle
    
    for i in range(0, N_files):
        
        print(i) # what file are we on? 
        
        tmp_ds = xr.open_dataset(filelist[i])
        
        tmp_rad_data = tmp_ds[keyword].values
        
        N_days = int(tmp_rad_data.shape[0] * (60*24)**(-1))
        
        tmp_rad = np.zeros(N_days*60*24) 
        
        if keyword == "F_solar": # keyword check to modify how this function works for synthetic data vs. obs data
            
            for k in range(0, tmp_rad_data.shape[0]): 
                tmp_rad[k] = tmp_rad_data[k][0]
        
        else:
            tmp_rad = tmp_rad_data
        
        #print(tmp_rad)
        
        N_days = int(len(tmp_rad) * (60*24)**(-1))
        
        #print(N_days)
        
        tmp_rad_daysplit = np.split(tmp_rad, N_days)
        
        for j in range(0, N_days):
            
            #print(j)
            
            if keyword == "BestEstimate_down_short_hemisp":
                tmp_rad_day_clean = getRemoveDefects(tmp_rad_daysplit[j], keyword) # day's worth of clean data from obs
            
            else:
                tmp_rad_day_clean = tmp_rad_daysplit[j]
            
            if len(tmp_rad_day_clean) > 0: # if getRemoveDefects didn't remove every data point, do:
                tmp_diurn = tmp_rad_day_clean[360:1080]
                daily_mean_array.append(np.mean(tmp_diurn, axis=0)) # compute mean parameter value for that day
                daily_max_array.append(np.amax(tmp_rad_day_clean, axis=0)) # report maximum param value for that day
                #daily_std_array.append(np.std(tmp_rad_day_clean, axis=0))
                #daily_skew_array.append(sts.skew(tmp_rad_day_clean, axis=0))
                
            else: # if getRemoveDefects found an entire day's worth of bad data, skip that day and continue looping
                continue
        
    #return daily_mean_array, daily_max_array, daily_std_array, daily_skew_array
    return daily_mean_array, daily_max_array

# N dimensional generalization of getPrecipEventMeansXARRAY_1D, i.e., if the files have more than one dimension to them (multiple summers, etc...)

def getPrecipEventMeansXARRAY_ND(filearray, threshold):
    
    N_files = len(filearray) # number of files
    
    print(N_files)
    
    precip_all_array = []
    event_array = []
    tmp_event = []
    
    # loop through all files and find the following things:
    
    for i in range(0,N_files):
        
        tmp_dataset = xr.open_dataset(filearray[i]) # open a file's worth of data
        
        tmp_precip = tmp_dataset["precip"].values # extract precip values as a numpy array EDIT KEYWORD WHEN YOU KNOW MORE
    
        #print(tmp_precip)
        
        N_sum = tmp_precip.shape[1] # number of summers
        
        print(N_sum)
        
        for k in range(0, N_sum): # loop through the summers
            
            index_array = np.where(tmp_precip[:,k] > threshold)[0] # find within this numpy array where an event has occured relative to some threshold

            N_mins = tmp_precip.shape[0]
            
            print(N_mins)
        
            for j in range(0, N_mins): 
                
                # loop through minute indices. if index j is in index array, then precip event happened at
                # that index. add those precip values to a temporary event array. 

                if j in index_array:
                    tmp_event.append(float(tmp_precip[j,k])) # append precipitation amount for an event

                    if j+1 in index_array: # if the next member of the array is also part of the event, keep looping
                        continue

                    else:
                        event_array.append(tmp_event) # if the next member is not part of the event, append the event to a list of all events
                        tmp_event = [] # clear temp event array

                else:
                    continue 
                    
    
    N_events = len(event_array)
    mean_array = []
    
    for i in range(0, N_events):
        tmp_mean = np.mean(event_array[i])
        mean_array.append(tmp_mean) # take mean of individual events and add them to an array
    
    return np.mean(mean_array) # mean of all events ("mean of the mean of individual events")

"""
# testing bits of this code 

rad_filelist=glob.glob("/data/keeling/a/adammb4/SGP_proj_2021/DATA/*_f.nc")
precip_filelist = glob.glob("/data/keeling/a/adammb4/SGP_proj_2021/DATA/*_p.nc")
threshold = 0.1


stats_array = getStats(rad_filelist, precip_filelist, threshold)


print(stats_array)


# These were used to make test1_f and test1_p. left here to check by hand calcs done by code above

F_solar_1 = [0,1,1,2,1,0]
f2 = [0,2,3,2,1,2]

p1 = [0,1,2,0,1,0]
p2 = [0,1,2,3,0,1]

F_solar = np.stack((F_solar_1, f2)).T
precip = np.stack((p1,p2)).T


# These were used to make test2_f and test2_p. left here to check by hand calcs done by code above

F1 = [0,4,5,0,0,0]
F2 = [0,1,3,0,0,1]

P1 = [0,.4,2,0,2,0]
P2 = [1,1,4,3,0,0]

F_solar2 = np.stack((F1, F2)).T
precip2 = np.stack((P1,P2)).T
"""