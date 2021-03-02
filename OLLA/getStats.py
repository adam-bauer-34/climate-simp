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


def getStats(rad_file_array, precip_file_array, threshold):
    
    rad_mean, rad_std, rad_max, rad_min = getMeanStdMaxMin(rad_file_array, "BestEstimate_down_short_hemisp") 
    
    #temp_mean, temp_std, temp_max, temp_min = getMeanStdMaxMin(temp_file_array, "temperature")
    
    #hum_mean, hum_std, hum_max, hum_min = getMeanStdMaxMin(hum_file_array, "humidity")
    
    precip_mean, precip_std, precip_max, precip_min = getMeanStdMaxMin(precip_file_array, "precip")
    
    indiv_precip_event_mean = getPrecipEventMeansXARRAY_1D(precip_file_array, threshold) # threshold = precip to consider an event
    
    ###
    stats_array = np.zeros(7) # initialize stats array
    
    # [rad_anom_norm, sw_limit, sw_average, sw_min, temp_mean, temp_anom_norm, hum_mean, hum_anom_mean, mean_precip, std_precip, indiv_event_mean] in getForcingFunctions
    
    stats_array[0] = rad_std
    stats_array[1] = rad_max
    stats_array[2] = rad_mean
    stats_array[3] = rad_min
    
    #stats_array[4] = temp_mean
    #stats_array[5] = temp_std
    
    #stats_array[6] = hum_mean
    #stats_array[7] = hum_std
    
    stats_array[4] = precip_mean
    stats_array[5] = precip_std
    
    stats_array[6] = indiv_precip_event_mean
    
    return stats_array


# function to calculate the maximum, minimum, mean, and standard deviation of a set of data


def getMeanStdMaxMin(filearray, data_var):
    
    N = len(filearray) 
    param_all_array = []
    param_max_array = []
    param_min_array = []
    
    for i in range(0, N):
        tmp_dataset = xr.open_dataset(filearray[i]) # Open ith summer's dataset
        
        tmp_param = tmp_dataset[data_var].values # open the "data_var" values for all summers
        
        tmp_param_fixed = getRemoveDefects(tmp_param) # fix radiation values of -1000 and below
        
        #print(tmp_param_fixed)
        
        tmp_param_max = np.amax(tmp_param_fixed) # take the file max
        tmp_param_min = np.amin(tmp_param_fixed) # take the file min
        
        param_max_array.append(tmp_param_max) # add file max to list
        param_min_array.append(tmp_param_min) # add file min to list
        param_all_array.append(tmp_param_fixed) # add values to total array
    
    param_all_array = np.concatenate(param_all_array) # concatenate arrays
    
    param_max = np.amax(param_max_array) # determine max of all files
    param_min = np.amin(param_min_array) # determine min of all files
    
    param_mean = np.mean(param_all_array) # determine mean of all files
    param_std = np.std(param_all_array) # determine standard deviation of all files
    
    return param_mean, param_std, param_max, param_min


# function that gets the mean of the individual precipitation events FOR MULTIPLE SUMMER DATA!!!!

def getRemoveDefects(array):
    
    defect_indices = np.where(array < 0)[0]
    
    new_array = []
    
    for i in range(0, len(array)):
        
        if i in defect_indices:
            continue
            
        else:
            new_array.append(array[i])
    
    return new_array
    

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


# get precip event mean for ONE DIMENSIONAL (ONE SUMMER) DATA

def getPrecipEventMeansXARRAY_1D(filearray, threshold):
    
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
            
        index_array = np.where(tmp_precip > threshold)[0] # find within this numpy array where an event has occured relative to some threshold

        N_mins = len(tmp_precip)
            
        print(N_mins)
        
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

