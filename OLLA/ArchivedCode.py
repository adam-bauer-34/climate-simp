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