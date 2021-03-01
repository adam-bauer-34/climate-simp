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

