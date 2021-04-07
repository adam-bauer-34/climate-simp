#!/usr/bin/env python
# coding: utf-8

# # Master OLLA Integrator
# 
# This file has the integrator and actual OLLA processing code. OLLA takes a .nc file of forcing and outputs the temperature time series. 


import numpy as np 
import xarray as xr 
from getStats import getRemoveDefects


# Get integration -- i.e., this is the model! 


def getOLLASyn(filename):
    
    forcing = xr.open_dataset(filename[0])
    
    F_solar = forcing["F_solar"].values
    #t_upper = forcing["t_upper"].values
    #q_upper = forcing["q_upper"].values
    precip = forcing["precip_syn"].values
    #cf = forcing["c_f"].values
    
    ### NEED A FAT UNITS CHECK ON THE ABOVE DATA TO MAKE SURE WE'RE NOT MESSING UP UNITS ###
    
    temp_series, mois_series = getIntegrationSyn(F_solar, precip)
    
    return temp_series, mois_series

def getOLLAObs(radfile, precipfile, sumnum):
    
    radforcing = xr.open_dataset(radfile[sumnum])
    precipforcing = xr.open_dataset(precipfile[sumnum])
    
    F_solar_raw = radforcing["BestEstimate_down_short_hemisp"].values
    F_solar_fixed = getRemoveDefects(F_solar_raw,"BestEstimate_down_short_hemisp")
    F_solar_clean = np.convolve(F_solar_fixed, np.ones(60)/60, mode='same') # smooth data on hour window
    
    precip = getRemoveDefects(precipforcing["precip"].values, "precip")
    
    ### NEED A FAT UNITS CHECK ON THE ABOVE DATA TO MAKE SURE WE'RE NOT MESSING UP UNITS ###
    
    temp_series, mois_series = getIntegrationObs(F_solar_clean, precip)
    
    return temp_series, mois_series


# Newtonian integrator for model


def getIntegrationSyn(F_solar, precip):
    
    N_minutes = F_solar.shape[0]
    N_summers = F_solar.shape[1]
    
    temp_series = np.zeros((N_minutes, N_summers))
    mois_series = np.zeros((N_minutes, N_summers))
    
    for i in range(0, N_summers): # loop through summers
        
        for j in range(0, N_minutes): # loop through minutes
            
            temp_series[j,i] = temp_series[j-1,i] + getTempFlux(temp_series[j-1,i], mois_series[j-1,i], F_solar[j-1,i])
            mois_series[j,i] = mois_series[j-1,i] + getMoisFlux(temp_series[j-1,i], mois_series[j-1,i], precip[j-1,i])
    
    return temp_series, mois_series

# for obs data 
def getIntegrationObs(F_solar, precip):
    
    N_minutes = F_solar.shape[0]
    
    temp_series = np.zeros(N_minutes)
    mois_series = np.zeros(N_minutes)
    

    for j in range(0, N_minutes): # loop through minutes

        temp_series[j] = temp_series[j-1] + getTempFlux(temp_series[j-1], mois_series[j-1], F_solar[j-1])
        mois_series[j] = mois_series[j-1] + getMoisFlux(temp_series[j-1], mois_series[j-1], precip[j-1])
    
    return temp_series, mois_series

# Get temperature flux for the integration


def getTempFlux(T, m, F):
    
    L = 2.5 * 10**6 # latent heat of evaporation
    T_0 = 273.15 # freezing temperature
    R_w = 461.25 # specific gas constant for water 
    e_s_0 = 6.11 # hPa; constant in front of clausius clapeyron relation
    P_surf = 1013.25 #hPa; surface pressure 
    C = 4180. # heat capacity of water
    q = 0.01
    
    alpha = 10. # radiative feedback
    v = 10**(-2) # density of air divided by surface resistance
    
    exp_arg = L * R_w**(-1) * (T_0**(-1) - (T+305.)**(-1))
    #print(exp_arg)
    sat_humidity = 0.622 * e_s_0 * P_surf**(-1) * np.exp(exp_arg)
    
    VPD = sat_humidity - q
    
    if VPD < 0:
        VPD = 0
    
    temp_flux = C**(-1) * (F - alpha * (T - T_0) - L * v * m * VPD)
    
    return temp_flux


# Get moisture flux for integration


def getMoisFlux(T, m, precip):
    
    L = 2.5 * 10**6 # latent heat of evaporation
    T_0 = 273.15 # freezing temperature
    R_w = 461.25 # specific gas constant for water 
    e_s_0 = 6.11 # hPa; constant in front of clausius clapeyron relation
    P_surf = 1013.25 #hPa; surface pressure 
    mu = 50. # density of water * soil column depth * porisity of soil
    q = 0.01
    
    alpha = 10. # radiative feedback
    v = 10**(-2) # density of air divided by surface resistance
    
    exp_arg = L * R_w**(-1) * (T_0**(-1) - (T+305.)**(-1))
    sat_humidity = 0.622 * e_s_0 * P_surf**(-1) * np.exp(exp_arg)
    
    VPD = sat_humidity - q
    
    if VPD < 0:
        VPD = 0
    
    mois_flux = (precip - v * m * VPD) * mu**(-1)
    
    return mois_flux
    #return 0