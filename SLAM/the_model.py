#! /usr/bin/python

# SLAM(...,...) contains the model code itself and can be modified to return any energy/moisture flux or state varaible that your heart desires. The only other operative function here is q_s(T), which takes a temperature value in Kelvin and spits out an approximation of the saturation specific humidity. There are two other functions at the bottom that have to do with moisture fluxes through the soil, but this version of the model code doesn't use them. 

import numpy as np
import math as m

def q_s(T):
	
	# Function to calculate saturation mixing ratio

	T = T - 273.15
	
	e_s = 6.11*10**(7.5*T/(237.5+T))  # TEMPERATURES IN CELSIUS!!!

	return(e_s*0.61/1000.)  # assuming surface pressure = 1000 hPa and in units of kg/kg

def SLAM(F_solar,T_up,q_up,precip_ar,cloud_frac):

############ CONSTANTS ###########################################


# Basic Physical Constants and model geometry

	sig   = 5.67e-08 	# SB constant (W/(m^2K^4)) from Wikipedia
	rho_a = 1.225 		# density of air (kg/m^3) from Wikipedia
	rho_l = 1000.    	# density of liquid water (kg/m^3)
	hc    = 100.            # height of canopy (m)
	H     = 1000.           # height of upper layer
	d_s   = .1		# depth of surface layer (m)
	d_smm = 100		# depth of surface layer (mm)
	D     = 1.              # depth of overall variable land surface layer
	Lv    = 2257000.        # latent heat of vaporization (J/kg) from Wikipedia
	c_pa  = 1003.           # heat capacity of air (J/(kg K)) from Ohio Mechanical Engineering
	dt    = 60.		# time step (sec)
	T_abyss = 280		# bottom boundary condition (K)

# Soil thermal properties 

	rho_B = 900.    	# bulk density of dry soil including air filled pore space (kg/m^3)
	rho_l = 1000.		# density of water (kg/m^3)
	c_B  = 900.	    	# heat capacity of dry soil including air filled pore space (J/(kg K)) 
	c_l  = 4185		# heat capacity of water (J/kg K)
	K_Tmax = 2.        	# Thermal conductivity in W/(m K) 

	theta_ssat = .4   	# ceiling on m_s
	theta_rsat = .4  	# ceiling on m_r
	min_theta_s = 0.1 	# For index calculation (unitless)
	min_theta_r = 0.1	# """"""""""""""""""""""""""""""""

# More parameters

	f_v = .2		# vegetation fraction 
	f_r = .5		# deep root fraction
	r_aT_max = 120		# maximum turbulent resistance for heat fluxes (s/m)
	r_aq_max = 45		# '' 			        '' moisture fluxes (s/m)
	sigma = 1000 		# peak insolation for buoyancy correction (W/m^2)
	r_s = 75		# surface resistance (s/m)
	r_st = 75		# stomatal resistance (s/m)

# Longwave

	e_o = 0.5		# first coefficient in longwave expansion
	e_1 = 0.1		# second coefficient in longwave expansion
	e_cloud = 0.4		# cloud component of longwave expansion

############## INITIALIZING ARRAYS #############################################################

	time = np.linspace(0,len(F_solar) - 1,len(F_solar))
	times = len(time)

# Dynamic Variables

	T_r = np.zeros(times)
	T_s = np.zeros(times)
	T_1 = np.zeros(times)
	T_2 = np.zeros(times)
	q_1 = np.zeros(times)
	q_2 = np.zeros(times)
	theta_s = np.zeros(times)
	theta_r = np.zeros(times)
	
# Fluxes of interest 

	Kout = np.zeros(times)
	Ksr  = np.zeros(times)
	Msr  = np.zeros(times)
	Wsr  = np.zeros(times)
	Hs   = np.zeros(times)
	Es   = np.zeros(times)
	Tr   = np.zeros(times)
	Ts   = np.zeros(times)
	lw1  = np.zeros(times)
	lw2  = np.zeros(times)
	lws  = np.zeros(times)
	Hatm = np.zeros(times)
	qatm = np.zeros(times)
	Hout = np.zeros(times) 
	qout = np.zeros(times)
	I = np.zeros(times)
	Drain = np.zeros(times)

	ep_1 = np.zeros(times)
	ep_2 = np.zeros(times)
	ep_up = np.zeros(times)

# Soil Thermal Params

	rho_s = np.zeros(times)
	rho_r = np.zeros(times)

	c_s   = np.zeros(times)
	c_r   = np.zeros(times)

# Initial Conditions

	q_1[0] = .015   	# kg/kg
	q_2[0] = .015   	# kg/kg

	theta_s[0] = 0.1 	# m^3/m^3
	theta_r[0] = 0.1	# m^3/m^3
		
	T_r[0] = 283.  		# Kelvin 
	T_s[0] = 15+273.15 	# Kelvin 
	T_1[0] = 20+273.15 	# Kelvin
	T_2[0] = 20+273.15 	# Kelvin

	rho_s[0] = rho_B + theta_s[0]*rho_l
	rho_r[0] = rho_B + theta_r[0]*rho_l
	
	c_s[0] = c_B + theta_s[0]*rho_B*(c_l - c_B)/rho_l
	c_r[0] = c_B + theta_r[0]*rho_B*(c_l - c_B)/rho_l

############################### THE BIG LOOOOOOOOOOP ##############################################

	i = 0

	while i < times - 1:

	# EMISSIVITY AS A FUNCTION OF WATER VAPOR

		ep_2[i] = e_o + e_1*m.log10(q_2[i]*1000) 		# unitless
		ep_1[i] = e_o + e_1*m.log10(q_1[i]*1000) 		# unitless

	# CLOUD FORCING ON UPPER LEVEL EMISSIVITY

		ep_up[i] = e_o + e_1*m.log10(q_up[i]*1000) + e_cloud*cloud_frac[i] # unitless

		if ep_up[i] > 1:
			ep_up[i] = 1

	# WATER AVAILABILITY FUNCTIONS/SOIL WETNESS INDEX
	# CURRENTLY LINEAR, but these can be messed around with to your heart's content

		A = (theta_s[i] - min_theta_s)/(theta_ssat-min_theta_s) # unitless
		B = (theta_r[i] - min_theta_r)/(theta_rsat-min_theta_r) # unitless

	# THERMAL CONDUCTIVITY
	
		K_T_up = K_Tmax*np.log((m.e - 1)*(A+B)/2. + 1)  	# W/(m K)
		K_T_down = K_Tmax*np.log((m.e - 1)*B + 1) 		# W/(m K)

	# TURBULENT ENERGY FLUXES in W/m^2

		if T_s[i] > T_1[i]:
			Hs[i] = rho_a*c_pa*(T_s[i] - T_1[i])/(r_s)  	# W/m^2 

	# Conductive Heat Fluxes

		Ksr[i] = (T_s[i] - T_r[i])*(K_T_up/d_s)  		# W/m^2
		Kout[i] = (T_r[i] - T_abyss)*(K_T_down/(D-d_s)) 	# W/m^2

	# SATURATION SPECIFIC HUMIDITY

		qss = q_s(T_s[i]) 					# kg/kg
		qs1 = q_s(T_1[i])    					# kg/kg

	# EVAPORATION AND TRANSPIRATION BEFORE ACCOUNTING FOR MOISUTRE AVAILABILITY

		if qss > q_1[i]:
		
			Es_ideal = (1-f_v)*rho_a*(qss - q_1[i])/r_s # in kg/(m^2 sec)  Kg OF LIQUID WATER!

		else:
			Es_ideal = 0

		if qs1 > q_1[i]:
		
			Ec_ideal = f_v*rho_a*(qs1 - q_1[i])/r_st # in kg/(m^2 sec)   Kg OF LIQUID WATER
			
		else:
			Ec_ideal = 0
	

	# TOTAL EVAPOTRANSPIRATION

		ET = A*Es_ideal + ((1-f_r)*A + f_r*B)*Ec_ideal # in kg/(m^2 sec)
		Es[i] = A*Es_ideal # in kg/(m^2 sec)
		Ts[i] = A*(1-f_r)*Ec_ideal # in kg/(m^2 sec)
		Tr[i] = f_r*B*Ec_ideal # in kg/(m^2 sec)

#		delta_P = Lv*(Ts[i] + Tr[i]) - c_s[i]*T_s[i]*Ts[i] - c_r[i]*T_r[i]*Tr[i]

	# Boundary Layer Turbulence

		Total_flux = Hs[i] + Lv*(Tr[i] + Es[i] + Ts[i]) # in W/m^2

		r_aT_onetwo = r_aT_max*(1 - Total_flux/sigma)  # in s/m
		r_aq_onetwo = r_aq_max*(1 - Total_flux/sigma)  # in s/m

		Hatm[i] = rho_a*c_pa*(T_1[i] - T_2[i])/r_aT_onetwo # in W/m^2
		qatm[i] = rho_a*(q_1[i] - q_2[i])/r_aq_onetwo # in kg/m^2s

		Hout[i] = rho_a*c_pa*(T_2[i] - T_up[i])/r_aT_onetwo  # in W/m^2
		qout[i] = rho_a*(q_2[i] - q_up[i])/r_aq_onetwo       # in kg/(m^2 s)

		
	# LONGWAVE	

		lws[i] = sig*(-T_s[i]**4 + ep_1[i]*(T_1[i]**4) + (1-ep_1[i])*ep_2[i]*T_2[i]**4 + (1-ep_1[i])*(1-ep_2[i])*ep_up[i]*(T_up[i]**4))  	# in W/m^2 

		lw1[i] = ep_1[i]*sig*(T_s[i]**4 - 2*T_1[i]**4 + ep_2[i]*T_2[i]**4 + (1-ep_2[i])*ep_up[i]*(T_up[i])**4)  				# in W/m^2

		lw2[i] = ep_2[i]*sig*(ep_1[i]*T_1[i]**4 - 2*T_2[i]**4 + ep_up[i]*(T_up[i])**4 + (1-ep_1[i])*T_s[i]**4) 					# in W/m^2

	# VAPOR STEPPING

		dq1 = (1/(hc*rho_a))*(Es[i] + Tr[i] + Ts[i] - qatm[i]) 		# kg/kg per second

		dq2 = (1/((H-hc)*rho_a))*(qatm[i] - qout[i]) 			# kg/kg per second

		q_1[i+1] = q_1[i] + dq1*dt   					# in kg/kg

		q_2[i+1] = q_2[i] + dq2*dt   					# in kg/kg

	# LIQUID WATER STEPPING

		dtheta_s = (1/(d_s*rho_l))*(-Es[i] - Ts[i]) 				# in 1/s
	
		dtheta_r = (1/((D-d_s)*rho_l))*(-Tr[i])			 		# in 1/s

		theta_s[i+1] = precip_ar[i]/d_smm + dtheta_s*dt + theta_s[i]    # in m^3/m^3

		theta_r[i+1] = dtheta_r*dt + theta_r[i]				# in m^3/m^3

	# DRAINAGE AND INFILTRATION

		if theta_s[i+1] > theta_ssat:
			I[i] = (theta_s[i+1] - theta_ssat)*d_s*rho_l/dt  	# in kg/m^2s
			theta_r[i+1] = theta_r[i+1] + dt*I[i]/(rho_l*(D-d_s))	# in m^3/m^3
			theta_s[i+1] = theta_ssat				# in m^3/m^3 

		if theta_r[i+1] > theta_rsat:
			Drain[i] = (theta_r[i+1] - theta_rsat)*(D-d_s)*rho_l/dt # in kg/m^2s
			theta_r[i+1] = theta_rsat				# in m^3/m^3

	# TEMPERATURE STEPPING

		dTr = (1/(rho_r[i]*(D-d_s)*c_r[i]))*(Ksr[i] - Kout[i])						# in K/s

		dTs = (1/(rho_s[i]*d_s*c_s[i]))*(F_solar[i] + lws[i] - Hs[i] - Lv*Es[i] - Ksr[i]) 		# in K/s

		dT1 = (1/(rho_a*hc*c_pa))*(lw1[i] + Hs[i] - Hatm[i])						# in K/s

		dT2 = (1/(rho_a*(H-hc)*c_pa))*(lw2[i] + Hatm[i] - Hout[i])					# in K/s

		T_s[i+1] = dTs*dt + T_s[i]									# in K
	
		T_1[i+1] = dT1*dt + T_1[i]									# in K
	
		T_2[i+1] = dT2*dt + T_2[i]									# in K

		T_r[i+1] = dTr*dt + T_r[i]									# in K


	# SOIL HEAT CAPACITY AND DENSITY AS A FUNCTION OF THETA

		rho_s[i+1] = rho_B + theta_s[i+1]*rho_l
		rho_r[i+1] = rho_B + theta_r[i+1]*rho_l
	
		c_s[i+1] = c_B + theta_s[i+1]*rho_l*(c_l - c_B)/rho_s[i+1]
		c_r[i+1] = c_B + theta_r[i+1]*rho_l*(c_l - c_B)/rho_s[i+1]

		i = i + 1

	return(T_s)

# Neither of these functions are currently used because I got rid of moisture exchange between two layers except for when the top bucket overfills

#def K_theta(theta_dex,max_K):
# Function to get K (hydraulic conductivity) (m/sec)  according to RAWLS 1982
#
#	if theta_dex <= 0:
#		K = 0
#	else:
#		
#		K = max_K*(m.log(theta_dex) + 1)
#
#	if K < 0:
#		K = 0
#	
#	return(K)
# Function to get moisture diffusivity (m^2/s) according to Libardi et al. 1982. 
# Currently not used because I got rid of moisture exchange between two layers except for when the top bucket overfills
#
#def D_theta(theta_dex):
#	# Function to get D_theta (moisture diffusivity) (m^2/s) according to Libardi et al. 1982
#	
#	D = 10**(theta_dex*3 - 8) 	
#	D = 0
#	return(D)
