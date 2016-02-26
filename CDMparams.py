"""
Some relevant CDM parameters and constants. This can act as the parameter file.

rewrite of CDM by EBG 2/16

"""
import numpy as np
import CDMfunctions as cdmfns


#time scales
startstep= 0;
starttime= 0.0;
dt_max = 1000	# time step (constant, max is not important) (sec)
Nt = 100000 	# number of iterations (such that total time is Nt*dt_max)
timefrac = 1; # fraction of the year wind is above threshold (only used for time scales) such that real time is Nt*dt_max / wind.frac
saveevery  = 100 	# nuber of iterations between saved files (such that the number of saved files are Nt/save.every)
timestep=dt_max;

#space scales
dx= 1;  			# (m)
nx= 100; 			# number of grid columns
ny= 1; 			# number of grid lines (c1: must be power of 2; c2: 4 is for 2D, higher for 3D)
length= nx*dx;
width= ny*dx;
threed= ny > 1;

#wind (5 types of wind total:constant, real, flatrand, sine, bidirectional )
constwind = 0.2828 # shear velocity (m/s) (usually from 0.2 [~transport threshold] to 0.5 [very large winds]); 0.2828 is ~2*transport threshold
ustar = constwind;
wdir = 0; # onshore wind


# VEgetation parameters
Lveg=30; #distance for vegetation growth from shore (m)
Tveg=1; #characteristic growth timescale (in days)
Zmin=0; # threshold elevation for veg. growth: Z_veg (m) (DM15)
sens=1; #sensitivity to burial
Hveg=1; #characteristic vegetation height
rho_max=1; #max vegetation density
sensParam=sens/Hveg #sensitivityParameter from DM13
xmin = Lveg/dx; #calculate Lveg in model domain

## wind reduction due to vegetation (from Duran et al. 2008, geomorphology):
vegetsigma = 0.80 	# ratio of plant basal to frontal area
vegetbeta 	= 150	# ratio of plant drag coefficitent to bare sand drag
vegetm 	= 0.16	# reduction parameter (empirical fitting parameter)
#combined Roughness from Duran & Moore (2013)
vegetgamma= (vegetm*vegetbeta)/vegetsigma
# Roughness factor from Duran & Moore (2013) =veget.m * veget.beta / veget.sigma

rho_fluid= 1.225; #density of air (kg.m^3)
rho_grains= 2650; #density of grains (kg.m^3)
rho_sand= rho_grains * 0.6226; # packing

# Shore parameters
shore_HMWL =  0.3; 		# MHWL relative to the watertable
shore_watertable =  0.0; #sea level
beach_angle =  2; #degrees
beach_tau_t_L = 0.05 	# ( h_w from Duran & Moore, 2013)

#saving
#save.dir    = ./Test # name of directory where output data is stored


repose_stat= 34; 	#static angle of repose
repose_dyn= 33; #dynamic angle of repose

#periodic_x = false;
#periodic_y = false;
# Periodic boundary condition in x direction requires NX to be power of 2!
# Periodic boundary condition in y direction requires NY to be power of 2!

#Units conversion*/
secday = 60*60*24; # seconds in a day
secmonth = 60*60*24*30; # seconds in a month (30 days)
secyear = 60*60*24*365; # seconds in a year (365 days)
Tvegs = Tveg * secday; #Tveg in seconds

#Parameters for flux calculation*/
d_grain=  250e-6; #grain size in meters
fluid_viscosity=  1.8e-5;
fluid_viscosity_kin = fluid_viscosity / rho_fluid;	#kinematic viscosity
g=  9.8; #gravity

s = rho_grains / rho_fluid;

alpha = 0.43 #parteli et al 2009
gamma = 0.2 #parteli et al 2009
u_star_t = 0.22 #parteli et al 2009
