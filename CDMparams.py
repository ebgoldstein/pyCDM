"""
Some relevant CDM parameters and constants. This can act as the parameter file.

rewrite of CDM by EBG 2/16

"""
import numpy as np


#time scales
startstep= 0;
starttime= 0.0;
dt_max = 1000	# time step (constant, max is not important) (sec)
Nt = 100000 	# number of iterations (such that total time is Nt*dt_max)
rtime = 1; # fraction of the year wind is above threshold (only used for time scales) such that real time is Nt*dt_max / wind.frac
saveevery  = 100 	# nuber of iterations between saved files (such that the number of saved files are Nt/save.every)

#space scales
dx= 1;  			# (m)
nx= 100; 			# number of grid columns
ny= 4; 			# number of grid lines (c1: must be power of 2; c2: 4 is for 2D, higher for 3D)
length= nx*dx;
width= ny*dx;
threed= ny > 1;

#wind
constwind = 0.2828 # shear velocity (m/s) (usually from 0.2 [~transport threshold] to 0.5 [very large winds]); 0.2828 is ~2*transport threshold

# VEgetation parameters
Lveg=15; #distance for vegetation growth from shore (m)
Tveg=1; #characteristic growth timescale (in days)
Zmin=0; # threshold elevation for veg. growth: Z_veg (m) (DM15)
sens=1; #sensitivity to burial
Hveg=1; #characteristic vegetation height
rho_max=1; #max vegetation density

## wind reduction due to vegetation (from Duran et al. 2008, geomorphology):
vegetsigma = 0.80 	# ratio of plant basal to frontal area
vegetbeta 	= 150	# ratio of plant drag coefficitent to bare sand drag
vegetm 	= 0.16	# reduction parameter (empirical fitting parameter)
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

#Parameters for flux calculation*/
d_grain=  250e-6; #grain size in meters
fluid_viscosity=  1.8e-5;
fluid_viscosity_kin = fluid_viscosity / rho_fluid;	#kinematic viscosity
g=  9.8; #gravity

s = rho_grains / rho_fluid;
Aft = 0.11; #Bagnold-shields_parameter
At = 0.8 * Aft;
u_star_ft = Aft * np.sqrt(g*d_grain*(s-1.));
u_star_t = (At/Aft) * u_star_ft;

r= 0; #salt.Z/zm
#z0eff= (1e-3/*1e-3*/);

#Auxiliar quantities */
t_v = np.exp(1./3. * np.log(fluid_viscosity_kin / (g * g)));	#characteristic time
l_v = np.exp(1./3. * np.log(fluid_viscosity_kin * fluid_viscosity_kin / (At*At * g*(s-1)))); #characteristic lenght
S = 0.25 * d_grain / fluid_viscosity_kin * np.sqrt(g*d_grain*(s-1.));

#Parameters
C_diff=  0.0 ; #salt.D
beta=  5.7e-4 ; #beta
gamma=  0.1;	#gamma

#derived quantities:
z0= d_grain / 20.;
zm= 14./(1.+1.4*r) * u_star_t * t_v;
z1= 35. * l_v;

alpha= 0.18 * (d_grain / l_v);


Ad = 0.95;
Bd = 5.12; #natural sand
C_drag = 4./3. * (Ad + np.sqrt(2 * alpha) * Bd / S)*(Ad + np.sqrt(2 * alpha) * Bd / S);

DeltaU = np.sqrt(4.0/3.0 * g * (s-1) * d_grain / (C_drag * 2 * alpha));

z0zm= z0/zm;
z1zm= z1/zm;
log_z1z0= np.log(z1/z0);

if( r < z0zm ):
    r = z0zm;

log_r = np.log(r/z0zm);
b_inf = -0.5772 - np.log(z0zm) + z0zm;

if r<2:
    b_Z = log_r-(r - z0zm)+0.25*(r*r - z0zm*z0zm)-0.042*(r*r*r - z0zm*z0zm*z0zm)
else:
    b_Z =b_inf - np.exp((-r)/r);
if log_r>0:
    M= log_r/b_Z
else:
    M= 1;

if z1zm < 2:
    b_z1 =  log_z1z0-(z1zm - z0zm)+0.25*(z1zm*z1zm - z0zm*z0zm)-0.042*(z1zm*z1zm*z1zm - z0zm*z0zm*z0zm)
else:
    b_z1 = b_inf - np.exp(-z1zm)/z1zm;

m_Cd_dx= C_diff / (dx*dx);
tau_ft= rho_fluid * u_star_ft * u_star_ft;
tau_t= rho_fluid * u_star_t * u_star_t;
twoalpha_g= 2.0*alpha / g;
gamma_2alpha_g2= gamma / (twoalpha_g * twoalpha_g);
beta_2alpha_g_tau_ft = beta / (twoalpha_g * tau_ft);


def TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny):
    #how long is foreshore
    ForeshoreDistance=(shore_HMWL- shore_watertable)/(np.tan(np.deg2rad(beach_angle)));
    ForeshoreCells=ForeshoreDistance/dx;
    #inclined foreshore (DM13)
    Foreshore=np.arange(0,(shore_HMWL- shore_watertable),((shore_HMWL- shore_watertable)/dx));
    #flat backshore (DM13)
    BackshoreLength=nx-len(Foreshore);
    Backshore=np.array([1,BackshoreLength]);
    Backshore.fill(shore_HMWL- shore_watertable);
    #merge the foreshore to the backshore
    CrossShoreProfile=np.concatenate((Foreshore,Backshore))
    #tile it up to make the 2D domain
    Topo=np.tile(CrossShoreProfile,(ny,1))
    return Topo

def VegDomain(nx,ny):
    Veg=np.zeros((ny,nx))
    return Veg
