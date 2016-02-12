"""
Some relevant CDM parameters and constants. This can act as the parameter file.

rewrite of CDM by EBG 2/16

"""


m_startstep= 0;
m_starttime= 0.0;
m_dx= 1;  			# (m)
m_nx= 100; 			# number of grid columns
m_ny= 4; 			# number of grid lines (c1: must be power of 2; c2: 4 is for 2D, higher for 3D)

m_length= m_nx*m_dx;
m_width= m_ny*m_dx;
m_3d= m_ny > 1;

m_z0eff= (1e-3/*1e-3*/);

m_repose_stat= 34; 	#static angle of repose
m_repose_dyn= 33; #dynamic angle of repose

m_periodic_x = false;
m_periodic_y = false;
# Periodic boundary condition in x direction requires NX to be power of 2!
# Periodic boundary condition in y direction requires NY to be power of 2!

# VEgetation parameters
m_Lveg=15; #distance for vegetation growth from shore (m)
m_Tveg=1; #characteristic growth timescale (in days)
m_sens=1; #sensitivity to burial
m_Hveg=1; #characteristic vegetation height
m_rho_max=1; #max vegetation density

m_rho_fluid= 1.225; #density of air (kg.m^3)
m_rho_grains= 2650; #density of grains (kg.m^3)
m_rho_sand= m_rho_grains * 0.6226; # packing

# Shore parameters
m_shore_HMWL =  0.3; 		# MHWL relative to the watertable
m_shore_watertable =  0.0; #sea level
m_beach_angle =  1; #degrees
m_beach_slope = tan(m_beach_angle * 3.14 / 180.);

# fraction of the year wind is above threshold (only used for time scales) such that real time is Nt*dt_max / wind.frac
m_rtime = 1;

#Units conversion*/
m_secday = 60*60*24; # seconds in a day
m_secmonth = 60*60*24*30; # seconds in a month (30 days)
m_secyear = 60*60*24*365; # seconds in a year (365 days)

#Parameters for flux calculation*/
m_d_grain=  250e-6; #grain size in meters
m_fluid_viscosity=  1.8e-5;
fluid_viscosity_kin = m_fluid_viscosity / m_rho_fluid;	#kinematic viscosity
m_g=  9.8; #gravity

s = m_rho_grains / m_rho_fluid;
Aft = 0.11; #Bagnold-shields_parameter
At = 0.8 * Aft;
m_u_star_ft = Aft * sqrt(m_g*m_d_grain*(s-1.));
m_u_star_t = (At/Aft) * m_u_star_ft;

r= 0; #salt.Z/zm

#Auxiliar quantities */
t_v = exp(1./3. * log(fluid_viscosity_kin / (m_g * m_g)));	#characteristic time
l_v = exp(1./3. * log(fluid_viscosity_kin * fluid_viscosity_kin / (At*At * m_g*(s-1)))); #characteristic lenght
S = 0.25 * m_d_grain / fluid_viscosity_kin * sqrt(m_g*m_d_grain*(s-1.));

#Parameters
m_C_diff=  0.0 ; #salt.D
m_beta=  5.7e-4 ; #beta
m_gamma=  0.1;	#gamma

#derived quantities:
m_z0= m_d_grain / 20.;
m_zm= 14./(1.+1.4*r) * m_u_star_t * t_v;
m_z1= 35. * l_v;

m_alpha= 0.18 * (m_d_grain / l_v);


Ad = 0.95, Bd = 5.12; //natural sand
m_C_drag = 4./3. * (Ad + sqrt(2 * m_alpha) * Bd / S)*(Ad + sqrt(2 * m_alpha) * Bd / S);

m_DeltaU = sqrt(4.0/3.0 * m_g * (s-1) * m_d_grain / (m_C_drag * 2 * m_alpha));

z0zm= m_z0/m_zm;
z1zm= m_z1/m_zm;
m_log_z1z0= log(m_z1/m_z0);

if( r < z0zm ):
    r = z0zm;

log_r = log(r/z0zm);
b_inf = -0.5772 - log(z0zm) + z0zm;

if r<2:
    b_Z = log_r-(r - z0zm)+0.25*(r*r - z0zm*z0zm)-0.042*(r*r*r - z0zm*z0zm*z0zm)
else:
    b_Z =b_inf - exp(-r)/r);
if log_r>0:
    m_M= log_r/b_Z
else:
    m_M= 1;

if z1zm <
    m_b_z1 =  m_log_z1z0-(z1zm - z0zm)+0.25*(z1zm*z1zm - z0zm*z0zm)-0.042*(z1zm*z1zm*z1zm - z0zm*z0zm*z0zm)
else:
    m_b_z1 = b_inf - exp(-z1zm)/z1zm);

m_Cd_dx= m_C_diff / (m_dx*m_dx);
m_tau_ft= m_rho_fluid * m_u_star_ft * m_u_star_ft;
m_tau_t= m_rho_fluid * m_u_star_t * m_u_star_t;
m_2alpha_g= 2.0*m_alpha / m_g;
m_gamma_2alpha_g2= m_gamma / (m_2alpha_g * m_2alpha_g);
m_beta_2alpha_g_tau_ft = m_beta / (m_2alpha_g * m_tau_ft);
m_datadir= DAT;	#saving directory
