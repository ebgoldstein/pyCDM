"""
Some relevant CDM parameters and constants

rewrite of CDM by EBG 2/16

"""
#time constants
m_secday = 60*60*24; #seconds in a day
m_secmonth = 60*60*24*30; #seconds in a month (30 days)
m_secyear = 60*60*24*365; #seconds in a year (365 days)



#Flux constants
m_g= 9.8 #gravity
m_d_grain = 250e-6; #grain size (m) (250e-6 is fine sand)
m_rho_fluid = 1.225; #fluid density (kg/m^3)
m_rho_grains = 2650; #grain density (kg/m^3)
s = m_rho_grains / m_rho_fluid;
Aft = 0.11 #Bagnold-shields_parameter
At = 0.8 * Aft;
m_u_star_ft = Aft * sqrt(m_g*m_d_grain*(s-1.));
m_u_star_t = (At/Aft) * m_u_star_ft;
r=0; #salt.Z/zm from CDM


m_b_z1 = duneglobals::b_z1();

m_2alpha_g = duneglobals::alpha2_g();
m_tau_t = duneglobals::tau_t();
m_u_star_t = duneglobals::u_star_t();
m_M = duneglobals::M();
m_C_drag
fluid_viscosity_kin
t_v



l_v = exp(1./3. * log(fluid_viscosity_kin * fluid_viscosity_kin / (At*At * m_g*(s-1))));	//characteristic lenght
m_z0= m_d_grain / 20.;
m_zm= 14./(1.+1.4*r) * m_u_star_t * t_v;
m_z1= 35. * l_v;
z0zm= m_z0/m_zm;
z1zm= m_z1/m_zm;
m_log_z1z0= log(m_z1/m_z0);

if( r < z0zm ):
	r = z0zm;

m_alpha= 0.18 * (m_d_grain / l_v);
m_DeltaU = sqrt(4.0/3.0 * m_g * (s-1) * m_d_grain / (m_C_drag * 2 * m_alpha));




