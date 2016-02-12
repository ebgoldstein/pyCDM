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