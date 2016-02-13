from CDMparams import *


#Vegetation parameters*/
m_xmin = m_Lveg/m_dx; #calculate Lveg in model domain
m_Tvegs = m_Tveg* secday #Tveg in seconds
m_sensParam=m_sens/m_Hveg #sensitivityParameter from DM13

#time conversion
m_wind_factor =  timefrac;


def vegevol(rho_veget, time, timestep, shoreline, h, dh_dt):

    #grow species
    V_gen = 1/m_Tvegs;


    #ADD Line for rho competition > rho max, then rho max

    #Heaviside for Lveg
    shorefactor = (x < shoreline + m_xmin ? 0 : 1);

    #growth
    dV = ((1 - rho_veget ) * V_gen * shorefactor)- (abs(dhdt) * m_sens;);

    #cover fraction evolves
    m_veget= m_veget + (timestep * dV * m_wind_factor);

    #limiting conditions
    if m_veget>1:
        m_veget=1
    if m_veget<0:
        m_veget=0

    return m_veget
