from CDMparams import *


#Vegetation parameters*/
xmin = Lveg/dx; #calculate Lveg in model domain
Tvegs = Tveg* secday #Tveg in seconds
sensParam=sens/Hveg #sensitivityParameter from DM13

#time conversion
wind_factor =  timefrac;


def vegevol(rho_veget, time, timestep, shoreline, h, dh_dt):

    #grow species
    V_gen = 1/Tvegs;


    #ADD Line for rho competition > rho max, then rho max

    #Heaviside for Lveg
    shorefactor = (x < shoreline + xmin ? 0 : 1);

    #growth
    dV = ((1 - rho_veget ) * V_gen * shorefactor)- (abs(dhdt) * sens;);

    #cover fraction evolves
    veget= veget + (timestep * dV * wind_factor);

    #limiting conditions
    if veget>1:
        veget=1
    if veget<0:
        veget=0

    return veget
