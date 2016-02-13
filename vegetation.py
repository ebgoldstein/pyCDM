from CDMparams import *

#def vegevol(Veg, time, timestep, shoreline, h, dh_dt, sensParam, Tvegs, xmin, timefrac):
def vegevol(Veg, Topo, dhdt):
    #grow species
    V_gen = 1/Tvegs;

    #ADD Line for rho competition > rho max, then rho max

    #Heaviside for Lveg
    #calculate Lveg in model domain; so moving boundary will be hard...
    xmin = Lveg/dx;
    NoVeg=np.full((xmin),0);
    GrowVeg=np.full((nx-xmin),1);
    shorefactorprofile=np.concatenate((NoVeg,GrowVeg));
    shorefactor=np.tile(shorefactorprofile,(ny,1))

    # make dhdt grid

    #growth
    dV = ((1 - Veg ) * V_gen * shorefactor)- (abs(dhdt) * sensParam);

    #cover fraction evolves (timefrac is the rescaled wind time)
    Veg= Veg + (dt_max * dV * timefrac);

    #limiting conditions (can't have cover density greater than rhomax or less than 0))
    Veg[Veg > 1] = 1;
    Veg[Veg < 0] = 0;

    return Veg
