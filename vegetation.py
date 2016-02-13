from CDMparams import *

def vegevol(Veg, time, timestep, shoreline, h, dh_dt, sensParam, Tvegs, xmin, timefrac):

    #grow species
    V_gen = 1/Tvegs;

    #ADD Line for rho competition > rho max, then rho max

    #Heaviside for Lveg
    #calculate Lveg in model domain; so moving boundary will be hard...
    xmin = Lveg/dx;
    NoVeg=np.full((1,xmin),0);
    GrowVeg=np.full((1,nx-xmin),1);
    shorefactorprofile=np.concatenate((NoVeg,GrowVeg));
    shorefactor=np.tile(shorefactorprofile,(ny,1))

    # make dhdt grid

    #growth
    dV = ((1 - Veg ) * V_gen * shorefactor)- (abs(dhdt) * sensParam);

    #cover fraction evolves (timefrac is the rescaled wind time)
    Veg= Veg + (timestep * dV * timefrac);

    #limiting conditions (NEED TO VECTORIZE/ logical statements
    if Veg>1:
        Veg=1
    if Veg<0:
        Veg=0

    return veget
