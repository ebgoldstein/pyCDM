import numpy as np
from CDMparams import *

def TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny):
    #how long is foreshore (round up to nearest meter)

    ForeshoreDistance=(shore_HMWL- shore_watertable)/(np.tan(np.deg2rad(beach_angle)));
    ForeshoreCells=np.round(ForeshoreDistance/dx,0);
    #inclined foreshore (DM13)
    Foreshore=np.linspace(0,(shore_HMWL- shore_watertable),ForeshoreCells)
    #flat backshore (DM13)
    BackshoreLength=nx-len(Foreshore);
    Backshore=np.full((BackshoreLength), (shore_HMWL- shore_watertable))
    #merge the foreshore to the backshore
    CrossShoreProfile=np.concatenate((Foreshore,Backshore))
    #tile it up to make the 2D domain
    Topo=np.tile(CrossShoreProfile,(ny,1))
    return Topo

def VegDomain(nx,ny):
    Veg=np.zeros((ny,nx))
    return Veg

def TauDomain(nx,ny):
    Tau=np.zeros((ny,nx))
    return Tau

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
    Veg= Veg + (timestep * dV * timefrac);

    #limiting conditions (can't have cover density greater than rhomax or less than 0))
    Veg[Veg > 1] = 1;
    Veg[Veg < 0] = 0;

    return Veg

def u_star_at(ustar):
	u_t = (M * u_star_t) - ((M - 1.) * u_star);
    return u_t

def Satflux_upwind(u_star):
    u_s = 2.5 *u_star *(log_z1z0 - (1.- (u_star_at(u_star)/u_star)*b_z1)) - DeltaU;
    flux = twoalpha_g * tau_t * u_s * (((u_star/u_star_t)*(u_star/u_star_t)) - 1);
    return flux

def set_ustar(u_star ):
    u_star= u_star;
    dTau0= u_star * u_star * rho_fluid;
    return dTau0

def calcshear(Topo,Veg):

    stall = stall-h_limit;

    pSepBub->Calc(hSepBub, stall, Topo);

    hSepBub_aux = hSepBub;
    if (addsealevel):
        if (hSepBub(x,y) < sealevel):
            hSepBub_aux(x,y) = sealevel;

    #calc shear stress pertubation
    L = CalcPertTau(hSepBub_aux, TauP); 


    #tau (matrix)
    factor= (1 - (vegetm * Veg))*(1 + (vegetm*(vegetbeta/vegetsigma) * Veg));
    if factor > 0:
        factor= 1.0/factor
    else:
        factor= 1e20);

    tauZero= factor;
    tauOne=factor;

    slope= tan(duneglobals::repose_dyn()*PI/180.0)*duneglobals::dx();
    delta=1.0/(tau_sepbub*slope);

    #Reducing the shear stress below separation surface (to mimic the
    #turbulence effects)
    h_delta= 1 - delta*(hSepBub(x,y)-h(x,y));
    if(h_delta < 0):
        h_delta= 0.0;
    elif (h_delta > 1.0):
        h_delta= 1.0;


    # Acotation
    if(TauP(x,y)[0] > 10):
            TauP(x,y)[0] = 10;
    elif(TauP(x,y)[0] < -10):
        TauP(x,y)[0] = -10;

    if(TauP(x,y)[1] > 10):
        TauP(x,y)[1] = 10;
    elif(TauP(x,y)[1] < -10):
        TauP(x,y)[1] = -10;

    # Shear stress calculation
    tauZero *= dTau0 * fabs( 1 + 0.5*TauP(x,y)[0])*( 1 + 0.5*TauP(x,y)[0]) * h_delta;
    tauOne *= bTauY * dTau0 * fabs( 1 + 0.5*TauP(x,y)[0])*TauP(x,y)[1] * h_delta;

    if(tauZero < 0):
        tauZero = 0;


    if duneglobals::messages():
        if( duneglobals::timing():
            clocktime= clock() - clocktime;
            cout << ", tau: " << clocktime*1000/CLOCKS_PER_SEC << "ms) ";
        }
        else:
            cout << ", tau) ";

    return L
