"""
rewrite of CDM by EBG 2/16

The CDM operates by iterating over 5 operations with 2 states:

The States are:
 	-'Topo' (the topographic domain)
	-'Veg' (the vegetation field)

The 5 operations are:
	1) Compute the shear stress field. The shear stress field is computed by
		using the u* defined for each time iteration and the surface topography.
		The surface topography sets the shear stress through a linearized RANS
		formulation (Weng et al. 1991). This is an analytical solution utilizing
		a fourier transform. The topographic surface is
		modified by the existence of a seperation bubble when there is a dune
		brink. The resultant shear stress, from the linearized RANS formulation,
		is then modified by the vegetation field
	2) A sand flux field is computed from the shear stress field. Flux in the
		separation zone is set to 0.
	3) The flux divergence is computed for a given time interval (forward Euler)
		and this sets the change in topography
	4) The new topo field is polled for slopes greater than the angle of repose
		(34 degrees). For slopes greater than this, avalanches occur down the
		steepest descent gradient
	5) The vegetation field is modified based on growth equations, which
		encompass the accretion/erosion rate of the sand surface.

"""

import numpy as np
from CDMparams import *
import CDMfunctions as cdmfns

"""
Build the initial grids
"""
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


"""
Step 1: shear stress field
"""

def shearfield(Topo, Veg, Tau):

    #find the local slope
    localslope=np.diff(Topo)

    #if (slope of 20 deg or greater it's a brink; Duran and Moore, 2013)
        #separation bubble; 3rd order polynomial
        #make new combined surface
        Topowind=np.maximum(SepBubble,Topo)
    #else
        Topowind =Topo

    #calculate tau perterbation

    #calculate combined tau

    #calculate reduction factor for Vegetation

    return Tau


"""
Step 2: sand flux  field
"""
def sandfluxfield(Topo, Tau):

return Flux

"""
Step 3: flux divergence and building topo
"""
def fluxgradient(Topo,Flux):

    return Topo, dhdt

"""
Step 4: avalanche
"""
def avalanche(Topo,dhdt):

    return Topo, dhdt

"""
Step 5: shear stress field
"""

def vegevol(Veg, Topo, dhdt):
    #grow species
    V_gen = 1/Tvegs;

    #ADD Line for rho competition > rho max, then rho max

    #Heaviside for Lveg
    #calculate Lveg in model domain; probably won;t work with a moving boundary yet
    xmin = Lveg/dx;
    NoVeg=np.full((xmin),0);
    GrowVeg=np.full((nx-xmin),1);
    shorefactorprofile=np.concatenate((NoVeg,GrowVeg));
    shorefactor=np.tile(shorefactorprofile,(ny,1))

    #growth
    dV = ((1 - Veg ) * V_gen * shorefactor)- (abs(dhdt) * sensParam);

    #cover fraction evolves (timefrac is the rescaled wind time)
    Veg= Veg + (timestep * dV * timefrac);

    #limiting conditions (can't have cover density greater than rhomax or less than 0))
    Veg[Veg > 1] = 1;
    Veg[Veg < 0] = 0;

    return Veg

"""
Implement all 5 steps in a row, with a time step..
This is the main script, which runs the model
"""

def runCDM:

    #MAKE THE INITIAL CONDIITONS!
    # (should make these 3D for time)
    #sand domain
    Topo=TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny)
    #vegetation domain
    Veg=VegDomain(nx,ny)
    #shear stress domain
    Tau=TauDomain(nx,ny)


    #for initial time: timestep: final time
	   #1. Compute the stress field

       #2. Compute the flux field

       #3. Compute flux divergence and update topo

       #4. Poll the grid for avalanches

       #5. Compute the vegetation field
       Veg = vegevol(Veg, Topo, dhdt);

       #save topo and veg each increment

    return Topo, Veg
