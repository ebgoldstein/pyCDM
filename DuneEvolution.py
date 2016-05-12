"""
rewrite of CDM by EBG 2/16
The MIT License (MIT)
Copyright (c) 2016 Evan B. Goldstein

The CDM operates by iterating over 5 operations with 2 states:

The States are:
 	-'Topo' (the topographic domain)
	-'Veg' (the vegetation field)

The 5 operations are:
	1) Compute the shear stress field. The shear stress field is computed by
		using the u* defined for each time iteration and the surface topography.
		The surface topography sets the shear stress through a analytical
		formulation (Weng et al. 1991). This is an analytical solution utilizing
		a fourier transform. The topographic surface is
		modified by the existence of a seperation bubble when there is a dune
		brink. The resultant shear stress is then modified by the vegetation field
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
    Tau=np.ones((ny,nx))
    return Tau


"""
Step 1: shear stress field

TO DO:
    seperation bubble is only in 2D
    No wind routine yet (weng )
"""

def shearfield(Topo, Veg, Tau):

    #find the local slope (offset by 1 in the x direction)
    localslope=100*np.diff(Topo)*(1/dx)

    # if anwhere in the domain there is a slope of -20 deg or less, it's a brink
    if np.amin(localslope)<-20:

        #make a new array highlighting the cell of the brink
        #(need to add a column at the end b/c slope domain is 1 less col. than actual domain)
        SepBubble=np.concatenate(localslope,np.zeros((ny,1)))
        SepBubble[SepBubble >= -20] = 0
        SepBubble[SepBubble < -20] = 1

        #build a separation bubble from brink to surface;
        #
        #3rd order polynomial (clamped spline)
        #MAX slope is 14 degrees or 0.244346 radians (Kroy et al 2002)
        #find index of the brink
        Brink=np.argwhere(SepBubble == 1);

        #note that this technique works for 1 Sep Bubble
        if Brink.size>1:
            print ('more than 1 brink!!')
        #return the topography after the Brink
        #note that this next line works for only a 1D array currently
        TunderSB=Topo[Brink[0]:Topo.size+1]
        xx=dx*np.arange(1,TunderSB.size)
        y=dx*Topo[Brink[0]]
        yy=dx*TunderSB[1:TunderSB.size+1]
        SBt=0.5
        SBa=-(yy-y)
        SBb=+(yy-y)
        #solve the derivative of the clamped spline equation to get max slope at midpoint
        #https://en.wikipedia.org/wiki/Spline_interpolation
        slope=np.arctan(((yy-y)/(xx)) + ((0.25)*((SBb-SBa)/(xx))) )
        #pick first value of mid point slope below 0.14
        reattachment=np.argwhere(slope > -0.245)

        #interpolate to make the seperation bubble
        SepBubbleL=dx*np.arange(1,reattachment[0])
        yy=Topo[Brink[0]+reattachment[0]]
        SBa=-(yy-y)
        SBb=+(yy-y)
        SBt=SepBubbleL/(reattachment[0])
        SepBubbleH=((1-SBt)*y)+ (SBt*yy) + (SBt*(1-SBt)*((SBa*(1-SBt))+(SBb*SBt)))

        #use that index, and add it to the Brink
        SepBubble[Brink[0]:reattachment[0]+1]=SepBubbleH

        #make new combined surface
        Topowind=np.maximum(SepBubble,Topo)

    else:
        Topowind =Topo

    #calculate tau perterbation using Weng et al 1991
    if ny=1: #Weng et al 1991 EQN 2.8
        H
        L
        U
        l
        sigmazero
        znaught
        zeta
        Knaught
        i
        k
        zetanaught

        Wtermone = -2*(H/L)/(U*U*l)
        Wparen=
        #calculate the perturbation
        TauPert= Wtermone*sigmazero*Wparen

    else:   #Weng et al 1991 EQN 2.14a,b


    #calculate total tau by adding the perturbation
    Ttau=Tau+((np.absolute(Tau))*TauPert)

    #all locations of seperation bubble sites have zero tau
    Stau=np.where((SepBubble==1),0,Ttau)

    #calculate reduction factor for vegetation
    #(could replace with Okin)
    Vtau=Tau/(1+(vegetgamma*Veg))
    Tau=np.where((Veg > 0), Vtau, Stau)

    return Tau


"""
Step 2: sand flux  field

to do:
    -mean grain velocity at saturation
    -calculate local (height integrated) flux
"""
def sandfluxfield(Topo, Tau):

    ustarfield= np.sqrt(Tau/rho_fluid)

    #mean grain velocity at saturation
    vs=

    #saturated length
    ls=((2*vs*vs*alpha)/(gamma*g)) / (np.power((ustarfield/u_star_t),2)-1)

    #saturated flux
    qs=((2*vs*alpha)/(g)) * rho_fluid * u_star_t * u_star_t * (np.power((ustarfield/u_star_t),2)-1)

    #calculate local height integrated flux
    Flux=

return Flux

"""
Step 3: flux divergence and building topo

to do:
    -check the signs
    -front of the beach should not erode (infinite source; Duran and Moore 2013 )

"""
def fluxgradient(Topo,Flux):

    dhdt = - (Flux / rho_sand)

    Topo = Topo + dhdt

    return Topo, dhdt

"""
Step 4: avalanche

to do:
    -All of it

"""
def avalanche(Topo,dhdt):

    return Topo, dhdt

"""
Step 5: Vegetation growth

    This is done
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

    To do:
        -Fill in the master loop and add a counter for time steps
"""

def runCDM:

    #MAKE THE INITIAL CONDIITONS!
    # (should make these 3D for time)
    #sand domain
    Topo=TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny)
    #vegetation domain
    Veg=VegDomain(nx,ny)

    #for initial time: timestep: final time

       #0. Set the wind for the timestep
       Tau=np.ones((ny,nx))*(ustar*ustar*rho_fluid)

       #1. Compute the stress field
       Tau=shearfield(Topo, Veg, Tau)

       #2. Compute the flux field
       Flux=sandfluxfield(Topo, Tau)

       #3. Compute flux divergence and update topo
       Topo, dhdt = fluxgradient(Topo,Flux)

       #4. Poll the grid for avalanches
       Topo,dhdt = avalanche(Topo,dhdt)

       #5. Compute the vegetation field
       Veg = vegevol(Veg, Topo, dhdt);

       #save topo and veg each increment

    return Topo, Veg
