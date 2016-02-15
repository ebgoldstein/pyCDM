"""
rewrite of CDM by EBG 2/16

This is the function which computes the actual change in the surface
profile.

The method for time evolution is forward Euler.
	1) The shear stress is computed with the help of m_calcshear.
	2) The the sand flux is calculated with m_calcflux.
	3) The rate of change of the height profile is calculated as the divergence of the
		flux divided by the sand density.
	4) m_h is changed according to this rate with a time step which is at most m_dtmax
		(parameter salt.dt_h_max in the .par file) and which is small enough that the
		change in height is at most m_dhmax (parameter salt.dh_max).
"""

import numpy as np
from CDMparams import *
import CDMfunctions as cdm


def step_implementation(Veg,Topo,Tau,ustar):
	ustar0=ustar;

	if ustar0 > 0.6*u_star_ft:
		if ustar0 > u_star_ft:
			Satflux_upwind =  cdm.Satflux_upwind(ustar0);
		else:
			Satflux_upwind = 0.0;

        dTau0=cdm.set_ustar(ustar0);
		halfmeanLength = cdm.calcshear( Topo, Tau, Veg );

        influx->set(flux_in, flux, Satflux_upwind, angle );
        gamma.SetAll(0.0);
        calcflux( flux_in, flux, Topo, h_nonerod, tau, gamma );

        Topoprev = Topo; #copy previous profile;
        timestep= update_height(halfmeanLength);

		calcavalanche(Topo, h_nonerod);

        dhdt = (Topo - Topoprev) / dtmax;

    #update vegetation
    Veg = cdm.vegevol(Veg, Topo, dhdt);


    #PROCESS DATA
	#steps = time() / timestep;
	#process = (steps % 50 == 0 ? 1 : 0);
	#analyzeCalc(steps, time(), m_qin/m_Satflux_upwind/duneglobals::ny(), m_qout/m_Satflux_upwind/duneglobals::ny(), m_h, m_rho_veget);
	return Veg,Topo




#Main pyCDM script. runs the model

#Make/get the initial conditions (should make these 3D for time)
#sand domain
Topo=cdm.TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny)
#vegetation domain
Veg=cdm.VegDomain(nx,ny)
#shear stress domain
Tau=cdm.TauDomain(nx,ny)


#for initial time: timestep: final time
  #step_implementation(Veg,Topo,ustar)
  #save each increment
