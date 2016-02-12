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
import CDMparams
import wind
import shear
import fluxstationary
import vegetation


def step_implementation():
	m_ustar0=m_ustar;

	if m_ustar0 > 0.6*u_star_ft:
		if m_ustar0 > u_star_ft:
			if m_ustar0 > u_star_ft:
				m_Satflux_upwind =  Satflux_upwind(m_ustar0);
			else:
				m_Satflux_upwind = 0.0;

        set_ustar(m_ustar0);
		halfmeanLength = calcshear( m_h, m_tau, m_rho_veget );


        m_influx->set(m_flux_in, m_flux, m_Satflux_upwind, angle );
        m_gamma.SetAll(0.0);
        calcflux( m_flux_in, m_flux, m_h, m_h_nonerod, m_tau, m_gamma );

        m_hprev = m_h; #copy previous profile;
        timestep= update_height(halfmeanLength);

		calcavalanche(m_h, m_h_nonerod);

        m_dh_dt = (m_h - m_hprev) / m_dtmax;

    #update vegetation
    m_veget_X0 = vegevol(m_rho_veget, time(), timestep, m_shoreline, m_h, m_dh_dt);


    #PROCESS DATA
	#steps = time() / timestep;
	#process = (steps % 50 == 0 ? 1 : 0);
	#analyzeCalc(steps, time(), m_qin/m_Satflux_upwind/duneglobals::ny(), m_qout/m_Satflux_upwind/duneglobals::ny(), m_h, m_rho_veget);
