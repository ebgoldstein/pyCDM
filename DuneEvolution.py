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

m_ustar0=m_ustar;

def step_implementation()
	if m_ustar0 > 0.6*u_star_ft):
		if m_ustar0 > u_star_ft():
			m_Satflux_upwind = (m_ustar0 > u_star_ft() ? m_calcflux->Satflux_upwind(m_ustar0) : 0.0);
        
        m_calcshear->set_ustar(m_ustar0);
        if(m_calc_veget) halfmeanLength = m_calcshear->Calc( m_h, m_tau, m_rho_veget );
        else halfmeanLength = m_calcshear->Calc( m_h, m_tau);
        
        m_influx->set(m_flux_in, m_flux, m_Satflux_upwind, angle );
        m_gamma.SetAll(0.0);
//       m_calcflux->waterlevel_factor(evolution::time()); // get temporal variation of wet-dry level
        m_calcflux->calc( m_flux_in, m_flux, m_h, m_h_nonerod, m_tau, m_gamma );
        
        m_hprev = m_h; // copy previous profile;
        timestep= update_height(halfmeanLength);
        m_avalanche->calc(m_h, m_h_nonerod);
        
        update_dhdt();
        
    }else{
