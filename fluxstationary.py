
import CDMparams

def u_star_at(ustar)
	u_t = m_M * m_u_star_t - (m_M - 1.) * u_star;
    return(u_t)

def Satflux_upwind(u_star)
    u_s = 2.5 *u_star *(m_log_z1z0 - (1.- u_star_at(u_star)/u_star)*m_b_z1) - m_DeltaU;
    flux = m_2alpha_g * m_tau_t*((u_star/m_u_star_t)*(u_star/m_u_star_t) - 1) * u_s;
    return(flux)
    
