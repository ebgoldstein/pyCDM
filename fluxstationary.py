
from CDMparams import *

def u_star_at(ustar):
	u_t = (M * u_star_t) - ((M - 1.) * u_star);
    return u_t

def Satflux_upwind(u_star):
    u_s = 2.5 *u_star *(log_z1z0 - (1.- (u_star_at(u_star)/u_star)*b_z1)) - DeltaU;
    flux = twoalpha_g * tau_t * u_s * (((u_star/u_star_t)*(u_star/u_star_t)) - 1);
    return flux
