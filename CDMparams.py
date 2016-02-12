"""
Some relevant CDM parameters and constants. This can act as the parameter file. 

rewrite of CDM by EBG 2/16

"""


m_startstep= 0;
m_starttime= parameters.getdefault("time0", 0.0);
m_dx= parameters.getrequired<double>("dx");
m_nx= parameters.getrequired<int>( "NX" );
m_length= m_nx*m_dx;
m_ny= parameters.getrequired<int>( "NY" );
m_width= m_ny*m_dx;
m_3d= m_ny > 1;

m_z0eff= parameters.getdefault("hlr.z0", 1e-3/*1e-3*/);
    
m_repose_stat= parameters.getdefault("aval.angle_repose_stat",34);
m_repose_dyn= parameters.getdefault("aval.angle_repose_dyn",33);
    
m_periodic_x= parameters.getdefault( "calc.x_periodic", false );

# Periodic boundary condition in x direction requires NX to be power of 2!  
# Periodic boundary condition in y direction requires NY to be power of 2!  
    

    
m_rho_fluid= parameters.getdefault( "fluid_density", 1.225);
m_rho_grains= parameters.getdefault( "grain_density", 2650);
m_rho_sand= m_rho_grains * parameters.getdefault( "packing", 0.6226); //1650
    
    # Shore parameters*/
m_shore_HMWL = parameters.getdefault("shore.MHWL", 0.0);
m_shore_watertable = parameters.getdefault("shore.sealevel", 0.0);
m_beach_angle = parameters.getdefault("beach.angle", 0.0);
m_beach_slope = tan(m_beach_angle * M_PI / 180.);

# fraction of the year wind is above threshold (only used for time scales) such that real time is Nt*dt_max / wind.frac
m_rtime = 1;

#Units conversion*/
m_secday = 60*60*24; // seconds in a day
m_secmonth = 60*60*24*30; // seconds in a month (30 days)
m_secyear = 60*60*24*365; // seconds in a year (365 days)
    
#Parameters for flux calculation*/
m_d_grain= parameters.getdefault( "salt.d", 250e-6);
m_fluid_viscosity= parameters.getdefault( "fluid_viscosity", 1.8e-5);
fluid_viscosity_kin = m_fluid_viscosity / m_rho_fluid;
m_g= parameters.getdefault( "salt.g", 9.8 );
    
s = m_rho_grains / m_rho_fluid;
Aft = parameters.getdefault( "Bagnold-shields_parameter", 0.11);
At = 0.8 * Aft;
m_u_star_ft = Aft * sqrt(m_g*m_d_grain*(s-1.));
m_u_star_t = (At/Aft) * m_u_star_ft;

        
r= parameters.getdefault( "salt.Z/zm", 0.);
    
#Auxiliar quantities */
t_v = exp(1./3. * log(fluid_viscosity_kin / (m_g * m_g)));	//characteristic time
l_v = exp(1./3. * log(fluid_viscosity_kin * fluid_viscosity_kin / (At*At * m_g*(s-1))));	//characteristic lenght
S = 0.25 * m_d_grain / fluid_viscosity_kin * sqrt(m_g*m_d_grain*(s-1.));
    
#Parameters
m_C_diff= parameters.getdefault( "salt.D", 0.0 );
m_beta= parameters.getdefault( "beta", 5.7e-4 );
m_gamma= parameters.getdefault( "gamma", 0.1);
        
#derived quantities:
m_z0= m_d_grain / 20.;
m_zm= 14./(1.+1.4*r) * m_u_star_t * t_v;
m_z1= 35. * l_v;

m_alpha= 0.18 * (m_d_grain / l_v);

    
Ad = 0.95, Bd = 5.12; //natural sand
m_C_drag = 4./3. * (Ad + sqrt(2 * m_alpha) * Bd / S)*(Ad + sqrt(2 * m_alpha) * Bd / S);

    
m_DeltaU = sqrt(4.0/3.0 * m_g * (s-1) * m_d_grain / (m_C_drag * 2 * m_alpha));
    
z0zm= m_z0/m_zm;
z1zm= m_z1/m_zm;
m_log_z1z0= log(m_z1/m_z0);
    
    if( r < z0zm ):
    	r = z0zm;
log_r = log(r/z0zm);    
b_inf = -0.5772 - log(z0zm) + z0zm;
b_Z = (r < 2. ? log_r-(r - z0zm)+0.25*(r*r - z0zm*z0zm)-0.042*(r*r*r - z0zm*z0zm*z0zm) : b_inf - exp(-r)/r);
m_M = (log_r > 0 ? log_r/b_Z : 1.);

    
m_b_z1 = (z1zm < 2. ? m_log_z1z0-(z1zm - z0zm)+0.25*(z1zm*z1zm - z0zm*z0zm)-0.042*(z1zm*z1zm*z1zm - z0zm*z0zm*z0zm) : b_inf - exp(-z1zm)/z1zm);
m_Cd_dx= m_C_diff / (m_dx*m_dx);
m_tau_ft= m_rho_fluid * m_u_star_ft * m_u_star_ft;
m_tau_t= m_rho_fluid * m_u_star_t * m_u_star_t;
m_2alpha_g= 2.0*m_alpha / m_g;
m_gamma_2alpha_g2= m_gamma / (m_2alpha_g * m_2alpha_g);
m_beta_2alpha_g_tau_ft = m_beta / (m_2alpha_g * m_tau_ft);
    
m_datadir= parameters.getdefault( "save.dir", string("DAT") );


