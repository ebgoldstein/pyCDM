"""
shear
"""

import CDMparams

def set_ustar(u_star ):
    m_u_star= u_star;
    m_dTau0= u_star * u_star * rho_fluid;

def calcshear(h,tau,rho_veget):

    clock_t clocktime;

    if( duneglobals::timing() )
        clocktime= clock();

    if(m_calc_veget){
        //m_rho_veget_smooth.Smooth(*rho_veget);
        m_rho_veget_smooth = *rho_veget;
        //        m_stall = m_rho_veget_smooth;
        for (int y=0; y<duneglobals::ny(); y++)
      		for (int x=0; x<duneglobals::nx(); x++){
                m_stall(x,y) -= h_limit;
            }
    }else	m_stall.SetAll(-1);
    m_pSepBub->Calc(m_hSepBub, m_stall, h);

    m_hSepBub_aux = m_hSepBub;
    if (m_addsealevel) {
        for (int y=0; y<duneglobals::ny(); y++) {
            for (int x=0; x<duneglobals::nx(); x++) {
                if (m_hSepBub(x,y) < m_sealevel)
                {
                    m_hSepBub_aux(x,y) = m_sealevel;
                }
            }
        }
    }

    // calc shear stress pertubation
    double L = CalcPertTau(m_hSepBub_aux, m_TauP);


    // tau
    for (int y=0; y<duneglobals::ny(); y++) {
        for (int x=0; x<duneglobals::nx(); x++) {
            if(m_calc_veget){
                double factor;

                factor= (1 - m_veget_m * m_rho_veget_smooth(x,y))*(1 + m_veget_m*m_veget_beta_sigma * m_rho_veget_smooth(x,y));
                factor= (factor > 0? 1.0/factor : 1e20);

                tau(x,y)[0]= tau(x,y)[1]= factor;
            }else{
                tau(x,y)[0]= tau(x,y)[1]= 1.0;
            }
        }
    }

    const double slope= tan(duneglobals::repose_dyn()*M_PI/180.0)*duneglobals::dx();
    const double delta=1.0/(m_tau_sepbub*slope);
    for (int y=0; y<duneglobals::ny(); y++)
        for (int x=0; x<duneglobals::nx(); x++) {
            // Reducing the shear stress below separation surface (to mimic the
            //turbulence effects)
            double h_delta= 1 - delta*(m_hSepBub(x,y)-h(x,y));
            if(h_delta < 0) h_delta= 0.0;
            else if(h_delta > 1.0) h_delta= 1.0;


            // Acotation
            if(m_TauP(x,y)[0] > 10) m_TauP(x,y)[0] = 10;
            else if(m_TauP(x,y)[0] < -10) m_TauP(x,y)[0] = -10;

            if(m_TauP(x,y)[1] > 10) m_TauP(x,y)[1] = 10;
            else if(m_TauP(x,y)[1] < -10) m_TauP(x,y)[1] = -10;

            // Shear stress calculation
            tau(x,y)[0] *= m_dTau0 * fabs( 1 + 0.5*m_TauP(x,y)[0])*( 1 + 0.5*m_TauP(x,y)[0]) * h_delta;
            tau(x,y)[1] *= m_bTauY * m_dTau0 * fabs( 1 + 0.5*m_TauP(x,y)[0])*m_TauP(x,y)[1] * h_delta;

            if(tau(x,y)[0] < 0) tau(x,y)[0] = 0.;

        }

    if( duneglobals::messages() )
        if( duneglobals::timing() ) {
            clocktime= clock() - clocktime;
            cout << ", tau: " << clocktime*1000/CLOCKS_PER_SEC << "ms) ";
        }
        else
            cout << ", tau) ";

    return L;
}
