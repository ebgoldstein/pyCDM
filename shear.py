"""
shear
"""

import CDMparams

def set_ustar(u_star ):
    m_u_star= u_star;
    m_dTau0= u_star * u_star * rho_fluid;

def calcshear(h,tau,rho_veget):

    clock_t clocktime;

    if duneglobals::timing():
        clocktime= clock();

    m_rho_veget_smooth = *rho_veget;

    m_stall = m_stall-h_limit;

    m_pSepBub->Calc(m_hSepBub, m_stall, h);

    m_hSepBub_aux = m_hSepBub;
    if (m_addsealevel):
        if (m_hSepBub(x,y) < m_sealevel):
            m_hSepBub_aux(x,y) = m_sealevel;

    #calc shear stress pertubation
    L = CalcPertTau(m_hSepBub_aux, m_TauP);


    #tau (matrix)
    factor= (1 - m_veget_m * m_rho_veget_smooth(x,y))*(1 + m_veget_m*m_veget_beta_sigma * m_rho_veget_smooth(x,y));
    if factor > 0:
        factor= 1.0/factor
    else:
        factor= 1e20);

    tauZero= factor;
    tauOne=factor;

    slope= tan(duneglobals::repose_dyn()*M_PI/180.0)*duneglobals::dx();
    delta=1.0/(m_tau_sepbub*slope);
    for (int y=0; y<duneglobals::ny(); y++)
        for (int x=0; x<duneglobals::nx(); x++) {
            // Reducing the shear stress below separation surface (to mimic the
            //turbulence effects)
            h_delta= 1 - delta*(m_hSepBub(x,y)-h(x,y));
            if(h_delta < 0):
                 h_delta= 0.0;
            elif (h_delta > 1.0):
                 h_delta= 1.0;


            // Acotation
            if(m_TauP(x,y)[0] > 10):
                 m_TauP(x,y)[0] = 10;
            elif(m_TauP(x,y)[0] < -10):
                 m_TauP(x,y)[0] = -10;

            if(m_TauP(x,y)[1] > 10):
                 m_TauP(x,y)[1] = 10;
            elif(m_TauP(x,y)[1] < -10):
                 m_TauP(x,y)[1] = -10;

            // Shear stress calculation
            tauZero *= m_dTau0 * fabs( 1 + 0.5*m_TauP(x,y)[0])*( 1 + 0.5*m_TauP(x,y)[0]) * h_delta;
            tauOne *= m_bTauY * m_dTau0 * fabs( 1 + 0.5*m_TauP(x,y)[0])*m_TauP(x,y)[1] * h_delta;

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
