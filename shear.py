"""
shear
"""

from CDMparams import *

def set_ustar(u_star ):
    u_star= u_star;
    dTau0= u_star * u_star * rho_fluid;

def calcshear(h,tau,rho_veget):

    clock_t clocktime;

    if duneglobals::timing():
        clocktime= clock();

    rho_veget_smooth = *rho_veget;

    stall = stall-h_limit;

    pSepBub->Calc(hSepBub, stall, h);

    hSepBub_aux = hSepBub;
    if (addsealevel):
        if (hSepBub(x,y) < sealevel):
            hSepBub_aux(x,y) = sealevel;

    #calc shear stress pertubation
    L = CalcPertTau(hSepBub_aux, TauP);


    #tau (matrix)
    factor= (1 - veget_m * rho_veget_smooth(x,y))*(1 + veget_m*veget_beta_sigma * rho_veget_smooth(x,y));
    if factor > 0:
        factor= 1.0/factor
    else:
        factor= 1e20);

    tauZero= factor;
    tauOne=factor;

    slope= tan(duneglobals::repose_dyn()*PI/180.0)*duneglobals::dx();
    delta=1.0/(tau_sepbub*slope);
    for (int y=0; y<duneglobals::ny(); y++)
        for (int x=0; x<duneglobals::nx(); x++) {
            // Reducing the shear stress below separation surface (to mimic the
            //turbulence effects)
            h_delta= 1 - delta*(hSepBub(x,y)-h(x,y));
            if(h_delta < 0):
                 h_delta= 0.0;
            elif (h_delta > 1.0):
                 h_delta= 1.0;


            // Acotation
            if(TauP(x,y)[0] > 10):
                 TauP(x,y)[0] = 10;
            elif(TauP(x,y)[0] < -10):
                 TauP(x,y)[0] = -10;

            if(TauP(x,y)[1] > 10):
                 TauP(x,y)[1] = 10;
            elif(TauP(x,y)[1] < -10):
                 TauP(x,y)[1] = -10;

            // Shear stress calculation
            tauZero *= dTau0 * fabs( 1 + 0.5*TauP(x,y)[0])*( 1 + 0.5*TauP(x,y)[0]) * h_delta;
            tauOne *= bTauY * dTau0 * fabs( 1 + 0.5*TauP(x,y)[0])*TauP(x,y)[1] * h_delta;

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
