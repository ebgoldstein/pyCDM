
"""
    m_Sepbub    = P.getdefault("sepbub.parabolic", false);
    m_Smooth    = P.getdefault("sepbub.smooth", 6);
    m_reattach_length    = P.getdefault("bubble.length", 0);
    m_Slope     = P.getdefault("sepbub.slope", 0.2 /*0.2*/);

    m_x_periodic= P.getdefault("calc.x_periodic", 0);
    m_y_periodic= P.getdefault("calc.y_periodic", 0);


    m_sepslope = tan(M_PI*P.getdefault("sep.angle", 20.0)/180.0);

    m_pS=np.zeros((ny,nx))
    m_pMask=np.zeros((ny,nx))
"""
def Calc(h_sepbub, stall, h):

    #Auxiliar const
    #double dhdxb, hb, C, B, l; //, h_plain=1e2;

    # Detect flow separation
    GradMid=np.gradient(h)
    #grad_h_down.GradUpWind(h, grad_h);

    h_plain=np.zeros((ny,nx))
    grad_h_x=np.zeros((ny,nx))

    #matrix
    #if we are not at the back edge, AND stall at the forward cell is less than 1 AND stall at the given cell is >1
    #
    if(x < duneglobals::nx()-1 && stall(x+1,y) < 0 && stall(x,y) > 0):
        stall(x+1,y)= 0;

    (*m_pMask)(x,y) = (stall(x,y) < 0 && grad_h_down(x,y)[0] < 0 && vabs(grad_h_down(x,y)) > m_sepslope ? 1.0 : -1.0);
    grad_h_x(x,y)= grad_h(x,y)[0];

    /* Smoothing h_plain */
    for (int i=0; i<10 +0*abs(m_Smooth); i++) {
        m_pS->Smooth(grad_h_x);
        grad_h_x.Smooth(*m_pS);
    }
    for (int y=0; y< duneglobals::ny(); y++){
    	for (int x=duneglobals::nx()-1; x > -1; x--){
            if((*m_pMask)(x,y) > 0 && (*m_pMask)(x,y-1) < 0 && (*m_pMask)(x,y+1) < 0)	stall(x,y) = -1.0;
            else	stall(x,y) = (*m_pMask)(x,y);
		}
	}
    for (int y=0; y< duneglobals::ny(); y++){
    	for (int x=duneglobals::nx()-1; x > -1; x--)
            if(grad_h_x(x,y)>= 0 && grad_h_x((!x ? duneglobals::nx()-2 : x-1),y)<0):
                	h_plain(x,y) = h(x,y);
    		else:
                 h_plain(x,y)= h_plain((x==duneglobals::nx()-1?1:x+1),y);
            if(h_plain(x,y) < 0):
                	h_plain(x,y)= 0;
            if(h_plain(x,y)>=h(x,y)):
                	h_plain(x,y)= h(x,y);

        if(m_x_periodic):
            h_plain(duneglobals::nx()-1,y)= h_plain(0,y);


    # calculation of the separation bubble at the brink,
    # sliced model, (2d), y is considered as a parameter!for (int i=0; i<abs(m_Smooth); i++) {

    for (int y=0; y<Ny; y++)
    {
        x= 0,
        x_sep= 0,
        x_brink= 0,
        x_slope= 0,
        x_next= 0,
        x_index;
        while(x < Nx){
            // Assume no separation
            h_sepbub(x,y) = h(x,y);
            # Looking for separation
            if(x<4):
                if(m_x_periodic):
                    	x_brink = Nx-1+(x-1);
            else:
                x_brink = x-1;
            if(x==Nx-1):
                x_next=(m_x_periodic ? 0 : x);
            else:
                x_next= x+1;
            if((x>3 || m_x_periodic) && stall(x,y)>0 && stall(x_brink,y)<0 && h(x,y) > h(x_next,y)):
            
                // x:   first slip face point
                // x_brink=x-1: crest/brink point (not used due to high
                //      fluctuations while moving along the grid.
                // Xb=x-2: start of separation bubble (x0)
                if(x<4){
                    if(m_x_periodic){
                        x_sep = (x>1 ? x-2 : Nx-1+(x-2));
                        x_slope = Nx-1+(x-4);
                    }
                }else{
                    x_sep = x-2;
                    x_slope = x-4;
                }
                hb= h(x_sep,y);
                dhdxb= 0.5*(hb-h(x_slope,y))/dx;
                dhdxb= (dhdxb > 0.64 ? 0.64 : dhdxb);
                hb-= h_plain(x_sep,y);
                if(!m_reattach_length){
                    double a= dhdxb/m_Slope;
                    l=hb*1.5*(1.+0.25*a+0.125*a*a)/m_Slope;
                }else	l=hb*m_reattach_length;
                if(l<1.5*hb) l=1.5*hb;
                if(m_Sepbub)
                    C=-(hb+dhdxb*l)/(l*l);
                else{
                    B=(2*hb+dhdxb*l)/(l*l*l);
                    C=-(3*hb+2*dhdxb*l)/(l*l);
                }
                int Xb = x-2;
                while(x< Xb + (int)(l/dx)){
                    double X=(x-Xb)*dx;
                    if(x > Nx-1)
                        if(m_x_periodic)	x_index = x-Nx;
                        else break;
                        else	x_index = x;
                    h_sepbub(x_index,y)= h_plain(x_sep,y) + (m_Sepbub ? C*X*X : B*X*X*X+C*X*X)+dhdxb*X+hb;
                    if(x > Xb+2 && h_sepbub(x_index,y) < h(x_index,y)){
                        h_sepbub(x_index,y)= h(x_index,y);
                        break;
                    }
                    x++;
                }
            }	// if (sepbubble)
            x++;
            // scan after the separation bubble for the next one
        } // for x
    } // for y(!x ? duneglobals::nx()-2 : x-1)


    // Smooth separation bubble in order to have a smooth object in
    // the lateral direction.

    if (abs(m_Smooth) > 0) {
        // We only want to smooth the sep. bubble (h_sepbub - h), NOT the surface!

        for (int y=0; y<Ny; y++) {
            for (int x=0; x<Nx; x++) {
                if (h_sepbub(x,y) - h(x,y) > 1e-5) {
                    (*m_pMask)(x,y) = 1.;
                } else {
                    (*m_pMask)(x,y) = -1.;
                }
            }
        }

        for (int i=0; i<abs(m_Smooth); i++) {
            m_pS->Smooth(h_sepbub);
            h_sepbub.Smooth(*m_pS);

            for (int y=0; y<Ny; y++) {
                for (int x=0; x<Nx; x++) {
                    if ((*m_pMask)(x,y) < 0. || h_sepbub(x,y)<h(x,y)) {
                        h_sepbub(x,y) = h(x,y);
                    }
                }
            }
        }

    }


    stall = grad_h;
    return
