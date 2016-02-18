"""
Get it ('runCDM' not 'RUN DMC')

"""

import numpy as np
from CDMparams import *
import CDMfunctions as cdmfns
import DuneEvolution as devol


#Main pyCDM script. runs the model

#MAKE THE INITIAL CONDIITONS!
# (should make these 3D for time)
#sand domain
Topo=cdmfns.TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny)
#vegetation domain
Veg=cdmfns.VegDomain(nx,ny)
#shear stress domain
Tau=cdmfns.TauDomain(nx,ny)


#for initial time: timestep: final time

  # move for 1 step
  #devol.step_implementation

  #save topo and veg each increment
