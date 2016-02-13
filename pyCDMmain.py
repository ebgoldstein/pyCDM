#Main pyCDM script. runs the model
from CDMparams import *

#Make/get the initial conditions (should make these 3D for time)
#sand surface
Topo=TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny)
#vegetation surface
Veg=VegDomain(nx,ny)


#for initial time: timestep: final time
  #DuneEvolution
  #save each increment
